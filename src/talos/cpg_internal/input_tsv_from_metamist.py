#!/usr/bin/env python3

"""
Build the Talos input TSV for a set of CPG Metamist cohorts. This is the first process of cpg_main.nf.

Each cohort becomes one row, with the usual input TSV columns plus `cpg_dataset` (where its outputs go) and
`sv_annotated` (GATK-SV's annotated SV VCF, which skips Talos's own SV annotation). The rules:

- a cohort's sequencing groups must all be in the cohort's own project, and that project must match the
  access level (`<dataset>-test` at `test`, the main project otherwise)
- at most one cohort per dataset per run, as each dataset's outputs are laid out by cohort, and its history
  is the dataset's latest results
- only active, short-read SGs of the requested sequencing type are used. Sample IDs are SG IDs throughout
- the small-variant VCF, and a pedigree entry for every SG, are required. All missing inputs, across all
  cohorts, are reported together in one error. SV (genome only), mito, STR, history and the seqr map are
  optional, and left blank when absent

The generated pedigree, external ID map, seqr map and merged config are written to the cohort's dated
analysis directory, `gs://cpg-<dataset>-{main|test}/<sequencing_type>/talos/<cohort>/analysis_<run_date>/inputs/`,
alongside a one-row copy of the cohort's TSV. The TSV holds full URIs to all of them.

Only the `_query` and `fetch_*` functions talk to Metamist, so everything else can be tested against fixture
responses. Those import `metamist` lazily, as it's in the `cpg` extra, which only the container installs.

    python -m talos.cpg_internal.input_tsv_from_metamist --cohort_ids COH0001 COH0002 --access_level test \\
        --sequencing_type genome --talos_config gs://.../talos.toml --run_date 20261030 --output input.tsv
"""

import copy
import csv
import json
import re
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from typing import Any

import toml
from cloudpathlib.anypath import to_anypath
from loguru import logger

ACCESS_LEVELS = ('test', 'standard', 'full')
SEQUENCING_TYPES = ('genome', 'exome')
TECHNOLOGY = 'short-read'

TSV_COLUMNS = [
    'cohort',
    'path',
    'type',
    'pedigree',
    'config',
    'history',
    'ext_ids',
    'seqr_map',
    'mito',
    'sv',
    'sv_annotated',
    'str',
    'cpg_dataset',
]

HPO_KEY = 'HPO Terms (present)'
HPO_RE = re.compile(r'HP:\d{7}')

SEQR_URL = 'https://seqr.populationgenomics.org.au'


@dataclass(frozen=True)
class AnalysisQuery:
    """One kind of Metamist analysis to look up for a cohort's dataset."""

    type: str
    stage: str | None = None
    # most analyses carry meta.sequencing_type, but GATK-SV's annotated VCF exists only for genomes, and doesn't
    by_sequencing_type: bool = True

    def meta(self, sequencing_type: str) -> dict[str, str]:
        meta = {'stage': self.stage} if self.stage else {}
        if self.by_sequencing_type:
            meta['sequencing_type'] = sequencing_type
        return meta

    def __str__(self) -> str:
        return f'"{self.type}" analysis' + (f' from {self.stage}' if self.stage else '')


SMALL_VARIANT_VCF = AnalysisQuery(type='dataset_vcf')

CNV_ANNOTATED_VCF = AnalysisQuery(type='single_dataset_cnv_annotated', by_sequencing_type=False)
SV_ANNOTATED_VCF = AnalysisQuery(type='single_dataset_sv_annotated', by_sequencing_type=False)
MITO_VCF = AnalysisQuery(type='vcf', stage='GenerateMitoJointCall')
STR_VCF = AnalysisQuery(type='vcf', stage='MakeStripyJointCall')
HISTORY = AnalysisQuery(type='aip-results')

COHORTS_QUERY = """
query TalosCohorts($ids: [String!]!) {
    cohorts(id: {in_: $ids}) {
        id
        name
        project {
            name
            meta
        }
        sequencingGroups(activeOnly: {eq: true}) {
            id
            type
            technology
            sample {
                project {
                    name
                }
                participant {
                    externalId
                    phenotypes
                }
            }
        }
    }
}
"""

# the whole project's pedigree, keyed on participant external IDs, with family external IDs
PEDIGREE_QUERY = """
query TalosPedigree($project: String!) {
    project(name: $project) {
        pedigree(replaceWithFamilyExternalIds: true)
    }
}
"""

# DatasetStage outputs are linked to SGs, not cohorts, so they're found on the project
ANALYSES_QUERY = """
query TalosAnalyses($project: String!, $type: String!, $meta: JSON) {
    project(name: $project) {
        analyses(type: {eq: $type}, status: {eq: COMPLETED}, active: {eq: true}, meta: $meta) {
            outputs
            timestampCompleted
        }
    }
}
"""


class MetamistInputError(Exception):
    """Cohorts, or their inputs, that Talos can't run on."""


@dataclass(frozen=True)
class SequencingGroup:
    id: str
    # the participant's external ID, which the Metamist pedigree is keyed on. None if the sample has no participant
    participant: str | None
    hpo_terms: tuple[str, ...]


@dataclass(frozen=True)
class Cohort:
    id: str
    dataset: str
    # the Metamist project: the dataset name, plus `-test` at the test access level
    project: str
    project_meta: dict[str, Any]
    sgs: tuple[SequencingGroup, ...]


@dataclass
class CohortInputs:
    """Everything generated for, or found for, one cohort, ready to be written out."""

    cohort: Cohort
    pedigree: list[list[str]]
    ext_ids: list[tuple[str, str]]
    seqr_map: dict[str, str] | None
    config: dict[str, Any]
    vcf: str
    sv_annotated: str | None
    mito: str | None
    str_vcf: str | None
    history: str | None


# ---- Metamist I/O ----


def _query(query_string: str, variables: dict[str, Any]) -> dict[str, Any]:
    from metamist.graphql import gql, query  # noqa: PLC0415 - the cpg extra is only installed in the container

    return query(gql(query_string), variables=variables)


def fetch_cohorts(cohort_ids: list[str]) -> list[dict[str, Any]]:
    return _query(COHORTS_QUERY, {'ids': cohort_ids})['cohorts']


def fetch_pedigree(project: str) -> list[dict[str, Any]]:
    return _query(PEDIGREE_QUERY, {'project': project})['project']['pedigree']


def fetch_analyses(project: str, analysis: AnalysisQuery, sequencing_type: str) -> list[dict[str, Any]]:
    variables = {'project': project, 'type': analysis.type, 'meta': analysis.meta(sequencing_type)}
    return _query(ANALYSES_QUERY, variables)['project']['analyses']


def fetch_seqr_family_guids(project: str, sequencing_type: str) -> dict[str, str]:
    """{family external ID: seqr family GUID}, from seqr via Metamist."""
    from metamist.apis import WebApi  # noqa: PLC0415 - the cpg extra is only installed in the container

    return WebApi().get_seqr_family_guid_map(sequencing_type, project=project)


# ---- transforms ----


def normalise_cohort_ids(values: list[str]) -> list[str]:
    """Accept cohort IDs as separate values, comma-delimited, or both. Order is kept, and repeats dropped."""
    ids = [part.strip() for value in values for part in value.split(',')]
    return list(dict.fromkeys(cohort_id for cohort_id in ids if cohort_id))


def parse_hpo_terms(phenotypes: dict[str, Any] | None) -> tuple[str, ...]:
    return tuple(sorted(set(HPO_RE.findall((phenotypes or {}).get(HPO_KEY, '')))))


def parse_cohorts(
    response: list[dict[str, Any]],
    cohort_ids: list[str],
    access_level: str,
    sequencing_type: str,
) -> list[Cohort]:
    """
    Turn the cohorts query response into Cohorts, in the order requested. Every unusable cohort is reported in
    one error, not just the first.
    """
    problems: list[str] = []
    by_id = {cohort['id']: cohort for cohort in response}

    if missing := [cohort_id for cohort_id in cohort_ids if cohort_id not in by_id]:
        problems.append(f'Cohorts missing from Metamist, or not visible to this run: {missing}')

    cohorts: list[Cohort] = []
    for cohort_id in cohort_ids:
        if (raw := by_id.get(cohort_id)) is None:
            continue

        project = raw['project']['name']
        # a test run must never publish to a -main bucket, nor a production run to a -test one
        if project.endswith('-test') != (access_level == 'test'):
            problems.append(f'{cohort_id}: project {project} does not match access level {access_level}')
            continue

        if foreign := sorted({sg['sample']['project']['name'] for sg in raw['sequencingGroups']} - {project}):
            problems.append(f'{cohort_id}: belongs to {project}, but has sequencing groups from {foreign}')
            continue

        sgs: list[SequencingGroup] = []
        for sg in raw['sequencingGroups']:
            if sg['type'] != sequencing_type or sg['technology'] != TECHNOLOGY:
                logger.warning(f'{cohort_id}: skipping {sg["id"]} ({sg["technology"]} {sg["type"]})')
                continue
            participant = sg['sample']['participant']
            sgs.append(
                SequencingGroup(
                    id=sg['id'],
                    participant=participant['externalId'] if participant else None,
                    hpo_terms=parse_hpo_terms(participant['phenotypes'] if participant else None),
                ),
            )

        if not sgs:
            problems.append(f'{cohort_id}: no active {TECHNOLOGY} {sequencing_type} sequencing groups')
            continue

        sgs_by_participant: dict[str | None, list[str]] = defaultdict(list)
        for sg in sgs:
            if sg.participant:
                sgs_by_participant[sg.participant].append(sg.id)
        if repeated := {participant: ids for participant, ids in sgs_by_participant.items() if len(ids) > 1}:
            problems.append(f'{cohort_id}: participants with more than one sequencing group: {repeated}')
            continue

        cohorts.append(
            Cohort(
                id=cohort_id,
                # from the project name, not Metamist's `dataset` field, which carries `-test` for some test projects
                dataset=project.removesuffix('-test'),
                project=project,
                project_meta=raw['project']['meta'] or {},
                sgs=tuple(sgs),
            ),
        )

    cohorts_by_dataset: dict[str, list[str]] = defaultdict(list)
    for cohort in cohorts:
        cohorts_by_dataset[cohort.dataset].append(cohort.id)
    problems.extend(
        f'Dataset {dataset} has more than one cohort in this run ({ids}), run them separately'
        for dataset, ids in cohorts_by_dataset.items()
        if len(ids) > 1
    )

    if problems:
        raise MetamistInputError('\n'.join(problems))
    return cohorts


def build_pedigree(cohort: Cohort, pedigree: list[dict[str, Any]]) -> tuple[list[list[str]], list[str]]:
    """
    A 7-column pedigree for just this cohort's SGs, from the project's Metamist pedigree: family, SG, father,
    mother, sex, affected, and comma-delimited HPO terms. Parents outside the cohort become 0.

    Returns:
        the pedigree rows, and a description of each SG with no pedigree entry
    """
    by_participant = {entry['individual_id']: entry for entry in pedigree}
    sg_for_participant = {sg.participant: sg.id for sg in cohort.sgs if sg.participant}

    rows: list[list[str]] = []
    missing: list[str] = []
    for sg in cohort.sgs:
        if sg.participant is None or (entry := by_participant.get(sg.participant)) is None:
            missing.append(f'{cohort.id}: no pedigree entry for {sg.id} (participant {sg.participant})')
            continue
        rows.append(
            [
                entry['family_id'],
                sg.id,
                sg_for_participant.get(entry['paternal_id'], '0'),
                sg_for_participant.get(entry['maternal_id'], '0'),
                str(entry['sex'] or 0),
                str(entry['affected'] or 0),
                ','.join(sg.hpo_terms),
            ],
        )
    return rows, missing


def build_ext_ids(cohort: Cohort) -> list[tuple[str, str]]:
    """SG ID -> participant external ID, used to label the report."""
    return [(sg.id, sg.participant) for sg in cohort.sgs if sg.participant]


def seqr_project_guid(cohort: Cohort, sequencing_type: str) -> str | None:
    """The cohort's seqr project for this sequencing type, if its Metamist project is synced to seqr."""
    if not cohort.project_meta.get('is_seqr', False):
        return None
    return cohort.project_meta.get(f'seqr-project-{sequencing_type}') or None


def build_seqr_map(pedigree: list[list[str]], family_guids: dict[str, str]) -> dict[str, str]:
    """SG ID -> seqr family GUID, joined on family external ID. SGs whose family isn't in seqr are left out."""
    return {row[1]: family_guids[row[0]] for row in pedigree if row[0] in family_guids}


def build_config(
    base: dict[str, Any],
    dataset: str,
    sequencing_type: str,
    seqr_project: str | None,
) -> dict[str, Any]:
    """
    One cohort's Talos config: every section of the CPG talos.toml except `cohorts`, with that dataset's
    `[cohorts.<dataset>]` overrides applied, and seqr hyperlinks if the dataset is in seqr.
    """
    config = copy.deepcopy({key: value for key, value in base.items() if key != 'cohorts'})
    config['dataset'] = dataset
    config['sequencing_type'] = sequencing_type

    overrides = base.get('cohorts', {}).get(dataset, {})

    panel_section = config.setdefault('GeneratePanelData', {})
    for key, dest in (('forbidden', 'forbidden_genes'), ('forced_panels', 'forced_panels'), ('blacklist', 'blacklist')):
        if key in overrides:
            panel_section[dest] = overrides[key]

    if 'solved_cases' in overrides:
        config.setdefault('ValidateMOI', {})['solved_cases'] = overrides['solved_cases']

    html_section = config.setdefault('CreateTalosHTML', {})
    if 'external_labels' in (seq_type_overrides := overrides.get(sequencing_type, {})):
        html_section['external_labels'] = seq_type_overrides['external_labels']

    if seqr_project:
        html_section['hyperlinks'] = {
            'template': f'{SEQR_URL}/project/{seqr_project}/family_page/{{sample}}',
            'variant_template': f'{SEQR_URL}/variant_search/variant/{{variant}}/family/{{sample}}',
        }
    elif 'hyperlinks' in overrides:
        html_section['hyperlinks'] = overrides['hyperlinks']

    return config


def output_path(outputs: Any) -> str | None:
    """The file an analysis points at. `outputs` is a path, a {'path': ...} dict, or an {'output': [...]} dict."""
    if isinstance(outputs, str):
        return outputs
    if isinstance(outputs, dict):
        if isinstance(path := outputs.get('path'), str):
            return path
        if isinstance(paths := outputs.get('output'), list) and paths:
            return paths[0]
    return None


def latest_output(analyses: list[dict[str, Any]]) -> str | None:
    """The newest analysis's output, by completion time, skipping partial-dataset (`_families-`) outputs."""
    by_timestamp: dict[str, str] = {}
    for analysis in analyses:
        if (path := output_path(analysis['outputs'])) is None:
            logger.debug(f'Skipping analysis with unexpected outputs: {analysis["outputs"]}')
            continue
        if '_families-' in path:
            logger.debug(f'Skipping partial-dataset output {path}')
            continue
        by_timestamp[analysis['timestampCompleted']] = path
    # ISO timestamps, so they sort as strings
    return by_timestamp[max(by_timestamp)] if by_timestamp else None


def cohort_root(dataset: str, cohort_id: str, access_level: str, sequencing_type: str) -> str:
    """Where everything for this cohort is published. cpg_main.nf must agree."""
    bucket = f'gs://cpg-{dataset}-{"test" if access_level == "test" else "main"}'
    return f'{bucket}/{sequencing_type}/talos/{cohort_id}'


# ---- gather and write ----


def gather_cohort_inputs(
    cohort_ids: list[str],
    access_level: str,
    sequencing_type: str,
    base_config: dict[str, Any],
) -> list[CohortInputs]:
    """Query Metamist for everything the cohorts need. Every missing required input is reported in one error."""
    cohorts = parse_cohorts(fetch_cohorts(cohort_ids), cohort_ids, access_level, sequencing_type)

    problems: list[str] = []
    gathered: list[CohortInputs] = []
    for cohort in cohorts:
        pedigree, missing = build_pedigree(cohort, fetch_pedigree(cohort.project))
        problems.extend(missing)

        def find(analysis: AnalysisQuery) -> str | None:
            return latest_output(fetch_analyses(cohort.project, analysis, sequencing_type))  # noqa: B023 - called at once

        if not (vcf := find(SMALL_VARIANT_VCF)):
            problems.append(f'{cohort.id}: no {SMALL_VARIANT_VCF} for {sequencing_type} in {cohort.project}')

        seqr_project = seqr_project_guid(cohort, sequencing_type)
        seqr_map = (
            build_seqr_map(pedigree, fetch_seqr_family_guids(cohort.project, sequencing_type)) if seqr_project else None
        )

        sv_vcf = None
        if sequencing_type == 'genome':
            sv_vcf = find(SV_ANNOTATED_VCF)
        elif sequencing_type == 'exome':
            sv_vcf = find(CNV_ANNOTATED_VCF)

        gathered.append(
            CohortInputs(
                cohort=cohort,
                pedigree=pedigree,
                ext_ids=build_ext_ids(cohort),
                seqr_map=seqr_map,
                config=build_config(base_config, cohort.dataset, sequencing_type, seqr_project),
                vcf=vcf or '',
                sv_annotated=sv_vcf,
                mito=find(MITO_VCF),
                str_vcf=find(STR_VCF),
                history=find(HISTORY),
            ),
        )

    if problems:
        raise MetamistInputError('Missing required inputs:\n' + '\n'.join(problems))
    return gathered


def write_cohort_inputs(inputs: CohortInputs, inputs_dir: str) -> dict[str, str]:
    """Write one cohort's generated files, and its one-row TSV, to inputs_dir. Returns the TSV row."""
    cohort_id = inputs.cohort.id
    inputs_dir = inputs_dir.rstrip('/')

    def write(name: str, content: str) -> str:
        path = f'{inputs_dir}/{cohort_id}_{name}'
        to_anypath(path).write_text(content, encoding='utf-8')
        return path

    pedigree = write('pedigree.ped', ''.join('\t'.join(row) + '\n' for row in inputs.pedigree))
    ext_ids = write('ext_ids.tsv', ''.join(f'{sg_id}\t{participant}\n' for sg_id, participant in inputs.ext_ids))
    config = write('config.toml', toml.dumps(inputs.config))
    seqr_map = write('seqr_map.json', json.dumps(inputs.seqr_map, indent=2)) if inputs.seqr_map is not None else ''

    row = {
        'cohort': cohort_id,
        'path': inputs.vcf,
        'type': 'vcf',
        'pedigree': pedigree,
        'config': config,
        'history': inputs.history or '',
        'ext_ids': ext_ids,
        'seqr_map': seqr_map,
        'mito': inputs.mito or '',
        'sv': '',
        'sv_annotated': inputs.sv_annotated or '',
        'str': inputs.str_vcf or '',
        'cpg_dataset': inputs.cohort.dataset,
    }
    with to_anypath(f'{inputs_dir}/{cohort_id}_input.tsv').open('w', encoding='utf-8') as handle:
        write_tsv([row], handle)
    return row


def write_tsv(rows: list[dict[str, str]], handle) -> None:
    writer = csv.DictWriter(handle, TSV_COLUMNS, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)


def cli_main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        '--cohort_ids',
        nargs='+',
        required=True,
        help='Metamist cohort IDs, as separate values and/or comma-delimited',
    )
    parser.add_argument('--access_level', required=True, choices=ACCESS_LEVELS)
    parser.add_argument('--sequencing_type', required=True, choices=SEQUENCING_TYPES)
    parser.add_argument('--talos_config', required=True, help='CPG talos.toml, with [cohorts.<dataset>] sections')
    parser.add_argument('--run_date', required=True, help="YYYYMMDD, the date in this run's analysis directories")
    parser.add_argument('--output', required=True, help='Where to write the combined input TSV')
    args = parser.parse_args()
    main(
        cohort_ids=normalise_cohort_ids(args.cohort_ids),
        access_level=args.access_level,
        sequencing_type=args.sequencing_type,
        talos_config=args.talos_config,
        run_date=args.run_date,
        output=args.output,
    )


def main(
    cohort_ids: list[str],
    access_level: str,
    sequencing_type: str,
    talos_config: str,
    run_date: str,
    output: str,
):
    """
    Build and write the input TSV for these cohorts, and each cohort's generated inputs.

    Args:
        cohort_ids (list[str]): Metamist cohort IDs
        access_level (str): test, standard or full
        sequencing_type (str): genome or exome
        talos_config (str): path to the CPG talos.toml
        run_date (str): YYYYMMDD, matching the analysis directories cpg_main.nf publishes to
        output (str): where to write the combined TSV
    """
    if not re.fullmatch(r'\d{8}', run_date):
        raise ValueError(f'run_date must be YYYYMMDD, got {run_date!r}')
    if not cohort_ids:
        raise ValueError('No cohort IDs provided')

    base_config = toml.loads(to_anypath(talos_config).read_text(encoding='utf-8'))

    # gather everything before writing anything, so a failed run leaves nothing behind
    gathered = gather_cohort_inputs(cohort_ids, access_level, sequencing_type, base_config)

    rows = []
    for inputs in gathered:
        root = cohort_root(inputs.cohort.dataset, inputs.cohort.id, access_level, sequencing_type)
        rows.append(write_cohort_inputs(inputs, f'{root}/analysis_{run_date}/inputs'))
        logger.info(f'{inputs.cohort.id} ({inputs.cohort.dataset}): {len(inputs.pedigree)} SGs, inputs in {root}')

    with open(output, 'w', encoding='utf-8') as handle:
        write_tsv(rows, handle)


if __name__ == '__main__':
    cli_main()
