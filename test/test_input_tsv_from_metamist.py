"""
Tests for talos.cpg_internal.input_tsv_from_metamist, against fixture Metamist responses in test/input/metamist.

The fetch_* functions are the only Metamist I/O, so the end-to-end tests patch those and nothing else.
"""

import copy
import csv
import json
from pathlib import Path

import pytest
import toml

from talos.cpg_internal import input_tsv_from_metamist as itm

FIXTURES = Path(__file__).parent / 'input' / 'metamist'
COHORTS = json.loads((FIXTURES / 'cohorts.json').read_text())
PEDIGREE_ALPHA = json.loads((FIXTURES / 'pedigree_alpha.json').read_text())
BASE_CONFIG = toml.loads((FIXTURES / 'talos.toml').read_text())


def _cohorts(
    *cohort_ids: str,
    access_level: str = 'test',
    sequencing_type: str = 'genome',
    response=None,
) -> list[itm.Cohort]:
    return itm.parse_cohorts(
        copy.deepcopy(COHORTS) if response is None else response,
        list(cohort_ids),
        access_level,
        sequencing_type,
    )


def test_normalise_cohort_ids():
    assert itm.normalise_cohort_ids(['COH1,COH2', 'COH1', ' COH3 ', ',']) == ['COH1', 'COH2', 'COH3']


def test_parse_cohorts():
    alpha, beta = _cohorts('COH0001', 'COH0002')

    assert (alpha.id, alpha.dataset, alpha.project) == ('COH0001', 'alpha', 'alpha-test')
    # the exome and long-read SGs are dropped
    assert [sg.id for sg in alpha.sgs] == ['SGPROBAND', 'SGFATHER', 'SGMOTHER']
    # only present terms, deduplicated and sorted
    assert alpha.sgs[0].hpo_terms == ('HP:0002779', 'HP:0004322')
    assert alpha.sgs[0].participant == 'PROBAND'

    assert beta.project_meta == {}
    assert beta.sgs[0].hpo_terms == ()


def test_parse_cohorts_dataset_from_project_name():
    # Metamist reports some test projects' dataset with the -test suffix; it must not leak into bucket paths
    response = copy.deepcopy(COHORTS)
    response[0]['project']['dataset'] = 'alpha-test'
    (alpha,) = _cohorts('COH0001', response=response)
    assert alpha.dataset == 'alpha'


def test_parse_cohorts_keeps_requested_order():
    assert [cohort.id for cohort in _cohorts('COH0002', 'COH0001')] == ['COH0002', 'COH0001']


def test_parse_cohorts_exome():
    (alpha,) = _cohorts('COH0001', sequencing_type='exome')
    assert [sg.id for sg in alpha.sgs] == ['SGPROBANDEXOME']


def test_parse_cohorts_reports_every_problem():
    response = copy.deepcopy(COHORTS)
    # a second alpha cohort, so alpha has two cohorts in this run
    second_alpha = copy.deepcopy(response[0])
    second_alpha['id'] = 'COH0003'
    # an SG from another project
    foreign = copy.deepcopy(response[1])
    foreign['id'] = 'COH0004'
    foreign['project'] = {'name': 'gamma-test', 'dataset': 'gamma', 'meta': {}}
    # nothing short-read genome
    exome_only = copy.deepcopy(response[1])
    exome_only['id'] = 'COH0005'
    exome_only['project'] = {'name': 'delta-test', 'dataset': 'delta', 'meta': {}}
    exome_only['sequencingGroups'][0]['type'] = 'exome'
    exome_only['sequencingGroups'][0]['sample']['project']['name'] = 'delta-test'
    # one participant, two SGs
    resequenced = copy.deepcopy(response[1])
    resequenced['id'] = 'COH0006'
    resequenced['project'] = {'name': 'epsilon-test', 'dataset': 'epsilon', 'meta': {}}
    resequenced['sequencingGroups'] = [copy.deepcopy(resequenced['sequencingGroups'][0]) for _ in range(2)]
    for idx, sg in enumerate(resequenced['sequencingGroups']):
        sg['id'] = f'SGRESEQ{idx}'
        sg['sample']['project']['name'] = 'epsilon-test'
    response.extend([second_alpha, foreign, exome_only, resequenced])

    with pytest.raises(itm.MetamistInputError) as error:
        _cohorts('COH0001', 'COH0003', 'COH0004', 'COH0005', 'COH0006', 'COH9999', response=response)

    message = str(error.value)
    assert "missing from Metamist, or not visible to this run: ['COH9999']" in message
    assert "COH0004: belongs to gamma-test, but has sequencing groups from ['beta-test']" in message
    assert 'COH0005: no active short-read genome sequencing groups' in message
    assert (
        "COH0006: participants with more than one sequencing group: {'SINGLETON': ['SGRESEQ0', 'SGRESEQ1']}" in message
    )
    assert "Dataset alpha has more than one cohort in this run (['COH0001', 'COH0003'])" in message


@pytest.mark.parametrize('access_level', ['standard', 'full'])
def test_parse_cohorts_test_project_at_production_access_level(access_level):
    with pytest.raises(itm.MetamistInputError, match=f'project alpha-test does not match access level {access_level}'):
        _cohorts('COH0001', access_level=access_level)


def test_parse_cohorts_main_project_at_test_access_level():
    response = copy.deepcopy(COHORTS)
    response[1]['project']['name'] = 'beta'
    response[1]['sequencingGroups'][0]['sample']['project']['name'] = 'beta'
    with pytest.raises(itm.MetamistInputError, match='project beta does not match access level test'):
        _cohorts('COH0002', response=response)
    assert _cohorts('COH0002', access_level='standard', response=response)[0].project == 'beta'


def test_build_pedigree():
    (alpha,) = _cohorts('COH0001')
    rows, missing = itm.build_pedigree(alpha, PEDIGREE_ALPHA)

    assert missing == []
    # SG IDs throughout, HPO terms in the 7th column, and nothing for the unsequenced FAM2 participant
    assert rows == [
        ['FAM1', 'SGPROBAND', 'SGFATHER', 'SGMOTHER', '2', '2', 'HP:0002779,HP:0004322'],
        ['FAM1', 'SGFATHER', '0', '0', '1', '1', ''],
        ['FAM1', 'SGMOTHER', '0', '0', '2', '1', ''],
    ]


def test_build_pedigree_parent_outside_cohort():
    response = copy.deepcopy(COHORTS)
    # mother is in the Metamist pedigree, but not this cohort
    response[0]['sequencingGroups'] = [sg for sg in response[0]['sequencingGroups'] if sg['id'] != 'SGMOTHER']
    (alpha,) = _cohorts('COH0001', response=response)
    rows, _missing = itm.build_pedigree(alpha, PEDIGREE_ALPHA)
    assert rows[0][:4] == ['FAM1', 'SGPROBAND', 'SGFATHER', '0']


def test_build_pedigree_missing_entries():
    response = copy.deepcopy(COHORTS)
    response[0]['sequencingGroups'][2]['sample']['participant'] = None
    (alpha,) = _cohorts('COH0001', response=response)
    rows, missing = itm.build_pedigree(alpha, [entry for entry in PEDIGREE_ALPHA if entry['individual_id'] != 'FATHER'])

    assert [row[1] for row in rows] == ['SGPROBAND']
    assert missing == [
        'COH0001: no pedigree entry for SGFATHER (participant FATHER)',
        'COH0001: no pedigree entry for SGMOTHER (participant None)',
    ]


def test_build_ext_ids():
    (alpha,) = _cohorts('COH0001')
    assert itm.build_ext_ids(alpha) == [('SGPROBAND', 'PROBAND'), ('SGFATHER', 'FATHER'), ('SGMOTHER', 'MOTHER')]


def test_seqr_project_guid():
    alpha, beta = _cohorts('COH0001', 'COH0002')
    assert itm.seqr_project_guid(alpha, 'genome') == 'R0001_alpha'
    assert itm.seqr_project_guid(alpha, 'exome') == 'R0002_alpha'
    assert itm.seqr_project_guid(beta, 'genome') is None

    alpha.project_meta['is_seqr'] = False
    assert itm.seqr_project_guid(alpha, 'genome') is None


def test_build_seqr_map():
    pedigree = [['FAM1', 'SGPROBAND'], ['FAM1', 'SGFATHER'], ['FAM2', 'SGOTHER']]
    assert itm.build_seqr_map(pedigree, {'FAM1': 'F000001_fam1', 'FAM9': 'F000009_fam9'}) == {
        'SGPROBAND': 'F000001_fam1',
        'SGFATHER': 'F000001_fam1',
    }


def test_build_config_with_seqr():
    config = itm.build_config(BASE_CONFIG, 'alpha', 'genome', 'R0001_alpha')

    assert 'cohorts' not in config
    assert config['dataset'] == 'alpha'
    assert config['sequencing_type'] == 'genome'
    assert config['GeneratePanelData'] == {
        'default_panel': 137,
        'forced_panels': [4059, 76],
        'forbidden_genes': ['GENE1'],
    }
    assert config['ValidateMOI'] == {'gnomad_max_af': 0.01, 'solved_cases': ['SGFATHER']}
    assert config['CreateTalosHTML'] == {
        'remove_solved_cases': True,
        'external_labels': 'gs://alpha/labels.json',
        'hyperlinks': {
            'template': 'https://seqr.populationgenomics.org.au/project/R0001_alpha/family_page/{sample}',
            'variant_template': 'https://seqr.populationgenomics.org.au/variant_search/variant/{variant}/family/{sample}',
        },
    }
    # the base config is shared between cohorts, and must not be modified
    assert BASE_CONFIG['GeneratePanelData']['forced_panels'] == [1]
    assert 'hyperlinks' not in BASE_CONFIG['CreateTalosHTML']


def test_build_config_without_seqr():
    beta = itm.build_config(BASE_CONFIG, 'beta', 'genome', None)
    # no overrides for beta, bar its own hyperlinks
    assert beta['GeneratePanelData'] == BASE_CONFIG['GeneratePanelData']
    assert beta['CreateTalosHTML']['hyperlinks'] == {'template': 'https://example.org/{sample}'}

    # external labels are per sequencing type
    assert 'external_labels' not in itm.build_config(BASE_CONFIG, 'alpha', 'exome', None)['CreateTalosHTML']

    # a dataset with no [cohorts.<dataset>] section gets the generic config
    gamma = itm.build_config(BASE_CONFIG, 'gamma', 'genome', None)
    assert gamma['GeneratePanelData'] == BASE_CONFIG['GeneratePanelData']
    assert 'hyperlinks' not in gamma['CreateTalosHTML']


@pytest.mark.parametrize(
    ('outputs', 'expected'),
    [
        ('gs://a/b.vcf.bgz', 'gs://a/b.vcf.bgz'),
        ({'path': 'gs://a/b.vcf.bgz', 'basename': 'b.vcf.bgz'}, 'gs://a/b.vcf.bgz'),
        ({'output': ['gs://a/b.vcf.bgz', 'gs://a/c.vcf.bgz']}, 'gs://a/b.vcf.bgz'),
        ({'output': []}, None),
        (None, None),
    ],
)
def test_output_path(outputs, expected):
    assert itm.output_path(outputs) == expected


def test_latest_output():
    analyses = [
        {'outputs': {'path': 'gs://a/old.vcf.bgz'}, 'timestampCompleted': '2026-01-01T00:00:00'},
        {'outputs': {'path': 'gs://a/new.vcf.bgz'}, 'timestampCompleted': '2026-03-01T00:00:00'},
        # newer, but partial-dataset
        {'outputs': {'path': 'gs://a/x_families-abc.vcf.bgz'}, 'timestampCompleted': '2026-04-01T00:00:00'},
        {'outputs': None, 'timestampCompleted': '2026-05-01T00:00:00'},
    ]
    assert itm.latest_output(analyses) == 'gs://a/new.vcf.bgz'
    assert itm.latest_output([]) is None


def test_analysis_query_meta():
    assert itm.MITO_VCF.meta('genome') == {'stage': 'GenerateMitoJointCall', 'sequencing_type': 'genome'}
    assert itm.HISTORY.meta('exome') == {'sequencing_type': 'exome'}
    assert itm.SV_ANNOTATED_VCF.meta('genome') == {}


@pytest.mark.parametrize(
    ('access_level', 'bucket'),
    [('test', 'gs://cpg-alpha-test'), ('standard', 'gs://cpg-alpha-main'), ('full', 'gs://cpg-alpha-main')],
)
def test_cohort_root(access_level, bucket):
    assert itm.cohort_root('alpha', 'COH0001', access_level, 'exome') == f'{bucket}/exome/talos/COH0001'


# ---- end to end, with only the fetch_* functions patched ----

ANALYSES = {
    ('alpha-test', itm.SMALL_VARIANT_VCF): 'gs://cpg-alpha-test/alpha.vcf.bgz',
    ('alpha-test', itm.SV_ANNOTATED_VCF): 'gs://cpg-alpha-test/alpha_sv.vcf.bgz',
    ('alpha-test', itm.MITO_VCF): 'gs://cpg-alpha-test/alpha_mito.vcf.bgz',
    ('alpha-test', itm.HISTORY): 'gs://cpg-alpha-test/history.json',
    ('beta-test', itm.SMALL_VARIANT_VCF): 'gs://cpg-beta-test/beta.vcf.bgz',
    ('beta-test', itm.STR_VCF): 'gs://cpg-beta-test/beta_str.vcf.bgz',
}


@pytest.fixture(name='fake_metamist')
def fixture_fake_metamist(monkeypatch):
    """patch the Metamist I/O. Returns the analyses lookup, and a log of analysis queries made"""
    analyses = dict(ANALYSES)
    queried: list[tuple[str, itm.AnalysisQuery, str]] = []

    def fetch_analyses(project, analysis, sequencing_type) -> list[dict]:
        queried.append((project, analysis, sequencing_type))
        if path := analyses.get((project, analysis)):
            return [{'outputs': {'path': path}, 'timestampCompleted': '2026-10-01T00:00:00'}]
        return []

    pedigrees = {
        'alpha-test': PEDIGREE_ALPHA,
        'beta-test': [
            {
                'family_id': 'FAM3',
                'individual_id': 'SINGLETON',
                'paternal_id': None,
                'maternal_id': None,
                'sex': 0,
                'affected': 2,
            },
        ],
    }
    monkeypatch.setattr(itm, 'fetch_cohorts', lambda _ids: copy.deepcopy(COHORTS))
    monkeypatch.setattr(itm, 'fetch_pedigree', lambda project: pedigrees[project])
    monkeypatch.setattr(itm, 'fetch_analyses', fetch_analyses)
    monkeypatch.setattr(itm, 'fetch_seqr_family_guids', lambda _project, _seq_type: {'FAM1': 'F000001_fam1'})
    return analyses, queried


def test_gather_cohort_inputs(fake_metamist):
    _analyses, queried = fake_metamist
    alpha, beta = itm.gather_cohort_inputs(['COH0001', 'COH0002'], 'test', 'genome', BASE_CONFIG)

    assert alpha.vcf == 'gs://cpg-alpha-test/alpha.vcf.bgz'
    assert alpha.sv_annotated == 'gs://cpg-alpha-test/alpha_sv.vcf.bgz'
    assert alpha.mito == 'gs://cpg-alpha-test/alpha_mito.vcf.bgz'
    assert alpha.str_vcf is None
    assert alpha.history == 'gs://cpg-alpha-test/history.json'
    assert alpha.seqr_map == {'SGPROBAND': 'F000001_fam1', 'SGFATHER': 'F000001_fam1', 'SGMOTHER': 'F000001_fam1'}

    assert beta.str_vcf == 'gs://cpg-beta-test/beta_str.vcf.bgz'
    assert beta.seqr_map is None
    assert beta.pedigree == [['FAM3', 'SGSINGLETON', '0', '0', '0', '2', '']]

    assert {sequencing_type for _project, _analysis, sequencing_type in queried} == {'genome'}


def test_gather_cohort_inputs_exome_skips_sv(fake_metamist):
    _analyses, queried = fake_metamist
    (alpha,) = itm.gather_cohort_inputs(['COH0001'], 'test', 'exome', BASE_CONFIG)

    assert [row[1] for row in alpha.pedigree] == ['SGPROBANDEXOME']
    assert alpha.sv_annotated is None
    assert itm.SV_ANNOTATED_VCF not in {analysis for _project, analysis, _seq_type in queried}


def test_gather_cohort_inputs_reports_every_missing_input(fake_metamist, monkeypatch):
    analyses, _queried = fake_metamist
    del analyses[('alpha-test', itm.SMALL_VARIANT_VCF)]
    del analyses[('beta-test', itm.SMALL_VARIANT_VCF)]
    without_mother = [entry for entry in PEDIGREE_ALPHA if entry['individual_id'] != 'MOTHER']
    monkeypatch.setattr(itm, 'fetch_pedigree', lambda project: without_mother if project == 'alpha-test' else [])

    with pytest.raises(itm.MetamistInputError) as error:
        itm.gather_cohort_inputs(['COH0001', 'COH0002'], 'test', 'genome', BASE_CONFIG)

    assert str(error.value).splitlines() == [
        'Missing required inputs:',
        'COH0001: no pedigree entry for SGMOTHER (participant MOTHER)',
        'COH0001: no "dataset_vcf" analysis for genome in alpha-test',
        'COH0002: no pedigree entry for SGSINGLETON (participant SINGLETON)',
        'COH0002: no "dataset_vcf" analysis for genome in beta-test',
    ]


@pytest.mark.usefixtures('fake_metamist')
def test_write_cohort_inputs(tmp_path):
    alpha, beta = itm.gather_cohort_inputs(['COH0001', 'COH0002'], 'test', 'genome', BASE_CONFIG)

    row = itm.write_cohort_inputs(alpha, f'{tmp_path}/')
    assert row == {
        'cohort': 'COH0001',
        'path': 'gs://cpg-alpha-test/alpha.vcf.bgz',
        'type': 'vcf',
        'pedigree': f'{tmp_path}/COH0001_pedigree.ped',
        'config': f'{tmp_path}/COH0001_config.toml',
        'history': 'gs://cpg-alpha-test/history.json',
        'ext_ids': f'{tmp_path}/COH0001_ext_ids.tsv',
        'seqr_map': f'{tmp_path}/COH0001_seqr_map.json',
        'mito': 'gs://cpg-alpha-test/alpha_mito.vcf.bgz',
        'sv': '',
        'sv_annotated': 'gs://cpg-alpha-test/alpha_sv.vcf.bgz',
        'str': '',
        'cpg_dataset': 'alpha',
    }
    assert (tmp_path / 'COH0001_pedigree.ped').read_text().splitlines()[0] == (
        'FAM1\tSGPROBAND\tSGFATHER\tSGMOTHER\t2\t2\tHP:0002779,HP:0004322'
    )
    assert (tmp_path / 'COH0001_ext_ids.tsv').read_text() == 'SGPROBAND\tPROBAND\nSGFATHER\tFATHER\nSGMOTHER\tMOTHER\n'
    assert json.loads((tmp_path / 'COH0001_seqr_map.json').read_text())['SGPROBAND'] == 'F000001_fam1'
    assert toml.loads((tmp_path / 'COH0001_config.toml').read_text()) == alpha.config

    # the one-row TSV reads back to the same row, empty cells included
    with (tmp_path / 'COH0001_input.tsv').open() as handle:
        assert list(csv.DictReader(handle, delimiter='\t')) == [row]

    # no seqr project, no seqr map file, and an empty cell
    beta_row = itm.write_cohort_inputs(beta, str(tmp_path))
    assert beta_row['seqr_map'] == ''
    assert not (tmp_path / 'COH0002_seqr_map.json').exists()


def test_main_writes_nothing_on_failure(fake_metamist, tmp_path, monkeypatch):
    analyses, _queried = fake_metamist
    del analyses[('beta-test', itm.SMALL_VARIANT_VCF)]
    written: list[str] = []
    monkeypatch.setattr(itm, 'write_cohort_inputs', lambda _inputs, inputs_dir: written.append(inputs_dir))

    with pytest.raises(itm.MetamistInputError, match='COH0002: no "dataset_vcf" analysis'):
        itm.main(
            ['COH0001', 'COH0002'], 'test', 'genome', str(FIXTURES / 'talos.toml'), '20261030', str(tmp_path / 'o')
        )
    assert written == []
    assert not (tmp_path / 'o').exists()


@pytest.mark.usefixtures('fake_metamist')
def test_main(tmp_path, monkeypatch):
    written: list[str] = []

    def fake_write(inputs, inputs_dir) -> dict[str, str]:
        written.append(inputs_dir)
        return {'cohort': inputs.cohort.id, 'cpg_dataset': inputs.cohort.dataset}

    monkeypatch.setattr(itm, 'write_cohort_inputs', fake_write)
    output = tmp_path / 'input.tsv'
    itm.main(['COH0001', 'COH0002'], 'test', 'genome', str(FIXTURES / 'talos.toml'), '20261030', str(output))

    assert written == [
        'gs://cpg-alpha-test/genome/talos/COH0001/analysis_20261030/inputs',
        'gs://cpg-beta-test/genome/talos/COH0002/analysis_20261030/inputs',
    ]
    lines = output.read_text().splitlines()
    assert lines[0].split('\t') == itm.TSV_COLUMNS
    assert [line.split('\t')[0] for line in lines[1:]] == ['COH0001', 'COH0002']


def test_main_rejects_bad_run_date(tmp_path):
    with pytest.raises(ValueError, match='YYYYMMDD'):
        itm.main(['COH0001'], 'test', 'genome', str(FIXTURES / 'talos.toml'), '2026-10-30', str(tmp_path / 'o'))
