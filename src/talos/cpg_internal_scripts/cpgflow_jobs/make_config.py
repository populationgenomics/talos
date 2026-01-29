import functools
import json

import toml
from loguru import logger

from cpg_flow import targets
from cpg_utils import Path, config
from metamist.apis import ProjectApi, WebApi
from metamist.exceptions import ForbiddenException

from talos.cpg_internal_scripts import cpg_flow_utils

PROJECT_API = ProjectApi()
WEB_API = WebApi()


@functools.lru_cache(maxsize=1)
def get_project_mapping() -> dict[str, dict]:
    """Retrieve the seqr project mapping for all datasets, indexed by project name."""

    json_response = PROJECT_API.get_seqr_projects()
    return {section['name']: section for section in json_response}


def get_seqr_project(dataset: str, seq_type: str) -> str | None:
    """Retrieve the seqr project ID for this combination of project, sequencing type, and long/short read - or None."""

    try:
        project_map = get_project_mapping()
    except ForbiddenException:
        logger.info(f'Forbidden access to seqr project mappings for {dataset}')
        return None

    # if this dataset is not in the project map, return None
    if not (section := project_map.get(dataset)):
        return None

    # if is_seqr is false, skip it
    if not section['meta'].get('is_seqr', False):
        return None
    seqr_key = f'seqr-project-{seq_type}'
    if config.config_retrieve(['workflow', 'long_read'], False):
        seqr_key = f'{seqr_key}-long-read'

    # if the seqr key for this combination is not present, return None
    if not (project_id := section['meta'].get(seqr_key, '')):
        return None
    return project_id


def get_hyperlink_section(cohort: targets.Cohort, seq_type: str, mapping_path: Path) -> dict | None:
    """
    Retrieve the hyperlink section for the cohort, including seqr project ID and mapping file.
    The Seqr project and SGID-SeqrID mapping is retrieved from Metamist
    """
    dataset = cohort.dataset.name
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in dataset:
        dataset += '-test'

    if project_id := get_seqr_project(dataset, seq_type=seq_type):
        # get the mapping of Family ID to Seqr ID
        mapping = WEB_API.get_seqr_family_guid_map(seq_type, project=dataset)

        # get the mapping of SG ID to Family ID
        sg_to_fam = cpg_flow_utils.query_for_sg_family_id_map(dataset)

        # mix them up to get SG ID: Seqr ID. Allow for missing samples in a fault tolerant way
        cpg_to_seqr_id = {}
        for sg in cohort.get_sequencing_groups():
            if sg.id in sg_to_fam and sg_to_fam[sg.id] in mapping:
                cpg_to_seqr_id[sg.id] = mapping[sg_to_fam[sg.id]]

        with mapping_path.open('w', encoding='utf-8') as file_handle:
            json.dump(cpg_to_seqr_id, file_handle, indent=2)

        project_template = f'https://seqr.populationgenomics.org.au/project/{project_id}/family_page/{{sample}}'
        variant_template = 'https://seqr.populationgenomics.org.au/variant_search/variant/{variant}/family/{sample}'
        return {
            'template': project_template,
            'variant_template': variant_template,
        }

    if section := config.config_retrieve(['cohorts', cohort.dataset.name, 'hyperlinks'], False):
        return section

    return None


def create_config(cohort: targets.Cohort, seqr_out: Path, config_out: Path):
    dataset = cohort.dataset.name
    # start off with a fresh config dictionary, including generic content
    new_config = {
        'GeneratePanelData': config.config_retrieve(['GeneratePanelData']),
        'RunHailFiltering': config.config_retrieve(['RunHailFiltering']),
        'RunHailFilteringSv': config.config_retrieve(['RunHailFilteringSv']),
        'ValidateMOI': config.config_retrieve(['ValidateMOI']),
        'HPOFlagging': config.config_retrieve(['HPOFlagging']),
        'CreateTalosHTML': {},  # populate from a separate part of config
        'splice_ai_ht': config.config_retrieve(['references', 'splice_ai_ht']),
        'dataset': dataset,
        'sequencing_type': config.config_retrieve(['workflow', 'sequencing_type']),
        'long_read': config.config_retrieve(['workflow', 'long_read'], False),
        'singletons': config.config_retrieve(['workflow', 'singletons'], False),
    }

    # pull the content relevant to this cohort + sequencing type (mandatory in CPG)
    seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
    if dataset_conf := config.config_retrieve(['cohorts', dataset], False):
        seq_type_conf = dataset_conf.get(seq_type, {})

        # forbidden genes and forced panels
        new_config['GeneratePanelData'].update(
            {
                'forbidden_genes': dataset_conf.get('forbidden', []),
                'forced_panels': dataset_conf.get('forced_panels', []),
                'blacklist': dataset_conf.get('blacklist', None),
            },
        )

        # optionally, all SG IDs to remove from analysis
        new_config['ValidateMOI']['solved_cases'] = dataset_conf.get('solved_cases', [])

        if 'external_labels' in seq_type_conf:
            new_config['CreateTalosHTML']['external_labels'] = seq_type_conf['external_labels']

    # adapt to new hyperlink config structure
    if hyperlinks := get_hyperlink_section(
        cohort=cohort,
        seq_type=seq_type,
        mapping_path=seqr_out,
    ):
        new_config['CreateTalosHTML']['hyperlinks'] = hyperlinks

    with config_out.open('w') as write_handle:
        toml.dump(new_config, write_handle)
