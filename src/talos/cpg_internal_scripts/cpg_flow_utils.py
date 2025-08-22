import functools

import loguru
from cpg_flow import workflow
from cpg_utils import Path, config, to_path

from metamist import graphql

LONG_READ_STRING = 'LongRead'
METAMIST_ANALYSIS_QUERY = graphql.gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)
METAMIST_FAMILY_SG_QUERY = graphql.gql(
    """
    query MyQuery($dataset: String!) {
        project(name: $dataset) {
            sequencingGroups {
                sample {
                    participant {
                        families {
                            externalId
                        }
                    }
                }
                id
            }
        }
    }
""",
)


@functools.cache
def generate_dataset_prefix(
    dataset: str,
    category: str | None = None,
    stage_name: str | None = None,
    hash_value: str | None = None,
) -> Path:
    """
    Generate a dictionary of prefixes for the current workflow and Stage.
    Needed because CPG-Flow currently lacks the granularity we need for both exome/genome and short/long read
    This is intended to generate the exact same prefix as CPG-Flow would generate, so that we continue previous work
    """

    # mandatory value in a cpg-flow config
    workflow_name = config.config_retrieve(['workflow', 'name'], 'talos')

    # generated from the included samples in the workflow
    # or passed directly if provided, e.g. to target outputs to a date-specific folder, not a callset-specific one
    hash_element = hash_value or workflow.get_workflow().output_version

    # allow subdivision by short/long read, and exome/genome
    # the current protocol here is to treat short read and genome as standard, and insert clarifying elements if needed
    exome_element = 'exome' if config.config_retrieve(['workflow', 'sequencing_type']) == 'exome' else None
    long_read_element = 'long_read' if config.config_retrieve(['workflow', 'long_read'], False) else None

    # line up all the elements into an ordered list, and then join the non-None elements
    suffix = '/'.join(
        [x for x in [long_read_element, exome_element, workflow_name, hash_element, stage_name] if isinstance(x, str)],
    )

    return to_path(config.dataset_path(suffix=suffix, dataset=dataset, category=category))


def query_for_sg_family_id_map(dataset: str) -> dict[str, str]:
    """
    Query for the mapping of Seqr IDs to CPG Family IDs for a given dataset.

    Args:
        dataset (str): project to query for
    Returns:
        dict[str, str]: mapping of Seqr IDs to CPG Family IDs
    """
    # swapping to a string we can freely modify
    loguru.logger.info(f'Querying for Seqr Family ID map in {dataset}')

    result = graphql.query(METAMIST_FAMILY_SG_QUERY, variables={'dataset': dataset})

    return {
        sg['id']: sg['sample']['participant']['families'][0]['externalId']
        for sg in result['project']['sequencingGroups']
        if sg['sample']['participant']['families']
    }


@functools.cache
def query_for_latest_analysis(
    dataset: str,
    analysis_type: str,
    sequencing_type: str = 'all',
    long_read: bool = False,
) -> str | None:
    """
    Query for the latest analysis object of a given type in the requested project.

    Analysis entries for Talos all have unique types, so we can use this generic query method

    Args:
        dataset (str):         project to query for
        analysis_type (str):   analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
                               seqr_loader used 'custom': using a config entry we can decide which type to use
        sequencing_type (str): optional, if set, only return entries with meta.sequencing_type == this
        long_read (bool):      if True, will skip over any entries that are not LongRead (SNPsIndels/SV)
    Returns:
        str, the path to the latest object for the given type, or log a warning and return None
    """

    # swapping to a string we can freely modify
    query_dataset = dataset
    if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    loguru.logger.info(f'Querying for {analysis_type} in {query_dataset}')

    result = graphql.query(METAMIST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': analysis_type})

    # get all the relevant entries, and bin by date
    analysis_by_date = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and (sequencing_type in {'all', analysis['meta'].get('sequencing_type')}):
            # skip over the partial-cohort AnnotateDataset objects
            if '_families-' in analysis['output']:
                loguru.logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It is a partial-cohort AnnotateDataset object',
                )
                continue

            # manually implementing an XOR check - long read (bool) and LongRead in output must match
            if long_read != (LONG_READ_STRING in analysis['output']):
                loguru.logger.debug(
                    f'Skipping analysis {analysis["output"]} for dataset {query_dataset}. '
                    f'It does not match query parameter long_read={long_read}',
                )
                continue

            analysis_by_date[analysis['timestampCompleted']] = analysis['output']

    if not analysis_by_date:
        loguru.logger.warning(f'No Analysis Entries found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return analysis_by_date[sorted(analysis_by_date)[-1]]
