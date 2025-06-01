import functools

import loguru

from cpg_utils import config
from metamist import graphql
from cpg_flow import targets


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


@functools.cache
def query_for_latest_analysis(
    dataset: str,
    analysis_type: str,
    sequencing_type: str = 'all',
) -> str | None:
    """
    Query for the latest analysis object of a given type in the requested project
    Analysis entries for Talos all have unique types, so we can use this generic query method

    Args:
        dataset (str):         project to query for
        analysis_type (str):   analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
                               seqr_loader used 'custom': using a config entry we can decide which type to use
        sequencing_type (str): optional, if set, only return entries with meta.sequencing_type == this
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
            analysis_by_date[analysis['timestampCompleted']] = analysis['output']

    if not analysis_by_date:
        loguru.logger.warning(f'No Analysis Entries found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return analysis_by_date[sorted(analysis_by_date)[-1]]
