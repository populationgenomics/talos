"""
track down the latest version of all reports
generate an index HTML page with links to all reports

Generate a second report for the latest variant only report
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from functools import lru_cache

import jinja2
from cloudpathlib.anypath import to_anypath

from metamist.graphql import gql, query

from talos.static_values import get_logger

DATE_REGEX = re.compile(r'(\d{4}-\d{2}-\d{2})')

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'
PROJECT_QUERY = gql(
    """
    query MyQuery {
        myProjects {
            dataset
        }
    }
    """,
)
REPORT_QUERY = gql(
    """
    query MyQuery($project: String!) {
        project(name: $project) {
            analyses(active: {eq: true}, type:  {eq: "aip-report"}) {
                output
                meta
                timestampCompleted
            }
        }
    }
    """,
)

WEB_BASE = 'gs://cpg-{}-main-web'
WEB_URL_BASE = 'https://main-web.populationgenomics.org.au/{}'
INDEX_HOME = 'gs://cpg-common-test-web/reanalysis/{}'


script_logger = get_logger(logger_name=__file__)


@dataclass
class Report:
    """
    generic object for storing report details
    """

    dataset: str
    address: str
    genome_or_exome: str
    date: str
    title: str


def get_my_projects() -> set[str]:
    """
    queries metamist for projects I have access to,
    returns the dataset names
    """
    response: dict[str, Any] = query(PROJECT_QUERY)
    all_projects = {dataset['dataset'] for dataset in response['myProjects']}
    script_logger.info(f'Running for projects: {", ".join(sorted(all_projects))}')
    return all_projects


def get_project_analyses(project: str) -> list[dict]:
    """
    find all the active analysis entries for this project
    Args:
        project (str): project to query for
    """

    response: dict[str, Any] = query(REPORT_QUERY, variables={'project': project})
    return response['project']['analyses']


@lru_cache(1)
def get_latest_analyses() -> dict[str, dict[str, str]]:
    """
    find the latest analysis entries for all projects

    Returns:
        dict[str, dict[str, str]]: key is project name, value is dict of sequencing type to output path
    """

    all_cohorts: dict[str, dict[str, str]] = {}

    for cohort in get_my_projects():
        for analysis in get_project_analyses(cohort):
            output_path = analysis['output']
            if 'sequencing_type' not in analysis['meta']:
                continue
            all_cohorts.setdefault(cohort, {})[analysis['meta']['sequencing_type']] = output_path
    return all_cohorts


def run_both():
    """
    run once for all main reports, then again for the latest-only reports
    """
    get_logger(__file__).info('Fetching main reports')
    main()
    get_logger().info('Fetching latest-only reports')
    main(latest=True)


def main(latest: bool = False):
    """
    finds all existing reports, generates an HTML file
    eventually we can latch onto the meta key `type:aip_output_html`
    but that won't be populated until we run more through the pipeline

    Args:
        latest (bool): whether to create the latest-only report
    """
    all_cohorts = get_latest_analyses()
    report_list: list[Report] = []

    for cohort, cohort_results in all_cohorts.items():
        for sequencing_type, output_path in cohort_results.items():
            this_file_name = Path(output_path).name
            trimmed_path = output_path.rstrip(this_file_name).rstrip('/')

            dir_contents = list(map(str, to_anypath(trimmed_path).glob('*.html')))

            for entry in filter(lambda x: (bool('latest' in x) == latest), dir_contents):
                report_address = entry.replace(WEB_BASE.format(cohort), WEB_URL_BASE.format(cohort))
                report_name = entry.split('/')[-1]
                if report_date := DATE_REGEX.search(report_address):
                    report_list.append(
                        Report(
                            dataset=cohort,
                            address=report_address,
                            genome_or_exome=sequencing_type,
                            date=report_date.group(1),
                            title=report_name,
                        ),
                    )

    # smoosh into a list for the report context - all reports sortable by date
    template_context = {'reports': report_list}

    # build some HTML
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('report_index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    write_index_to = to_anypath(INDEX_HOME.format('latest_aip_index.html' if latest else 'aip_index.html'))
    write_index_to.write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    run_both()
