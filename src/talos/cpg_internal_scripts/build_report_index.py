"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""

import re
from argparse import ArgumentParser
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import jinja2
from cloudpathlib.anypath import to_anypath
from loguru import logger

from metamist.graphql import gql, query

DATE_REGEX = re.compile(r'(\d{4}-\d{2}-\d{2})')

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent.parent / 'templates'
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
                outputs
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


@dataclass
class Report:
    """
    generic object for storing report details
    """

    dataset: str
    address: str
    is_exome: bool
    is_long_read: bool
    date: str
    title: str


@lru_cache(1)
def get_my_projects() -> set[str]:
    """
    Queries metamist for projects I have access to, returns the dataset names.
    """
    response: dict[str, Any] = query(PROJECT_QUERY)
    all_projects = {dataset['dataset'] for dataset in response['myProjects']}
    logger.info(f'Running for projects: {", ".join(sorted(all_projects))}')
    return all_projects


def get_project_analyses(project: str) -> dict[tuple[bool, bool], str]:
    """
    Find all the active analysis entries for this project, subdivide the analyses by long/short read and exome/genome

    Create a dictionary indexed on a double boolean key:
    - the file is long_read
    - the file is exome

    This was chosen as in the CPG infrastructure, 'genome' is the default, so 'exome' is in the path if relevant.
    Likewise for short-read vs. "long-read". This is not a perfect solution, but it works for now.

    Also... you can use a tuple as a dictionary key and that's cool.
    """

    project_reports: dict[tuple[bool, bool], str] = {}

    all_analyses = query(REPORT_QUERY, variables={'project': project})['project']['analyses']
    for analysis in all_analyses:
        # skip the much older analysis entry formats
        if not (outputs := analysis.get('outputs', None)):
            continue
        if isinstance(outputs, str):
            continue

        output_path = outputs['path']

        long = 'long_read' in output_path
        exome = 'exome' in output_path

        project_reports[(long, exome)] = output_path

    return project_reports


def main(dataset: str = 'aip') -> None:
    """
    Finds all existing reports, generates an HTML file as an index page.
    Args:
        dataset (str): The dataset to generate the index for, defaults to 'aip' for legacy reasons.
    """

    parsed_reports = {cohort: get_project_analyses(cohort) for cohort in get_my_projects()}

    report_list: list[Report] = []

    for cohort, cohort_results in parsed_reports.items():
        for (long_read, exome), report_path in cohort_results.items():
            # general - only one of these
            if report_path:
                this_file_name = Path(report_path).name
                trimmed_path = report_path.rstrip(this_file_name).rstrip('/')

                for entry in list(map(str, to_anypath(trimmed_path).glob('*.html'))):
                    report_address = entry.replace(WEB_BASE.format(cohort), WEB_URL_BASE.format(cohort))
                    report_name = entry.split('/')[-1]
                    if report_date := DATE_REGEX.search(report_address):
                        report_list.append(
                            Report(
                                dataset=cohort,
                                address=report_address,
                                is_exome=exome,
                                is_long_read=long_read,
                                date=report_date.group(1),
                                title=report_name,
                            ),
                        )

    html_from_reports(report_list, f'{dataset}_index.html')


def html_from_reports(reports: list[Report], title: str):
    """
    Build some HTML from the collection of reports we found from Metamist.
    """

    # smoosh into a list for the report context - all reports sortable by date
    template_context = {'reports': reports}

    # build some HTML
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR), autoescape=True)
    template = env.get_template('report_index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    write_index_to = to_anypath(INDEX_HOME.format(title))
    logger.info(f'Writing {title} to {write_index_to}')
    write_index_to.write_text('\n'.join(line for line in content.split('\n') if line.strip()))


def cli_main():
    """
    Command line interface for the script.
    """
    parser = ArgumentParser(description='Generate an index page for AIP reports')
    parser.add_argument('--dataset', help='Dataset for the index page', default='aip')
    args = parser.parse_args()
    main(dataset=args.dataset)


if __name__ == '__main__':
    logger.info('Fetching all reports')
    cli_main()
