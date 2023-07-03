#!/usr/bin/env python3


"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""

from dataclasses import dataclass
from os.path import join
from pathlib import Path
from typing import Any

import jinja2

from cpg_utils import to_path
from cpg_utils.config import get_config
from metamist.graphql import gql, query

JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'


# pylint: disable=unsubscriptable-object


@dataclass
class Report:
    """
    generic object for storing report details
    """

    dataset: str
    address: str
    genome_or_exome: str
    subtype: str
    date: str


def get_my_projects():
    """
    queries metamist for projects I have access to,
    returns the dataset names
    """
    project_query = gql(
        """
query MyQuery {
    myProjects {
        dataset
    }
}
    """
    )
    # validate(project_query)
    response: dict[str, Any] = query(project_query)
    return {dataset['dataset'] for dataset in response['myProjects']}


def get_project_analyses(project: str) -> list[dict]:
    """
    find all the active analysis entries for this project
    Args:
        project (str): project to query for
    """

    project_query = gql(
        """
    query MyQuery($project: String!) {
        project(name: $project) {
            analyses(active: {eq: true}, type:  {eq: "web"}) {
                output
                meta
                timestampCompleted
            }
        }
    }
    """
    )
    # validate(project_query)
    response: dict[str, Any] = query(project_query, variables={'project': project})
    return response['project']['analyses']


def main():
    """
    finds all existing reports, generates an HTML file
    """

    all_cohorts = {}

    for cohort in get_my_projects():
        if 'test' in cohort:
            continue

        for analysis in get_project_analyses(cohort):
            # only look for HTML reanalysis entries
            if 'reanalysis' not in analysis['output']:
                continue

            # pull the exome/singleton flags
            exome_output = analysis['meta'].get('is_exome', False)
            singleton_output = analysis['meta'].get('is_singleton', False)

            # incorporate that into a key when gathering
            all_cohorts[f'{cohort}_{exome_output}_{singleton_output}'] = Report(
                dataset=cohort,
                address=analysis['meta']['display_url'],
                genome_or_exome='Exome' if exome_output else 'Genome',
                subtype='Singleton' if singleton_output else 'Familial',
                date=analysis['timestampCompleted'].split('T')[0],
            )

    # smoosh into a list for the report context - all reports sortable by date
    template_context = {'reports': list(all_cohorts.values())}

    # build some HTML
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(JINJA_TEMPLATE_DIR),
    )
    template = env.get_template('index.html.jinja')
    content = template.render(**template_context)

    # write to common web bucket - either attached to a single dataset, or communal
    to_path(
        join(
            get_config()['storage']['common']['test']['web'],
            'reanalysis',
            'aip_index.html',
        )
    ).write_text('\n'.join(line for line in content.split('\n') if line.strip()))


if __name__ == '__main__':
    main()
