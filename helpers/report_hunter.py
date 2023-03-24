#!/usr/bin/env python3


"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""


from dataclasses import dataclass
from os.path import join
from pathlib import Path

import jinja2

from cpg_utils import to_path
from cpg_utils.config import get_config
from sample_metadata.apis import AnalysisApi, ProjectApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel


JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'


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


def main():
    """
    finds all existing reports, generates an HTML file
    """

    all_cohorts = {}

    for cohort in ProjectApi().get_my_projects():

        # find any previous AnalysisEntries... Update to active=False
        a_query_model = AnalysisQueryModel(
            projects=[cohort], type=AnalysisType('web'), active=True
        )
        for analysis in AnalysisApi().query_analyses(
            analysis_query_model=a_query_model
        ):
            # only look for reanalysis entries
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
                date=analysis['timestamp_completed'].split('T')[0],
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
