#!/usr/bin/env python3


"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""


from dataclasses import dataclass
from os.path import join
from pathlib import Path

import jinja2
from google.api_core.exceptions import NotFound
from google.cloud import storage

from cpg_utils import to_path
from cpg_utils.config import get_config
from sample_metadata.apis import ProjectApi


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

        try:
            root = join(get_config()['storage'][cohort]['web'], 'reanalysis')
            bucket_name, prefix = root.removeprefix('gs://').split('/', maxsplit=1)
            client = storage.Client()
            blobs = sorted(
                client.list_blobs(bucket_name, prefix=prefix), key=lambda x: x.updated
            )
        except NotFound:
            continue
        except KeyError:
            continue

        # iterate through, oldest to latest, replacing as new ones are found
        for blob in blobs:
            if not blob.name.endswith('html'):
                continue

            web_template = get_config()['storage'][cohort]['web_url']

            report_obj = Report(
                dataset=cohort,
                address=join(web_template, blob.name),
                genome_or_exome='Exome' if '/exomes/' in blob.name else 'Genome',
                subtype='Singleton' if 'singleton_' in blob.name else 'Familial',
                date=str(blob.updated).split()[0],
            )
            all_cohorts[
                f'{cohort}_{report_obj.genome_or_exome}_{report_obj.subtype}'
            ] = report_obj

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
