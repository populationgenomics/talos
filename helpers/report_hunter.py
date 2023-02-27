"""
track down the latest version of all reports
generate an index HTML page with links to all reports
"""


from dataclasses import dataclass
from os.path import join
from pathlib import Path

import jinja2
from google.cloud import storage

from cpg_utils import to_path
from cpg_utils.config import get_config


JINJA_TEMPLATE_DIR = Path(__file__).absolute().parent / 'templates'
TEMPLATE = 'gs://cpg-{dataset}-main-web/reanalysis'
WEB_TEMPLATE = 'https://main-web.populationgenomics.org.au'


cohorts = [
    'acute-care',
    'ag-hidden',
    'brain-malf',
    'broad-rgp',
    'circa',
    'epileptic-enceph',
    'heartkids',
    'hereditary-neuro',
    'ibmdx',
    'kidgen',
    'leukodystrophies',
    'mito-disease',
    'ohmr3-mendelian',
    'ohmr4-epilepsy',
    'perth-neuro',
    'ravenscroft-arch',
    'ravenscroft-rdstudy',
    'rdnow',
    'schr-neuro',
    'udn-aus',
]


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

    assert 'cohorts' in get_config(), 'Config file does not contain Cohorts'
    all_cohorts = {}

    for cohort in get_config()['cohorts'].keys():
        print(cohort)
        root = TEMPLATE.format(dataset=cohort)
        bucket_name, prefix = root.removeprefix('gs://').split('/', maxsplit=1)
        client = storage.Client()
        blobs = sorted(
            client.list_blobs(bucket_name, prefix=prefix), key=lambda x: x.updated
        )

        # iterate through, oldest to latest, replacing as new ones are found
        for blob in blobs:
            if not blob.name.endswith('html'):
                continue

            report_obj = Report(
                dataset=cohort,
                address=join(WEB_TEMPLATE, cohort, blob.name),
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

    # write this file to a web bucket - either attached to a single dataset, or communal
    output_path = join(TEMPLATE.format(dataset='common'), 'aip_index.html')
    output_path = output_path.replace('main', 'test')
    to_path(output_path).write_text(
        '\n'.join(line for line in content.split('\n') if line.strip())
    )


if __name__ == '__main__':
    main()
