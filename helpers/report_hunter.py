"""
track down the latest version of all reports
"""


from os.path import join
from google.cloud import storage
from cpg_utils import to_path


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

all_cohorts = {}

for cohort in cohorts:
    all_cohorts[cohort] = {}
    ROOT = TEMPLATE.format(dataset=cohort)
    bucket_name, prefix = ROOT.removeprefix('gs://').split('/', maxsplit=1)
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = sorted(
        client.list_blobs(bucket_name, prefix=prefix), key=lambda x: x.updated
    )

    # iterate through, oldest to latest, replacing as new ones are found
    for blob in blobs:
        if not blob.name.endswith('html'):
            continue
        if 'singleton_' in blob.name:
            if '/exomes/' in blob.name:
                all_cohorts[cohort]['singleton_exome'] = join(
                    WEB_TEMPLATE, cohort, blob.name
                )
            else:
                all_cohorts[cohort]['singleton_genome'] = join(
                    WEB_TEMPLATE, cohort, blob.name
                )
        else:
            if '/exomes/' in blob.name:
                all_cohorts[cohort]['exome'] = join(WEB_TEMPLATE, cohort, blob.name)
            else:
                all_cohorts[cohort]['genome'] = join(WEB_TEMPLATE, cohort, blob.name)

HTML = """
<!DOCTYPE html>

<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Latest AIP Reports</title>
    <meta name="description" content="Latest AIP Reports">
  </head>
  <body>
"""

for cohort, categories in all_cohorts.items():
    HTML += f'\t<h2>{cohort}</h2><br>\n'
    for cat in sorted(categories):
        HTML += f'\t\t<a href="{categories[cat]}" target="_blank">{cat}</a><br>\n'
HTML += '</body>'

# write this file to a web bucket - either attached to a single dataset, or communal
output_path = join(TEMPLATE.format(dataset='common'), 'aip_index.html')
output_path = output_path.replace('main', 'test')

with to_path(output_path).open('w') as handle:
    handle.writelines(HTML)
