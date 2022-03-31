#!/usr/bin/env python3


"""
Entrypoint for the interpretation pipeline process, runs the end-to-end
pipeline stages either directly or via Hail Batch(es)
 - Data extraction from PanelApp
 - Filtering and Annotation of variant data
 - Re-headering of resultant VCF

Steps are run only where the specified output does not exist
i.e. the full path to the output file is crucial, and forcing steps to
re-run currently requires the deletion of previous outputs

compound-het calculations moved to Hail, removed requirement for Slivar stage
"""


from typing import Any, Dict, Optional, Union

import json
import logging
import os
import sys

import click
from cloudpathlib import AnyPath
import hailtop.batch as hb

from cpg_utils.hail import init_batch, output_path, remote_tmpdir

from query_panelapp import main as panelapp_main
from hail_filter_and_categorise import main as category_main
from validate_classifications import main as validate_main


# static paths to write outputs
PANELAPP_JSON_OUT = output_path('panelapp_137_data.json')
HAIL_VCF_OUT = output_path('hail_categorised.vcf.bgz')
COMP_HET_JSON = output_path('hail_comp_het.json')
REHEADERED_OUT = output_path('hail_categories_reheadered.vcf.bgz')
MT_TMP = output_path('tmp_hail_table.mt', category='tmp')
RESULTS_JSON = output_path('summary_results.json')

# location of the CPG BCFTools image - to be removed with new cpg-utils package
AR_REPO = 'australia-southeast1-docker.pkg.dev/cpg-common/images'
BCFTOOLS_TAG = 'bcftools:1.10.2--h4f4756c_2'
BCFTOOLS_IMAGE = f'{AR_REPO}/{BCFTOOLS_TAG}'

# local script references
HAIL_FILTER = os.path.join(os.path.dirname(__file__), 'hail_filter_and_categorise.py')
RESULTS_SCRIPT = os.path.join(os.path.dirname(__file__), 'validate_categories.py')
DATAPROC_SETUP_SCRIPTS = [
    'gs://cpg-reference/hail_dataproc/install_common.sh',
    'gs://cpg-reference/vep/vep-GRCh38.sh',  # install & configure VEP 105
]


def read_json_dict_from_path(bucket_path: str) -> Dict[str, Any]:
    """
    take a path to a JSON file, read into an object
    :param bucket_path:
    """
    with open(AnyPath(bucket_path), encoding='utf-8') as handle:
        return json.load(handle)


def set_job_resources(job: Union[hb.batch.job.BashJob, hb.batch.job.Job]):
    """
    applied resources to the job
    :param job: apply resources to _this_ job
    """
    job.cpu(2)
    job.memory('standard')
    job.storage('20G')


def handle_hail_filtering(matrix_path: str, config: str):
    """
    hail-query backend version of the filtering implementation
    use the init query service instead of running inside dataproc

    :param matrix_path: path to annotated matrix table
    :param config: can this be used as a dict if already loaded?
    :return:
    """
    init_batch()
    category_main(
        mt_input=matrix_path,
        mt_tmp=MT_TMP,
        config_path=config,
        panelapp_path=PANELAPP_JSON_OUT,
        out_vcf=HAIL_VCF_OUT,
    )


def handle_reheader_job(
    batch: hb.Batch,
    local_vcf: str,
    config_dict: Dict[str, Any],
    prior_job: Optional[hb.batch.job.BashJob] = None,
) -> hb.batch.job.BashJob:
    """
    runs the bcftools re-header process
    :param batch:
    :param local_vcf:
    :param config_dict:
    :param prior_job:
    :return:
    """

    bcft_job = batch.new_job(name='bcftools_reheader_stage')
    set_job_resources(bcft_job)
    bcft_job.image(BCFTOOLS_IMAGE)

    if prior_job is not None:
        bcft_job.depends_on(prior_job)

    bcft_job.declare_resource_group(
        vcf={'vcf': '{root}.vcf.bgz', 'vcf.tbi': '{root}.vcf.bgz.tbi'}
    )

    # reheader the VCF using BCFtools and sed
    # replace the empty description with the full CSQ line from config
    desc = '##INFO=<ID=CSQ,Number=.,Type=String,Description="'

    # grotty string formatting to deliver the correct syntax to bcftools
    conf_csq = config_dict['variant_object'].get('csq_string').replace('|', r'\|')
    new_format = rf"Format: '{conf_csq}'"

    # only to reduce line length
    b_rh = 'bcftools reheader'

    bcft_job.command(
        'set -ex; '
        f'bcftools view -h {local_vcf} | sed \'s/'
        f'{desc}">/{desc}{new_format}">/\' > new_header; '
        f'{b_rh} -h new_header --threads 4 -o {bcft_job.vcf["vcf"]} {local_vcf}; '
        f'tabix {bcft_job.vcf["vcf"]}; '
    )
    return bcft_job


@click.command()
@click.option(
    '--matrix', 'matrix_path', help='variant matrix table to analyse', required=True
)
@click.option(
    '--config_json',
    help='dictionary of runtime settings',
    default='gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json',
)
@click.option(
    '--panelapp_version',
    help='panelapp current comparison with this earlier version',
    required=False,
)
@click.option(
    '--panel_genes',
    help='location of a Gene list for use in analysis',
    required=False,
)
@click.option('--pedigree', help='location of a PED file')
def main(
    matrix_path: str,
    config_json: str,
    panelapp_version: Optional[str],
    panel_genes: Optional[str],
    pedigree: str,
):
    """
    main method, which runs the full reanalysis process

    :param matrix_path: annotated input matrix table
    :param config_json:
    :param panelapp_version:
    :param panel_genes:
    :param pedigree:
    """

    logging.info('Starting the reanalysis batch')

    config_dict = read_json_dict_from_path(config_json)

    service_backend = hb.ServiceBackend(
        billing_project=os.getenv('HAIL_BILLING_PROJECT'),
        remote_tmpdir=remote_tmpdir(),
    )

    # -------------------------------- #
    # query panelapp for panel details #
    # -------------------------------- #
    # no need to launch in a separate batch, minimal dependencies
    panelapp_main(
        panel_id='137',
        out_path=PANELAPP_JSON_OUT,
        previous_version=panelapp_version,
        gene_list=panel_genes,
    )

    # ----------------------- #
    # run hail categorisation #
    # ----------------------- #
    # permit continuity if the hail job isn't required
    prior_job = None

    # only run if the output VCF doesn't already exist
    if not AnyPath(REHEADERED_OUT).exists():

        if not AnyPath(HAIL_VCF_OUT).exists():
            # do we need to run the full annotation stage?
            if not AnyPath(matrix_path.rstrip('/') + '/').exists():
                raise Exception(
                    f'Currently this process demands an annotated '
                    f'MatrixTable. The provided path "{matrix_path}" '
                    f'does not exist or is inaccessible'
                )

            handle_hail_filtering(matrix_path=matrix_path, config=config_json)

        # --------------------------------- #
        # bcftools re-headering of hail VCF #
        # --------------------------------- #
        reheader_batch = hb.Batch(
            name='run_BCFTools (re-header)',
            backend=service_backend,
            cancel_after_n_failures=1,
        )
        # copy the Hail output file into the remaining batch jobs
        hail_output_in_batch = reheader_batch.read_input_group(
            **{'vcf': HAIL_VCF_OUT, 'vcf.tbi': HAIL_VCF_OUT + '.tbi'}
        )

        # this is no longer explicitly required...
        # it was required to run slivar: geneId, consequences, and transcript
        # if we can avoid using slivar for comp-hets, this isn't required
        # when extracting the consequences in python we can use the config string
        # this would mean the VCF is mostly useless when separated from the
        # config file... retain for now at least
        bcftools_job = handle_reheader_job(
            batch=reheader_batch,
            local_vcf=hail_output_in_batch['vcf'],
            config_dict=config_dict,
            prior_job=prior_job,
        )
        reheader_batch.write_output(bcftools_job.vcf, REHEADERED_OUT)
        # run the batch, and wait, so that the result metadata updates
        reheader_batch.run(wait=True)

    # now utilise the compound-hets and categorised variant VCF to identify
    # plausibly pathogenic variants where the MOI is viable compared to the
    # PanelApp expectation
    validate_main(
        class_vcf=REHEADERED_OUT,
        comp_het=COMP_HET_JSON,
        config_path=config_dict,
        out_json=RESULTS_JSON,
        panelapp=PANELAPP_JSON_OUT,
        pedigree=pedigree,
    )


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    main()  # pylint: disable=E1120
