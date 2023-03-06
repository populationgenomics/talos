#!/usr/bin/env python3


"""
Entrypoint for clinvar summary generation
"""

from datetime import datetime

import click
from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job
from cpg_workflows.batch import get_batch
from reanalysis import summarise_clinvar_entries


@click.command
@click.option('--ht_out', help='Path to write the Hail table to')
@click.option('--date', help='date cut-off, optional', default=None)
def main(ht_out: str, date: str | None = None):
    """
    run the clinvar summary, output to defined path
    Args:
        ht_out ():
        date ():
    """

    clinvar_table_path = to_path(ht_out)
    clinvar_folder = clinvar_table_path.parent

    # create a bash job to copy data from remote
    bash_job = get_batch().new_bash_job(name='copy clinvar files to local')
    bash_job.image(get_config()['workflow']['driver_image'])

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(
        (
            f'wget -q {directory}{sub_file} -O {bash_job.subs} && '
            f'wget -q {directory}{var_file} -O {bash_job.vars}'
        )
    )

    # construct a date String to identify the origina date for these files
    today = datetime.now().strftime('%y-%m-%d')

    # write output files date-specific
    get_batch().write_output(bash_job.subs, str(clinvar_folder / f'{today}_{sub_file}'))
    get_batch().write_output(bash_job.vars, str(clinvar_folder / f'{today}_{var_file}'))

    # create a job to run the summary
    summarise = get_batch().new_job(name='summarise clinvar')
    summarise.cpu(2).image(get_config()['workflow']['driver_image']).storage('20G')
    authenticate_cloud_credentials_in_job(summarise)
    command_options = f'-s {bash_job.subs} -v {bash_job.vars} -o {ht_out}'
    if date:
        command_options += f' -d {date}'
    summarise.command(f'python3 {summarise_clinvar_entries.__file__} {command_options}')
    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
