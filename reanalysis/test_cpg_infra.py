#!/usr/bin/env python3
"""
    Simple script to test whether the CPG infrastructure and permissions are
    configured appropriately to permit running AIP.
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir
import hailtop.batch as hb


@click.command()
@click.option('--blob', help='thing to process', required=False)
def main(blob: str):
    """
    main
    """

    service_backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    batch = hb.Batch(
        name='AIP batch',
        backend=service_backend,
        cancel_after_n_failures=1,
        default_timeout=6000,
        default_memory='highmem',
    )

    j = batch.new_job(name='Write the file')
    j.command(f'echo "Hello World." > {j.ofile}')

    k = batch.new_job(name='Read the file')
    k.command(f'cat {j.ofile}')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120
