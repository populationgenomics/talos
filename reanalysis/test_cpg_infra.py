#!/usr/bin/env python3
'''Simple script to test whether the CPG infrastructure and permissions are configured appropriately to permit running AIP.'''

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir
import hailtop.batch as hb

# docker run -e CPG_CONFIG_PATH=/root/aip/reanalysis/master.toml -v $(realpath our-aip):/root/aip -v $(realpath $HOME/.hail):/root/.hail/ -w /root/aip -it hail

@click.command()
@click.option(
    '--blob', help='thing to process', required=True
)
def main(
    blob: str
):
    """
    main
    """

    # print("hello-olleh")
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

    j = batch.new_job()
    j.command('echo "Hello World."')

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120


batch.populationgenomics.org.au
