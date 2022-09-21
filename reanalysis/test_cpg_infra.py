'''Simple script to test whether the CPG infrastructure and permissions are configured appropriately to permit running AIP.'''

import click

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

    # if not AnyPath(blob).exists():
    #     raise Exception(
    #         f'The provided path "{blob}" does not exist or is inaccessible'
    #     )

    # service_backend = hb.ServiceBackend(
    #     billing_project=get_config()['hail']['billing_project'],
    #     remote_tmpdir=remote_tmpdir(),
    # )
    # batch = hb.Batch(
    #     name='AIP batch',
    #     backend=service_backend,
    #     cancel_after_n_failures=1,
    #     default_timeout=6000,
    #     default_memory='highmem',
    # )

    print(f"File to process: {blob}")

    # batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=E1120


