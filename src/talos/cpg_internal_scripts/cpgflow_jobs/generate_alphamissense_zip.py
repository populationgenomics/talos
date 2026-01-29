from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def encode_alphamissense(output_path: Path) -> 'BashJob':
    """Encode the alphamissense TSV as an Echtvar zip."""
    batch = hail_batch.get_batch()
    tsv_input = batch.read_input(config.config_retrieve(['references', 'alphaissense_tsv']))

    job = batch.new_bash_job('Encode Alphamissense TSV')
    job.image(config.config_retrieve(['workflow', 'driver_image'])).storage('10GiB')

    # parse the AM data as a VCF
    job.command(f'python -m talos.annotation_scripts.parse_alphamissense --input {tsv_input} --output {job.out}')

    # encode that VCF as a ZIP
    job.command(f'echtvar encode {job.zip} /talos/echtvar/am_config.json {job.out}')

    batch.write_output(job.zip, output_path)

    return job
