"""
Create Hail Batch jobs for joint genotyping.
"""

from cpg_workflows.resources import STANDARD
from cpg_workflows.utils import can_reuse
from cpg_utils.hail_batch import image_path, command
from cpg_utils import Path
import hailtop.batch as hb
from hailtop.batch.job import Job


def add_make_sitesonly_job(
    b: hb.Batch,
    input_vcf: hb.ResourceGroup,
    output_vcf_path: Path | None = None,
    job_attrs: dict | None = None,
    storage_gb: int | None = None,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Create sites-only VCF with only site-level annotations.
    Speeds up the analysis in the AS-VQSR modeling step.

    Returns: a Job object with a single output j.sites_only_vcf of type ResourceGroup
    """
    if output_vcf_path and can_reuse(output_vcf_path):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(output_vcf_path),
                'vcf.gz.tbi': str(output_vcf_path) + '.tbi',
            }
        )

    job_name = 'MakeSitesOnlyVcf'
    job_attrs = (job_attrs or {}) | {'tool': 'gatk MakeSitesOnlyVcf'}
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    res = STANDARD.set_resources(j, ncpu=2)
    if storage_gb:
        j.storage(f'{storage_gb}G')
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    assert isinstance(j.output_vcf, hb.ResourceGroup)
    j.command(
        command(
            f"""
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    MakeSitesOnlyVcf \\
    -I {input_vcf['vcf.gz']} \\
    -O {j.output_vcf['vcf.gz']}

    if [[ ! -e {j.output_vcf['vcf.gz.tbi']} ]]; then
        tabix -p vcf {j.output_vcf['vcf.gz']}
    fi
    """
        )
    )
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf
