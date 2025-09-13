from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def make_bcftools_mito_jobs(
    cohort_id: str,
    mito_vcf: str,
    output: Path,
    job_attrs: dict,
) -> 'BashJob':
    """"""

    batch_instance = hail_batch.get_batch('Talos')

    # localise the mito VCF, this is on chrM
    local_mito_vcf = batch_instance.read_input(mito_vcf)

    # get the GFF3 file required to generate consequences
    gff3_file = batch_instance.read_input(config.config_retrieve(['references', 'ensembl_gff3']))

    # get the fasta mito data was aligned/called on
    fasta = batch_instance.read_input(config.config_retrieve(['references', 'mito_fasta']))

    fasta_edit_job = batch_instance.new_bash_job('Edit mito reference contig name')
    fasta_edit_job.image(config.config_retrieve(['images', 'samtools']))

    fasta_edit_job.declare_resource_group(fa={'fasta': '{root}.fasta', 'fasta.fai': '{root}.fasta.fai'})

    # edit the reference and re-index - this should be done once and stored
    fasta_edit_job.command(f'sed "s/chrM/chrMT/" {fasta} > {fasta_edit_job.fa["fasta"]}')
    fasta_edit_job.command(f'samtools faidx {fasta_edit_job.fa["fasta"]}')

    csq_job = batch_instance.new_bash_job(
        name=f'AnnotateMitoConsequenceWithBcftools: {cohort_id}',
        attributes=job_attrs | {'tool': 'bcftools'},
    )
    csq_job.image(config.config_retrieve(['images', 'bcftools']))
    csq_job.cpu(4).storage('20G')

    # to shift contig names for BCFtools
    csq_job.command('echo "chrM chrMT" > chrM_chrMT.txt')
    # to reverse the change for Hail
    csq_job.command('echo "chrMT chrM" > chMT_chrM.txt')

    csq_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    csq_job.command(f'bcftools annotate --rename-chrs chrM_chrMT.txt -Oz -o renamed.vcf.bgz -W=tbi {local_mito_vcf}')

    # the echtvar image doesn't have a tool to index, so first add that
    # then run the csq command:
    # -g is the GFF3 file
    # -B 10 is where to truncate long protein changes
    # --local-csq indicates we want each variant annotated independently (haplotype unaware)
    csq_job.command(
        f"""
        bcftools csq --force -f {fasta_edit_job.fa['fasta']} \\
            --local-csq \\
            --threads 4 \\
            -g {gff3_file} \\
            -B 10 \\
            -Oz -o annotated_chrMT.vcf.bgz \\
            -W=tbi \\
            {local_mito_vcf}
        """,
    )

    csq_job.command(
        f'bcftools annotate --rename-chrs chrMT_chrM.txt -Oz -o annotated_chrMT.vcf.bgz -W=tbi annotated_chrMT.vcf.bgz'
    )

    batch_instance.write_output(csq_job.output, str(output).removesuffix('.vcf.bgz'))

    return csq_job
