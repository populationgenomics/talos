from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


CLINVAR_FTP = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited'


def run_clinvarbitration_in_full(
    clinvar_file_tmp: Path,
    decisions: Path,
    pm5: Path,
) -> 'BashJob':
    """
    gets the remote resources for submissions and variants

    Args:
        clinvar_file_tmp (Pathlike): temp directory for clinvar files
        decisions (Pathlike): ht for all decisions
        pm5 (Pathlike): ht for all PM5 annotations
    """

    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_bash_job('Run ClinvArbitration').storage('10Gi').memory('highmem').cpu('4')
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command('set -eo pipefail')

    submission_name = 'submission_summary.txt.gz'
    variant_name = 'variant_summary.txt.gz'

    # region: download clinvar source files
    # for each input, if the file exists in GCP temp, use it, otherwise download a fresh one, write it to temp and GCP
    for filename in [submission_name, variant_name]:
        if not (clinvar_file_tmp / filename).exists():
            # download with wget, write to stdout, use tee to write that to GCP temp, and a local file. Nice.
            job.command(
                f'wget -q {CLINVAR_FTP}/{filename} -O - | tee ${{BATCH_TMPDIR}}/{filename} | gcloud storage cp - {clinvar_file_tmp}/{filename}',  # noqa: E501
            )
        else:
            temp_file = batch_instance.read_input(clinvar_file_tmp / filename)
            job.command(f'mv {temp_file} ${{BATCH_TMPDIR}}/{filename}')
    # endregion

    # region: make new summary
    if sites_to_blacklist := config.config_retrieve(['workflow', 'site_blacklist'], []):
        blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
        blacklist_string = f' -b {blacklist_sites}'
    else:
        blacklist_string = ''

    job.command(f"""
        python3 -m clinvarbitration.scripts.resummarise_clinvar \\
            -v ${{BATCH_TMPDIR}}/{variant_name} \\
            -s ${{BATCH_TMPDIR}}/{submission_name} \\
            {blacklist_string} -o ${{BATCH_TMPDIR}}/clinvarbitration
    """)
    # copy up to GCP
    job.command(f'gcloud storage cp -r ${{BATCH_TMPDIR}}/clinvarbitration.ht {decisions}')
    # endregion

    # region: annotate SNVs
    # read in the input files required for bcftools csq
    gff3 = batch_instance.read_input(config.config_retrieve(['references', 'ensembl_gff3']))
    ref_fa = batch_instance.read_input(config.config_retrieve(['references', 'ref_fasta']))
    # -g is the GFF3 file, -f is the reference fasta
    # --local-csq is required to apply non-phase aware annotation
    # --force is required to use annotations without phase data
    job.command(f"""
        bcftools csq \\
        --force \\
        --local-csq \\
        -f {ref_fa} \\
        --unify-chr-names 'chr,-,chr' \\
        -g {gff3} \\
        "${{BATCH_TMPDIR}}/clinvarbitration.vcf.bgz" \\
        -o out.vcf
    """)

    # split the bcftools CSQ fields, filter to missense, and write out a tab-delimited file
    # -d - duplicate, writes each transcript onto a new line
    # -s :missense - only keep Consequence==missense variants
    # -f - format string - tab delimited, Transcript, Amino Acid Change, ClinVar allele ID, ClinVar gold stars
    job.command(
        """
        bcftools \\
            +split-vep out.vcf \\
            -d \\
            -s :missense \\
            -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \\
            > ${{BATCH_TMPDIR}}/clinvarbitration.annotated.tsv
        """,
    )
    # endregion

    # region: generate PM5 data
    # write both HT and TSV outputs to the same root location
    job.command("""
        python3 -m clinvarbitration.scripts.clinvar_by_codon \\
            -i ${BATCH_TMPDIR}/clinvarbitration.annotated.tsv \\
            -o ${BATCH_TMPDIR}/clinvarbitration.pm5
    """)
    # copy up
    job.command(f'gcloud storage cp -r ${{BATCH_TMPDIR}}/clinvarbitration.pm5.ht {pm5}')
    # endregion

    return job
