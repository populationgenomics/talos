#!/usr/bin/env python3


"""
Entrypoint for clinvar summary generation
"""

from datetime import datetime
from os.path import join

from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path, reference_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch, query_command

from reanalysis import clinvar_by_codon, seqr_loader, summarise_clinvar_entries
from reanalysis.static_values import get_logger
from reanalysis.vep_jobs import add_vep_jobs


def generate_clinvar_table(cloud_folder: Path, clinvar_outputs: str):
    """
    set up the job that does de novo clinvar summary

    Args:
        cloud_folder (Path): folder for this analysis
        clinvar_outputs (str): prefix for writing new files/dirs
    """

    bash_job = get_batch().new_bash_job(name='copy clinvar files to local')
    bash_job.image(config_retrieve(['workflow', 'driver_image']))

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(
        f'wget -q {directory}{sub_file} -O {bash_job.subs} && wget -q {directory}{var_file} -O {bash_job.vars}',
    )

    # write output files date-specific
    get_batch().write_output(bash_job.subs, str(cloud_folder / sub_file))
    get_batch().write_output(bash_job.vars, str(cloud_folder / var_file))

    # region: run the summarise_clinvar_entries script
    summarise = get_batch().new_job(name='summarise clinvar')
    summarise.depends_on(bash_job)

    summarise.cpu(2).image(config_retrieve(['workflow', 'driver_image'])).storage('20G')
    authenticate_cloud_credentials_in_job(summarise)
    summarise.command(
        f'python3 {summarise_clinvar_entries.__file__} '
        f'-s {bash_job.subs} '
        f'-v {bash_job.vars} '
        f'-o {clinvar_outputs}',
    )

    return summarise


def vep_json_to_ht(json_paths: list, output_ht: str):
    """
    call as a python Job. Take all the JSON files and make a Hail Table?

    Args:
        json_paths ():
        output_ht ():

    Returns:

    """
    import hail as hl

    hl.init()
    hl.default_reference('GRCh38')

    # this will need some tweaking
    json_schema = hl.dtype(
        """struct{
        minimised:int32,
        assembly_name:str,
        allele_string:str,
        ancestral:str,
        context:str,
        end:int32,
        id:str,
        input:str,
        seq_region_name:str,
        start:int32,
        strand:int32,
        transcript_consequences:array<struct{
            allele_num:int32,
            amino_acids:str,
            appris:str,
            biotype:str,
            canonical:int32,
            cdna_start:int32,
            cdna_end:int32,
            cds_end:int32,
            cds_start:int32,
            codons:str,
            consequence_terms:array<str>,
            distance:int32,
            exon:str,
            gene_id:str,
            gene_pheno:int32,
            gene_symbol:str,
            gene_symbol_source:str,
            hgnc_id:str,
            hgvsc:str,
            hgvsp:str,
            hgvs_offset:int32,
            impact:str,
            intron:str,
            minimised:int32,
            mirna:array<str>,
            protein_end:int32,
            strand:int32,
            swissprot:array<str>,
            transcript_id:str,
            trembl:array<str>,
            tsl:int32,
            uniparc:array<str>,
            uniprot_isoform:array<str>,
            variant_allele:str,
            source:str,
            flags:array<str>
        }>,
        variant_class:str
    }""",
    )
    ht = hl.import_table(paths=json_paths, no_header=True, types={'f0': json_schema})
    ht = ht.transmute(vep=ht.f0)

    # Can't use ht.vep.start for start because it can be modified by VEP (e.g. it
    # happens for indels). So instead parsing POS from the original VCF line stored
    # as ht.vep.input field.
    input_split = ht.vep.input.split('\t')
    start = hl.parse_int(input_split[1])
    chrom = ht.vep.seq_region_name
    ht = ht.annotate(locus=hl.locus(chrom, start), alleles=[input_split[3], input_split[4]])
    ht = ht.key_by(ht.locus, ht.alleles)
    ht.write(output_ht)


def generate_annotated_data(annotation_out: Path, snv_vcf: str, tmp_path: Path, dependency: Job | None = None) -> Job:
    """
    if the annotated data Table doesn't exist, generate it

    Args:
        annotation_out (str): MT path to create
        snv_vcf (str): path to a VCF file
        tmp_path (str): path to temp directory
        dependency (Job | None): optional job dependency

    Returns:
        The Job for future dependency setting
    """

    snv_vcf_in_batch = get_batch().read_input_group(**{'vcf.gz': snv_vcf, 'vcf.gz.tbi': f'{snv_vcf}.tbi'})['vcf.gz']

    # split the whole vcf into chromosomes
    output_json_files: list = []
    for chromosome in [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY']:
        # subset the whole VCF to this chromosome
        bcftools_job = get_batch().new_job(f'subset_{chromosome} with bcftools')
        if dependency:
            bcftools_job.depends_on(dependency)

        # set some resources
        bcftools_job.image(image_path('bcftools')).cpu(1).memory('8G')
        # create a vcf fragment
        bcftools_job.command(f'bcftools view -Oz -o {bcftools_job.fragment} -r {chromosome} {snv_vcf_in_batch}')
        # annotate that fragment, making a JSON output
        vep_job = get_batch().new_job(f'annotate_{chromosome} with VEP')
        # configure the required resources
        vep_job.image(image_path('vep_110')).cpu(1).memory('highmem')

        # gcsfuse works only with the root bucket, without prefix:
        vep_mount_path = to_path(reference_path('vep_110_mount'))
        data_mount = to_path(f'/{vep_mount_path.drive}')
        vep_job.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
        vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

        vep_job.command(
            f"""\
            FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz
            vep --format vcf --json -o {vep_job.output} -i {bcftools_job.fragment} \\
            --protein --species homo_sapiens --cache --offline --assembly GRCh38 \\
            --dir_cache {vep_dir}/vep/ --fasta $FASTA
            """,
        )

        get_batch().write_output(vep_job.output, str(tmp_path / f'{chromosome}.json'))
        output_json_files.append(vep_job.output)

    # call a python job to stick all those together?!
    json_to_mt_job = get_batch().new_python_job('aggregate JSON into MT')
    json_to_mt_job.call(vep_json_to_ht, output_json_files, str(annotation_out))
    return json_to_mt_job


def main():
    """
    run the clinvar summary, output to common path
    """

    cloud_folder = to_path(
        join(config_retrieve(['storage', 'common', 'analysis']), 'aip_clinvar_new', datetime.now().strftime('%y-%m')),
    )

    # clinvar VCF, decisions, annotated VCF, and PM5
    clinvar_output_path = join(str(cloud_folder), 'clinvar_decisions')
    clinvar_ht = f'{clinvar_output_path}.ht'
    snv_vcf = f'{clinvar_output_path}.vcf.bgz'
    clinvar_pm5_path = cloud_folder / 'clinvar_pm5.ht'
    annotated_clinvar = cloud_folder / 'annotated_clinvar.mt'

    # check if we can just quit already
    if all(this_path.exists() for this_path in [annotated_clinvar, clinvar_ht, clinvar_pm5_path]):
        get_logger().info('Clinvar data already exists, exiting')
        return

    temp_path = to_path(
        join(config_retrieve(['storage', 'common', 'tmp']), 'aip_clinvar_new', datetime.now().strftime('%y-%m')),
    )

    dependency = None

    # generate a new round of clinvar decisions
    if not all(to_path(output).exists() for output in [clinvar_ht, snv_vcf]):
        dependency = generate_clinvar_table(cloud_folder, clinvar_output_path)

    # create the annotation job(s)
    if not annotated_clinvar.exists():
        dependency = generate_annotated_data(annotated_clinvar, snv_vcf, temp_path, dependency=dependency)

    # region: run the clinvar_by_codon script
    if not clinvar_pm5_path.exists():
        clinvar_by_codon_job = get_batch().new_job(name='clinvar_by_codon')
        clinvar_by_codon_job.image(config_retrieve(['workflow', 'driver_image'])).cpu(2).storage('20G')
        authenticate_cloud_credentials_in_job(clinvar_by_codon_job)
        clinvar_by_codon_job.command(
            f'python3 {clinvar_by_codon.__file__} --mt_path {annotated_clinvar} --write_path {clinvar_pm5_path}',
        )
        if dependency:
            clinvar_by_codon_job.depends_on(dependency)
    # endregion

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
