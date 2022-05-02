"""
Creates a Hail Batch job to run the command line VEP tool.
"""


import logging
from typing import Literal, Optional

import hail as hl

# grouped imports from elastic-search-pipelines
from hail_scripts.computed_fields.variant_id import (
    get_expr_for_contig,
    get_expr_for_variant_id,
    get_expr_for_variant_ids,
    get_expr_for_xpos,
)

# custom imports to reduce line lengths
from hail_scripts.computed_fields.vep import (
    get_expr_for_vep_sorted_transcript_consequences_array as vep_array,
    get_expr_for_vep_transcript_ids_set as vep_set,
    get_expr_for_worst_transcript_consequence_annotations_struct as worst_csq,
    get_expr_for_vep_consequence_terms_set as vep_csq_set,
    get_expr_for_vep_gene_ids_set as vep_gene_id_set,
    get_expr_for_vep_protein_domains_set_from_sorted as vep_prot,
)
import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import images, utils as cpg_pipe_utils, Path, to_path
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.hb.command import wrap_command, python_command
from cpg_pipes.jobs import split_intervals
from cpg_pipes.jobs.vcf import gather_vcfs, subset_vcf
from cpg_pipes.query import vep as vep_module
from cpg_pipes.refdata import RefData
from cpg_pipes.types import SequencingType


logger = logging.getLogger(__file__)


def vep_jobs(  # pylint: disable=too-many-arguments
    b: Batch,
    vcf_path: Path,
    refs: RefData,
    hail_billing_project: str,
    hail_bucket: Path,
    tmp_bucket: Path,
    out_path: Path | None = None,
    overwrite: bool = False,
    scatter_count: int | None = RefData.number_of_vep_intervals,
    sequencing_type: SequencingType = SequencingType.GENOME,
    intervals_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    if out_path and cpg_pipe_utils.can_reuse(out_path, overwrite):
        return [b.new_job('VEP [reuse]', job_attrs)]

    scatter_count = scatter_count or RefData.number_of_vep_intervals
    jobs: list[Job] = []
    intervals_j, intervals = split_intervals.get_intervals(
        b=b,
        refs=refs,
        sequencing_type=sequencing_type,
        intervals_path=intervals_path,
        scatter_count=scatter_count,
    )
    jobs.append(intervals_j)

    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    parts_bucket = tmp_bucket / 'vep' / 'parts'
    part_files = []

    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        subset_j = subset_vcf(
            b,
            vcf=vcf,
            interval=intervals[idx],
            refs=refs,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        if to_hail_table:
            part_path = parts_bucket / f'part{idx + 1}.json_list'
        else:
            part_path = None
        # noinspection PyTypeChecker
        j = vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=part_path,
            refs=refs,
            job_attrs=(job_attrs or {}) | dict(part=f'{idx + 1}/{scatter_count}'),
            overwrite=overwrite,
        )
        jobs.append(j)
        if to_hail_table:
            part_files.append(part_path)
        else:
            part_files.append(j.output['vcf.gz'])

    if to_hail_table:
        gather_j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=part_files,
            hail_billing_project=hail_billing_project,
            hail_bucket=hail_bucket,
            out_path=out_path,
            job_attrs=job_attrs,
        )
    else:
        assert len(part_files) == scatter_count
        gather_j, _gather_vcf = gather_vcfs(
            b=b,
            input_vcfs=part_files,
            out_vcf_path=out_path,
        )
    gather_j.depends_on(*jobs)
    jobs.append(gather_j)
    return jobs


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: list[Path],
    hail_billing_project: str,
    hail_bucket: Path,
    out_path: Path,
    job_attrs: dict | None = None,
):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    j = b.new_job('VEP json to Hail table', job_attrs)
    j.image(images.DRIVER_IMAGE)
    cmd = python_command(
        vep_module,
        vep_module.vep_json_to_ht.__name__,
        [str(p) for p in vep_results_paths],
        str(out_path),
        setup_gcp=True,
        hail_billing_project=hail_billing_project,
        hail_bucket=str(hail_bucket),
        default_reference=RefData.genome_build,
    )
    j.command(cmd)
    return j


def vep_one(
    b: Batch,
    vcf: Path | hb.Resource,
    refs: RefData,
    out_path: Path | None = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Run a single VEP job.
    """
    j = b.new_job('VEP', job_attrs)
    if out_path and cpg_pipe_utils.can_reuse(out_path, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.VEP_IMAGE)
    STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    if not isinstance(vcf, hb.Resource):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(
            output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )
        output = j.output['vcf.gz']
    else:
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    base_bucket_name = refs.vep_bucket.drive
    data_mount = to_path(f'/{base_bucket_name}')
    j.cloudfuse(base_bucket_name, str(data_mount), read_only=True)
    vep_dir = data_mount / 'vep' / 'GRCh38'
    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
    }

    cmd = f"""\
    ls {vep_dir}
    ls {vep_dir}/vep

    LOFTEE_PLUGIN_PATH=/root/micromamba/share/ensembl-vep-105.0-1
    FASTA={vep_dir}/vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

    vep \\
    --format vcf \\
    --{out_format} {'--compress_output bgzip' if out_format == 'vcf' else ''} \\
    -o {output} \\
    -i {vcf} \\
    --everything \\
    --allele_number \\
    --minimal \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --dir_plugins $LOFTEE_PLUGIN_PATH \\
    --fasta $FASTA \\
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())}
    """
    if out_format == 'vcf':
        cmd += f'tabix -p vcf {output}'

    j.command(
        wrap_command(
            cmd,
            setup_gcp=True,
            monitor_space=True,
        )
    )
    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))
    return j


def apply_annotations(
    vcf_path: str,
    vep_ht_path: str,
    out_mt_path: str,
    overwrite: bool = False,
    genome_build: str = 'GRCh38',
    sequencing_type: str = 'WGS',
    checkpoints_bucket: Optional[str] = None,
):
    """
    Convert VCF to mt, add VEP annotations.
    """

    def _checkpoint(t: hl.Table, filename: str):
        """
        if a checkpoint location was provided, periodically dump the table to
        that location

        :param t:
        :param filename:
        :return:
        """
        if checkpoints_bucket:
            t = t.checkpoint(
                f'{checkpoints_bucket}/{filename}', _read_if_exists=not overwrite
            )
            logger.info(f'Wrote {checkpoints_bucket}/{filename}')
        return t

    mt = hl.import_vcf(
        str(vcf_path),
        reference_genome=genome_build,
        skip_invalid_loci=True,
        force_bgz=True,
    )
    logger.info(
        f'Importing VCF {vcf_path}, ' f'adding VEP annotations from {vep_ht_path}'
    )

    logger.info(f'Loading VEP Table from {vep_ht_path}')
    # Annotate VEP. Do before splitting multi, because we run VEP on unsplit VCF,
    # and hl.split_multi_hts can handle multiallelic VEP field.
    vep_ht = hl.read_table(str(vep_ht_path))
    logger.info(f'Adding VEP annotations into the Matrix Table from {vep_ht_path}')
    mt = mt.annotate_rows(vep=vep_ht[mt.locus].vep)

    # Splitting multi-allelics. We do not handle AS info fields here
    mt = hl.split_multi_hts(
        mt.annotate_rows(locus_old=mt.locus, alleles_old=mt.alleles)
    )
    mt = _checkpoint(mt, 'mt-vep-split.mt')

    # Add dummy GRCh37 loci
    logger.info('Adding dummy GRCh37 loci')
    mt = mt.annotate_rows(rg37_locus=mt.locus)

    seqr_ref_bucket = 'gs://cpg-seqr-reference-data'
    ref_ht_path = 'GRCh38/all_reference_data/v2/combined_reference_data_grch38-2.0.3.ht'
    ref_ht = hl.read_table(f'{seqr_ref_bucket}/{ref_ht_path}')
    clinvar_ht = hl.read_table(
        f'{seqr_ref_bucket}/GRCh38/clinvar/clinvar.GRCh38.2020-06-15.ht'
    )

    logger.info('Annotating with seqr-loader fields: round 1')
    mt = mt.annotate_rows(
        AC=mt.info.AC,
        AF=mt.info.AF[mt.a_index - 1],
        AN=mt.info.AN,
        aIndex=mt.a_index,
        wasSplit=mt.was_split,
        originalAltAlleles=get_expr_for_variant_ids(mt.locus_old, mt.alleles_old),
        sortedTranscriptConsequences=vep_array(mt.vep),
        variantId=get_expr_for_variant_id(mt),
        contig=get_expr_for_contig(mt.locus),
        pos=mt.locus.position,
        start=mt.locus.position,
        end=mt.locus.position + hl.len(mt.alleles[0]) - 1,
        ref=mt.alleles[0],
        alt=mt.alleles[1],
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(mt.locus) + hl.len(mt.alleles[0]) - 1,
        rg37_locus=mt.rg37_locus,
        clinvar_data=clinvar_ht[mt.row_key],
        ref_data=ref_ht[mt.row_key],
    )
    mt = _checkpoint(mt, 'mt-vep-split-round1.mt')

    logger.info(
        'Annotating with seqr-loader fields: round 2 '
        '(expanding sortedTranscriptConsequences, ref_data, clinvar_data)'
    )
    mt = mt.annotate_rows(
        domains=vep_prot(mt.sortedTranscriptConsequences),
        transcriptConsequenceTerms=vep_csq_set(mt.sortedTranscriptConsequences),
        transcriptIds=vep_set(mt.sortedTranscriptConsequences),
        mainTranscript=worst_csq(mt.sortedTranscriptConsequences),
        geneIds=vep_gene_id_set(mt.sortedTranscriptConsequences),
        codingGeneIds=vep_gene_id_set(
            mt.sortedTranscriptConsequences, only_coding_genes=True
        ),
        cadd=mt.ref_data.cadd,
        dbnsfp=mt.ref_data.dbnsfp,
        geno2mp=mt.ref_data.geno2mp,
        gnomad_exomes=mt.ref_data.gnomad_exomes,
        gnomad_exome_coverage=mt.ref_data.gnomad_exome_coverage,
        gnomad_genomes=mt.ref_data.gnomad_genomes,
        gnomad_genome_coverage=mt.ref_data.gnomad_genome_coverage,
        eigen=mt.ref_data.eigen,
        exac=mt.ref_data.exac,
        g1k=mt.ref_data.g1k,
        mpc=mt.ref_data.mpc,
        primate_ai=mt.ref_data.primate_ai,
        splice_ai=mt.ref_data.splice_ai,
        topmed=mt.ref_data.topmed,
        clinvar=hl.struct(
            **{
                'allele_id': mt.clinvar_data.info.ALLELEID,
                'clinical_significance': hl.delimit(mt.clinvar_data.info.CLNSIG),
                'gold_stars': mt.clinvar_data.gold_stars,
            }
        ),
    )
    mt = mt.annotate_globals(
        sourceFilePath=vcf_path,
        genomeVersion=genome_build.replace('GRCh', ''),
        sampleType=sequencing_type,
        hail_version=hl.version(),
    )
    logger.info('Done:')
    mt.describe()
    mt.write(str(out_mt_path), overwrite=overwrite)
    logger.info(f'Written final matrix table into {out_mt_path}')
