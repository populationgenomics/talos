"""
Hail Batch jobs to run VEP on a VCF in parallel.
"""

import logging
from enum import Enum
from typing import Literal, Optional, Union, Dict, Tuple, List

import hailtop.batch as hb
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_utils import to_path, Path
from cpg_utils.config import get_config
from cpg_utils.workflows.resources import STANDARD
from cpg_utils.workflows.utils import can_reuse
from cpg_utils.hail_batch import (
    copy_common_env,
    image_path,
    reference_path,
    authenticate_cloud_credentials_in_job,
    query_command,
    fasta_res_group,
)

from . import query

logger = logging.getLogger(__file__)


class SequencingType(Enum):
    """
    Type (scope) of a sequencing experiment.
    """

    GENOME = 'genome'
    EXOME = 'exome'


def vep_jobs(  # pylint: disable=too-many-arguments
    b: Batch,
    vcf_path: Path,
    tmp_prefix: Path,
    out_path: Optional[Path] = None,
    overwrite: bool = False,
    scatter_count: int = 50,
    sequencing_type: SequencingType = SequencingType.GENOME,
    job_attrs: Optional[dict] = None,
) -> List[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """
    if out_path and out_path.exists() and not overwrite:
        return []

    to_hail_table = out_path and out_path.suffix == '.ht'
    if not to_hail_table:
        assert str(out_path).endswith('.vcf.gz'), out_path

    jobs: List[Job] = []
    scatter_count = max(scatter_count, 2)
    parts_bucket = tmp_prefix / 'vep' / 'parts'
    part_files = []
    vcf = b.read_input_group(
        **{'vcf.gz': str(vcf_path), 'vcf.gz.tbi': str(vcf_path) + '.tbi'}
    )

    intervals_j, intervals = get_intervals(
        b=b,
        scatter_count=scatter_count,
        output_prefix=tmp_prefix / 'intervals',
        sequencing_type=sequencing_type,
    )
    if intervals_j:
        jobs.append(intervals_j)

    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        # find the eventual output path
        if to_hail_table:
            part_path = parts_bucket / f'part{idx + 1}.jsonl'
        else:
            part_path = parts_bucket / f'part{idx + 1}.vcf.gz'
        part_files.append(part_path)

        # here we assume that if the eventual path exists, the subset and
        # annotation were both done and can be re-used
        # the subset-vcf is not persisted, so we can skip either both jobs,
        # or neither
        if part_path and part_path.exists() and not overwrite:
            continue

        subset_j = subset_vcf(
            b,
            vcf=vcf,
            interval=intervals[idx],
            job_attrs=job_attrs or dict(part=f'{idx + 1}/{scatter_count}'),
        )
        jobs.append(subset_j)
        # To make mypy happy:
        assert isinstance(subset_j.output_vcf, hb.ResourceGroup)

        # noinspection PyTypeChecker
        vep_one_job = vep_one(
            b,
            vcf=subset_j.output_vcf['vcf.gz'],
            out_format='json' if to_hail_table else 'vcf',
            out_path=part_path,
            job_attrs=job_attrs or dict(part=f'{idx + 1}/{scatter_count}'),
            overwrite=overwrite,
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    if to_hail_table:
        gather_j = gather_vep_json_to_ht(
            b=b,
            vep_results_paths=part_files,
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
    if gather_j:
        gather_j.depends_on(*jobs)
        jobs.append(gather_j)
    return jobs


def get_intervals(
    b: hb.Batch,
    scatter_count: int,
    source_intervals_path: Optional[Path] = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
    job_attrs: Optional[dict] = None,
    output_prefix: Optional[Path] = None,
) -> tuple[Job | None, list[hb.ResourceFile]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param source_intervals_path: path to source intervals to split. Would check for
        config if not provided.
    @param sequencing_type: genome of exome,
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task:
    https://github.com/broadinstitute/warp/blob/
    bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredicted number of intervals. WDL can handle
    that, but Hail Batch is not dynamic and have to expect certain number of output
    files.
    """
    assert scatter_count > 0, scatter_count

    sequencing_type = sequencing_type or {
        'genome': SequencingType.GENOME,
        'exome': SequencingType.EXOME,
    }.get(get_config()['workflow']['sequencing_type'])

    source_intervals_path = source_intervals_path or reference_path(
        f'broad/{SequencingType.value}_calling_interval_lists'
    )

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(str(source_intervals_path))]

    if output_prefix and (output_prefix / '1.interval_list').exists():
        return None, [
            b.read_input(str(output_prefix / f'{idx + 1}.interval_list'))
            for idx in range(scatter_count)
        ]

    if (
        not source_intervals_path
        and (
            (
                existing_split_intervals_prefix := (
                    reference_path('intervals_prefix')
                    / sequencing_type
                    / f'{scatter_count}intervals'
                )
            )
            / '1.interval_list'
        ).exists()
    ):
        # We already have split intervals for this sequencing_type:
        return None, [
            b.read_input(
                str(existing_split_intervals_prefix / f'{idx + 1}.interval_list')
            )
            for idx in range(scatter_count)
        ]

    j = b.new_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | dict(tool='picard_IntervalListTools'),
    )
    j.image(image_path('picard'))
    j.memory('2G')
    j.storage('16G')
    j.cpu(2)

    break_bands_at_multiples_of = {
        SequencingType.GENOME: 100000,
        SequencingType.EXOME: 0,
    }[sequencing_type]

    cmd = f"""
    mkdir $BATCH_TMPDIR/out

    picard -Xms1000m -Xmx1500m \
    IntervalListTools \
    SCATTER_COUNT={scatter_count} \
    SUBDIVISION_MODE=INTERVAL_SUBDIVISION \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF={break_bands_at_multiples_of} \
    INPUT={b.read_input(str(source_intervals_path))} \
    OUTPUT=$BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out
    ls $BATCH_TMPDIR/out/*
    """
    for idx in range(scatter_count):
        name = f'temp_{str(idx + 1).zfill(4)}_of_{scatter_count}'
        cmd += f"""
        ln $BATCH_TMPDIR/out/{name}/scattered.interval_list {j[f'{idx + 1}.interval_list']}
        """

    j.command(cmd)
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(
                j[f'{idx + 1}.interval_list'],
                str(output_prefix / f'{idx + 1}.interval_list'),
            )

    intervals: list[hb.ResourceFile] = []
    for idx in range(scatter_count):
        interval = j[f'{idx + 1}.interval_list']
        assert isinstance(interval, hb.ResourceFile)
        intervals.append(interval)
    return j, intervals


def subset_vcf(
    b: hb.Batch,
    vcf: hb.ResourceGroup,
    interval: hb.Resource,
    job_attrs: Optional[Dict] = None,
    output_vcf_path: Optional[Path] = None,
) -> Job:
    """
    Subset VCF to provided intervals, and drop sample/genotype
    information (e.g. outputs sites-only VCF).
    """
    job_name = 'VEP: subset VCF'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    j.memory('16Gi')
    j.storage('50G')
    j.cpu(2)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
    gatk SelectVariants \\
    -R {fasta_res_group(b).base} \\
    -V {vcf['vcf.gz']} \\
    -L {interval} \\
    -O /io/batch/tmp.vcf.gz

    gatk MakeSitesOnlyVcf \\
    -I /io/batch/tmp.vcf.gz \\
    -O {j.output_vcf['vcf.gz']} \\
    --CREATE_INDEX
    """
    j.command(cmd)
    if output_vcf_path:
        b.write_output(j.output_vcf, str(output_vcf_path).replace('.vcf.gz', ''))
    return j


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceFile],
    overwrite: bool = True,
    out_vcf_path: Optional[Path] = None,
    site_only: bool = False,
    gvcf_count: Optional[int] = None,
    job_attrs: Optional[dict] = None,
) -> Tuple[Job | None, hb.ResourceGroup]:
    """
    Combines per-interval scattered VCFs into a single VCF.
    Saves the output VCF to a bucket `output_vcf_path`
    """
    if out_vcf_path and can_reuse(out_vcf_path, overwrite):
        return None, b.read_input_group(
            **{
                'vcf.gz': str(out_vcf_path),
                'vcf.gz.tbi': f'{out_vcf_path}.tbi',
            }
        )

    job_name = f'Gather {len(input_vcfs)} {"site-only " if site_only else ""}VCFs'
    j = b.new_job(job_name, job_attrs)
    j.image(image_path('gatk'))
    if gvcf_count:
        storage_gb = (1 if site_only else 2) * gvcf_count
        res = STANDARD.set_resources(j, fraction=1, storage_gb=storage_gb)
    else:
        res = STANDARD.set_resources(j, fraction=1)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    input_cmdl = ' '.join([f'--input {v}' for v in input_vcfs])
    assert isinstance(j.output_vcf, hb.ResourceGroup)
    cmd = f"""
    # --ignore-safety-checks makes a big performance difference so we include it in
    # our invocation. This argument disables expensive checks that the file headers
    # contain the same set of genotyped samples and that files are in order
    # by position of first record.
    gatk --java-options -Xms{res.get_java_mem_mb()}m \\
    GatherVcfsCloud \\
    --ignore-safety-checks \\
    --gather-type BLOCK \\
    {input_cmdl} \\
    --output {j.output_vcf['vcf.gz']}
    tabix -p vcf {j.output_vcf['vcf.gz']}
    """
    j.command(cmd)
    if out_vcf_path:
        b.write_output(j.output_vcf, str(out_vcf_path).replace('.vcf.gz', ''))
    return j, j.output_vcf


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    j = b.new_job('VEP json to Hail table', job_attrs)
    copy_common_env(j)
    j.image(get_config()['workflow']['driver_image'])
    cmd = query_command(
        query,
        query.vep_json_to_ht.__name__,
        [str(p) for p in vep_results_paths],
        str(out_path),
        setup_gcp=True,
        setup_hail=True,
    )
    j.command(cmd)
    return j


def vep_one(
    b: Batch,
    vcf: Union[Path, hb.ResourceFile],
    out_path: Optional[Path] = None,
    out_format: Literal['vcf', 'json'] = 'vcf',
    job_attrs: Optional[Dict] = None,
    overwrite: bool = False,
) -> Job | None:
    """
    Run a single VEP job.
    """
    if out_path and out_path.exists() and not overwrite:
        return None

    j = b.new_job('VEP', job_attrs)
    j.image(image_path('vep'))
    STANDARD.set_resources(j, storage_gb=50, mem_gb=50, ncpu=16)

    if not isinstance(vcf, hb.ResourceFile):
        vcf = b.read_input(str(vcf))

    if out_format == 'vcf':
        j.declare_resource_group(
            output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
        )
        assert isinstance(j.output, hb.ResourceGroup)
        output = j.output['vcf.gz']
    else:
        assert isinstance(j.output, hb.ResourceFile)
        output = j.output

    # gcsfuse works only with the root bucket, without prefix:
    vep_mount_path = reference_path('vep_mount')
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])
    loftee_conf = {
        'loftee_path': '$LOFTEE_PLUGIN_PATH',
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
    }

    authenticate_cloud_credentials_in_job(j)
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
        cmd += f'tabix -p vcf {output} '

    j.command(cmd)
    if out_path:
        b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))
    return j
