#!/usr/bin/env python3

"""
Creates a Hail Batch job to run the command line VEP tool.
"""

import hailtop.batch as hb
from hailtop.batch import Batch
from hailtop.batch.job import Job

from cpg_utils import Path, to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import (
    command,
    fasta_res_group,
    image_path,
    query_command,
    reference_path,
)


def subset_vcf(b: hb.Batch, vcf: hb.ResourceGroup, interval: hb.ResourceFile, job_attrs: dict | None = None) -> Job:
    """
    Subset VCF to provided intervals.
    """

    job_attrs = (job_attrs or {}) | {'tool': 'gatk SelectVariants'}
    j = b.new_job('Subset VCF', job_attrs)
    j.image(image_path('gatk'))
    j.storage('16Gi')
    j.memory('standard')
    j.cpu(2)

    j.declare_resource_group(output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    reference = fasta_res_group(b)
    assert isinstance(j.output_vcf, hb.ResourceGroup)

    cmd = f"""
    gatk SelectVariants \\
    -R {reference.base} \\
    -V {vcf['vcf.gz']} \\
    -L {interval} \\
    -O {j.output_vcf['vcf.gz']}
    """
    j.command(
        command(
            cmd,
            monitor_space=True,
        ),
    )
    return j


def scatter_intervals(
    b: hb.Batch,
    scatter_count: int,
    source_intervals_path: Path | None = None,
    job_attrs: dict[str, str] | None = None,
    output_prefix: Path | None = None,
) -> tuple[Job | None, list[hb.ResourceFile]]:
    """
    Add a job that splits genome/exome intervals into sub-intervals to be used to
    parallelize variant calling.

    @param b: Hail Batch object,
    @param scatter_count: number of target sub-intervals,
    @param source_intervals_path: path to source intervals to split. Would check for
        config if not provided.
    @param job_attrs: attributes for Hail Batch job,
    @param output_prefix: path optionally to save split subintervals.

    The job calls picard IntervalListTools to scatter the input interval list
    into scatter_count sub-interval lists, inspired by this WARP task :
    https://github.com/broadinstitute/warp/blob/bc90b0db0138747685b459c83ce52c8576ce03cd/tasks/broad/Utilities.wdl

    Note that we use the mode INTERVAL_SUBDIVISION instead of
    BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW. Modes other than
    INTERVAL_SUBDIVISION produce an unpredictable number of intervals. WDL can
    handle that, but Hail Batch is not dynamic and expects a certain number
    of output files.
    """
    assert scatter_count > 0, scatter_count
    sequencing_type = get_config()['workflow']['sequencing_type']
    source_intervals_path = source_intervals_path or reference_path(f'broad/{sequencing_type}_calling_interval_lists')

    if scatter_count == 1:
        # Special case when we don't need to split
        return None, [b.read_input(str(source_intervals_path))]

    if output_prefix and to_path(output_prefix / '1.interval_list').exists():
        return None, [b.read_input(str(output_prefix / f'{idx + 1}.interval_list')) for idx in range(scatter_count)]

    j = b.new_job(
        f'Make {scatter_count} intervals for {sequencing_type}',
        attributes=(job_attrs or {}) | {'tool': 'picard IntervalListTools'},
    )
    j.image(image_path('picard'))
    j.storage('16Gi')
    j.memory('2Gi')
    j.cpu(1)

    break_bands_at_multiples_of = {'genome': 100000, 'exome': 0}.get(sequencing_type, 0)

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

    j.command(command(cmd))
    if output_prefix:
        for idx in range(scatter_count):
            b.write_output(j[f'{idx + 1}.interval_list'], str(output_prefix / f'{idx + 1}.interval_list'))

    intervals: list[hb.ResourceFile] = []
    for idx in range(scatter_count):
        interval = j[f'{idx + 1}.interval_list']
        assert isinstance(interval, hb.ResourceFile)
        intervals.append(interval)
    return j, intervals


def add_vep_jobs(
    b: Batch,
    input_siteonly_vcf_path: Path,
    tmp_prefix: Path,
    scatter_count: int,
    out_path: Path | None = None,
    job_attrs: dict | None = None,
) -> list[Job]:
    """
    Runs VEP on provided VCF. Writes a VCF into `out_path` by default,
    unless `out_path` ends with ".ht", in which case writes a Hail table.
    """

    if out_path and to_path(out_path).exists():
        return []

    jobs: list[Job] = []
    siteonly_vcf = b.read_input_group(
        **{'vcf.gz': str(input_siteonly_vcf_path), 'vcf.gz.tbi': str(input_siteonly_vcf_path) + '.tbi'},
    )

    input_vcf_parts: list[hb.ResourceGroup] = []
    intervals_j, intervals = scatter_intervals(
        b=b,
        scatter_count=scatter_count,
        job_attrs=job_attrs,
        output_prefix=tmp_prefix / f'intervals_{scatter_count}',
    )
    if intervals_j:
        jobs.append(intervals_j)

    result_parts_bucket = tmp_prefix / 'vep' / 'parts'
    result_part_paths = []

    # Splitting variant calling by intervals
    for idx in range(scatter_count):
        result_part_path = result_parts_bucket / f'part{idx + 1}.jsonl'
        result_part_paths.append(result_part_path)
        if to_path(result_part_path).exists():
            continue

        subset_j = subset_vcf(
            b,
            vcf=siteonly_vcf,
            interval=intervals[idx],
            job_attrs=(job_attrs or {}) | {'part': f'{idx + 1}/{scatter_count}'},
        )
        jobs.append(subset_j)
        input_vcf_parts.append(subset_j.output_vcf)

        # noinspection PyTypeChecker
        vep_one_job = vep_one(
            b,
            vcf=input_vcf_parts[idx]['vcf.gz'],
            out_path=result_part_paths[idx],
            job_attrs=(job_attrs or {}) | {'part': f'{idx + 1}/{scatter_count}'},
        )
        if vep_one_job:
            jobs.append(vep_one_job)

    j = gather_vep_json_to_ht(
        b=b,
        vep_results_paths=result_part_paths,
        out_path=out_path,
        job_attrs=job_attrs,
        depends_on=jobs,
    )
    j.depends_on(*jobs)
    jobs.append(j)
    return jobs


def gather_vep_json_to_ht(
    b: Batch,
    vep_results_paths: list[Path],
    out_path: Path,
    job_attrs: dict | None = None,
    depends_on: list[hb.job.Job] | None = None,
) -> Job:
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table using a Batch job.
    """
    from reanalysis import vep

    j = b.new_job('VEP', job_attrs)
    j.image(get_config()['workflow']['driver_image'])
    j.command(
        query_command(
            vep,
            vep.vep_json_to_ht.__name__,
            [str(p) for p in vep_results_paths],
            str(out_path),
            setup_gcp=True,
        ),
    )
    if depends_on:
        j.depends_on(*depends_on)
    return j


def vep_one(b: Batch, vcf: Path | hb.ResourceFile, out_path: Path, job_attrs: dict | None = None) -> Job | None:
    """
    Run a single VEP job.
    """
    if out_path and to_path(out_path).exists():
        return None

    j = b.new_job('VEP', (job_attrs or {}) | {'tool': 'vep'})
    j.image(image_path('vep_110'))

    # vep is single threaded, with a middling memory requirement
    # tests have exceeded 8GB, so bump to ~13 (2 * highmem)
    j.memory('highmem')
    j.storage('10Gi')
    j.cpu(2)

    if not isinstance(vcf, hb.ResourceFile):
        vcf = b.read_input(str(vcf))

    # gcsfuse works only with the root bucket, without prefix:
    vep_mount_path = reference_path('vep_110_mount')
    data_mount = to_path(f'/{vep_mount_path.drive}')
    j.cloudfuse(vep_mount_path.drive, str(data_mount), read_only=True)
    vep_dir = data_mount / '/'.join(vep_mount_path.parts[2:])

    # assume VEP 110 has a standard install location
    loftee_conf = {
        'gerp_bigwig': f'{vep_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        'human_ancestor_fa': f'{vep_dir}/human_ancestor.fa.gz',
        'conservation_file': f'{vep_dir}/loftee.sql',
        'loftee_path': '$VEP_DIR_PLUGINS',
    }

    # sexy new plugin - only present in 110 build
    alpha_missense_plugin = f'--plugin AlphaMissense,file={vep_dir}/AlphaMissense_hg38.tsv.gz '
    j.command(
        f"""\
    FASTA={vep_dir}/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz
    vep \\
    --format vcf \\
    --json \\
    -o {j.output} \\
    -i {vcf} \\
    --everything \\
    --mane_select \\
    --allele_number \\
    --minimal \\
    --species homo_sapiens \\
    --cache --offline --assembly GRCh38 \\
    --dir_cache {vep_dir}/vep/ \\
    --fasta $FASTA \\
    {alpha_missense_plugin} \
    --plugin LoF,{','.join(f'{k}:{v}' for k, v in loftee_conf.items())}
    """,
    )

    b.write_output(j.output, str(out_path).replace('.vcf.gz', ''))

    return j
