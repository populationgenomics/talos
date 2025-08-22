"""
This is a central script for the Talos process, implemented at the CPG, using the cpg-flow workflow framework
"""

import datetime
import functools
from os.path import join
from random import randint
from typing import TYPE_CHECKING

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config, hail_batch, to_path
from loguru import logger

from talos.cpg_internal_scripts.annotation_stages import TransferAnnotationsToMt
from talos.cpg_internal_scripts.cpg_flow_utils import generate_dataset_prefix, query_for_latest_analysis
from talos.cpg_internal_scripts.cpgflow_jobs import MakeConfig
from talos.utils import get_granular_date

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


SV_ANALYSIS_TYPES = {
    'exome': 'single_dataset_cnv_annotated',
    'genome': 'single_dataset_sv_annotated',
}


def set_up_job_with_resources(
    name: str,
    memory: str | None = None,
    cpu: float | None = None,
    storage: str = '10Gi',
    image: str | None = None,
    attrs: dict | None = None,
) -> 'BashJob':
    """
    Wrapper to create a job with all elements set up
    Name is mandatory, the rest is optional

    Args:
        name (str): name of the Job
        memory (str): optional, defaults to 'standard'. Can be exact ("12G") or a category ("highmem")
        cpu (float): optional, defaults to 2.0
        storage (str): optional, defaults to 10Gi
        image (str): optional, full path to Docker image to use
        attrs (dict): optional, attributes to add to the job

    Returns:
        A job in the current Batch with all resources allocated
    """

    job = hail_batch.get_batch().new_job(name=name, attributes=attrs)
    if image:
        job.image(image)
    else:
        job.image(config.config_retrieve(['workflow', 'driver_image']))
    if cpu:
        job.cpu(cpu)
    if memory:
        job.memory(memory)
    if storage:
        job.storage(storage)

    return job


def tshirt_mt_sizing(sequencing_type: str, cohort_size: int) -> str:
    """
    Some way of taking the details we have (#SGs, sequencing type)
    and producing an estimate (with padding) of the MT size on disk
    used to determine VM provision

    This method is being copied here and modified - Talos is using minimised MTs, so the size is much smaller
    The current calculation method in cpg-workflows/flow is still useful in the Seqr pipeline

    Returns a String, the value for job.storage(X)
    """

    if (sequencing_type == 'genome' and cohort_size < 400) or (sequencing_type == 'exome' and cohort_size < 3000):
        return '50Gi'
    return '250Gi'


@functools.cache
def get_date_string() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses

    Returns:
        either an override in config, or the default (today, YYYY-MM-DD)
    """
    return config.config_retrieve(['workflow', 'date_folder_override'], get_granular_date())


@stage.stage(analysis_type='panelapp')
class DownloadPanelAppData(stage.MultiCohortStage):
    """
    runs a single instance of the stage which downloads the whole of PanelApp into a cached file
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return to_path(
            join(
                config.config_retrieve(['storage', 'common', 'analysis']),
                'panelapp_monthly',
                f'panelapp_data_{datetime.datetime.now().strftime("%y-%m")}.json',  # noqa: DTZ005
            ),
        )

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        inputs: stage.StageInput,
    ) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)

        # get the MANE json file - used to map gene Symbols <-> IDs
        mane_json = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'mane_1.4', 'json']))
        job = set_up_job_with_resources(name='DownloadPanelAppData', cpu=1)
        job.command(f'DownloadPanelApp --output {job.output} --mane {mane_json}')

        hail_batch.get_batch().write_output(job.output, output)

        return self.make_outputs(target=multicohort, data=output, jobs=job)


@stage.stage
class MakeRuntimeConfig(stage.CohortStage):
    """
    create a config file for this run,
    this new config should include all elements specific to this Project and sequencing_type
    this new unambiguous config file should be used in all jobs
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        prefix = generate_dataset_prefix(
            dataset=cohort.dataset.name,
            stage_name=self.name,
            hash_value=get_date_string(),
        )
        return {
            'config': prefix / 'config.toml',
            'seqr_lookup': prefix / 'seqr_lookup.json',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        expected_outputs = self.expected_outputs(cohort)

        MakeConfig.create_config(
            cohort=cohort,
            seqr_out=expected_outputs['seqr_lookup'],
            config_out=expected_outputs['config'],
        )
        return self.make_outputs(target=cohort, data=expected_outputs)


@stage.stage
class MakeHpoPedigree(stage.CohortStage):
    """Generate a pedigree from metamist - additional column for HPO terms"""

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'pedigree.ped'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        job = set_up_job_with_resources(name=f'MakeHpoPedigree: {cohort.id} ({cohort.dataset.name})', cpu=1)

        expected_out = self.expected_outputs(cohort)
        query_dataset = cohort.dataset.name
        if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
            query_dataset += '-test'

        # mandatory argument
        seq_type = config.config_retrieve(['workflow', 'sequencing_type'])

        # insert a little stagger
        job.command(f'sleep {randint(0, 180)}')

        tech_string = '--tech long-read' if config.config_retrieve(['workflow', 'long_read'], False) else ''
        job.command(
            f"""
            python -m talos.cpg_internal_scripts.MakeHpoPedigree \\
                --dataset {query_dataset} \\
                --output {job.output} \\
                --type {seq_type} {tech_string}
            """,
        )
        hail_batch.get_batch().write_output(job.output, expected_out)
        logger.info(f'Pedigree file for {cohort.id} ({cohort.dataset.name}) going to {expected_out!s}')

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(required_stages=[DownloadPanelAppData, MakeRuntimeConfig, MakeHpoPedigree])
class UnifiedPanelAppParser(stage.CohortStage):
    """
    Job to parse the PanelApp data, output specific to this Dataset
    """

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'panelapp_data.json'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        # create and resource a new job
        job = set_up_job_with_resources(name=f'UnifiedPanelAppParser: {cohort.id} ({cohort.dataset.name})', cpu=1)

        # read in the config made in MakeRuntimeConfig
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))

        # read in the output from DownloadPanelAppData
        panelapp_download = hail_batch.get_batch().read_input(
            inputs.as_path(workflow.get_multicohort(), DownloadPanelAppData),
        )

        local_ped = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeHpoPedigree))

        hpo_file = hail_batch.get_batch().read_input(config.config_retrieve(['GeneratePanelData', 'obo_file']))

        expected_out = self.expected_outputs(cohort)

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.UnifiedPanelAppParser \\
                --input {panelapp_download} \\
                --output {job.output} \\
                --pedigree {local_ped} \\
                --hpo {hpo_file}
            """,
        )

        hail_batch.get_batch().write_output(job.output, expected_out)

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[
        UnifiedPanelAppParser,
        MakeHpoPedigree,
        MakeRuntimeConfig,
        TransferAnnotationsToMt,
    ],
)
class RunHailFiltering(stage.CohortStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'hail_labelled.vcf.bgz'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        # integrate this into the earlier workflow
        input_mt = inputs.as_path(cohort, TransferAnnotationsToMt)

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))

        # Talos-prep workflow does a lot of compression, which is nice
        # MTs can vary from <10GB for a small exome, to 170GB for a larger one, Genomes ~200GB
        storage_estimate = tshirt_mt_sizing(
            sequencing_type=config.config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(cohort.get_sequencing_group_ids()),
        )

        job = set_up_job_with_resources(
            name=f'RunHailFiltering: {cohort.id} ({cohort.dataset.name})',
            storage=config.config_retrieve(['RunHailFiltering', 'storage', 'small_variants'], storage_estimate),
            cpu=16,
            memory='highmem',
        )
        job.command('set -eux pipefail')

        panelapp_json = hail_batch.get_batch().read_input(inputs.as_path(cohort, UnifiedPanelAppParser))

        # peds can't read cloud paths
        pedigree = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeHpoPedigree))
        expected_out = self.expected_outputs(cohort)

        # copy vcf & index out manually
        job.declare_resource_group(
            output={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # find the clinvar table, localise, and expand
        if not (clinvar_tar := config.config_retrieve(['workflow', 'clinvar_data'], None)):
            clinvar_tar = query_for_latest_analysis(
                dataset=config.config_retrieve(['workflow', 'dataset']),
                analysis_type='clinvarbitration',
            )
            if clinvar_tar is None:
                raise ValueError('No ClinVar data found')

        job.command(f'tar -xzf {hail_batch.get_batch().read_input(clinvar_tar)} -C $BATCH_TMPDIR')

        # read in the MT using gcloud, directly into batch tmp
        job.command(f'gcloud storage cp -r {input_mt!s} $BATCH_TMPDIR')

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.RunHailFiltering \\
                --input "${{BATCH_TMPDIR}}/{cohort.id}.mt" \\
                --panelapp {panelapp_json} \\
                --pedigree {pedigree} \\
                --output {job.output['vcf.bgz']} \\
                --clinvar "${{BATCH_TMPDIR}}/clinvarbitration_data/clinvar_decisions.ht" \\
                --pm5 "${{BATCH_TMPDIR}}/clinvarbitration_data/clinvar_decisions.pm5.ht" \\
                --checkpoint "${{BATCH_TMPDIR}}/checkpoint.mt"
            """,
        )
        hail_batch.get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(required_stages=[UnifiedPanelAppParser, MakeHpoPedigree, MakeRuntimeConfig])
class RunHailFilteringSv(stage.CohortStage):
    """
    hail job to filter & label the SV MT
    """

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        if (
            query_for_latest_analysis(
                dataset=cohort.dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
                long_read=config.config_retrieve(['workflow', 'long_read'], False),
            )
            is not None
        ):
            return (
                generate_dataset_prefix(
                    dataset=cohort.dataset.name,
                    stage_name=self.name,
                    hash_value=get_date_string(),
                )
                / 'labelled_svs.vcf.bgz'
            )
        return {}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        # early skip if the stage has nothing to run on
        expected_out = self.expected_outputs(cohort)
        if (
            path_or_none := query_for_latest_analysis(
                dataset=cohort.dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
                long_read=config.config_retrieve(['workflow', 'long_read'], False),
            )
        ) is None:
            logger.info(f'No SV MT found for {cohort.id} ({cohort.dataset.name}), skipping')
            return self.make_outputs(cohort, data=expected_out)

        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))
        panelapp_json = hail_batch.get_batch().read_input(inputs.as_path(cohort, UnifiedPanelAppParser))
        pedigree = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeHpoPedigree))

        cpu: int = config.config_retrieve(['RunHailFiltering', 'cores', 'sv'], 2)
        job = set_up_job_with_resources(
            name=f'RunHailFilteringSV: {cohort.id} ({cohort.dataset.name}), {path_or_none}',
            cpu=cpu,
            memory='highmem',
        )

        # manually extract the VCF and index after the job completes
        job.declare_resource_group(
            output={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # get the MANE json file - used to map gene Symbols <-> IDs
        mane_json = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'mane_1.4', 'json']))

        # copy the VCF in
        annotated_vcf = hail_batch.get_batch().read_input_group(
            **{
                'vcf.bgz': path_or_none,
                'vcf.bgz.tbi': f'{path_or_none}.tbi',
            },
        )['vcf.bgz']
        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.RunHailFilteringSv \\
                --input {annotated_vcf} \\
                --panelapp {panelapp_json} \\
                --pedigree {pedigree} \\
                --mane_json {mane_json} \\
                --output {job.output['vcf.bgz']}
            """,
        )
        hail_batch.get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[
        MakeHpoPedigree,
        UnifiedPanelAppParser,
        RunHailFiltering,
        RunHailFilteringSv,
        MakeRuntimeConfig,
    ],
)
class ValidateVariantInheritance(stage.CohortStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'summary_output.json'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        # resource consumption here has dropped hugely
        job = set_up_job_with_resources(name=f'ValidateMOI: {cohort.id} ({cohort.dataset.name})', cpu=2.0)

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))

        panelapp_data = hail_batch.get_batch().read_input(inputs.as_path(cohort, UnifiedPanelAppParser))

        pedigree = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeHpoPedigree))
        hail_inputs = inputs.as_path(cohort, RunHailFiltering)

        # either find a SV vcf, or None
        if (
            query_for_latest_analysis(
                dataset=cohort.dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
                long_read=config.config_retrieve(['workflow', 'long_read'], False),
            )
            is not None
        ):
            hail_sv_inputs = inputs.as_path(cohort, RunHailFilteringSv)
            labelled_sv_vcf = hail_batch.get_batch().read_input_group(
                **{'vcf.bgz': hail_sv_inputs, 'vcf.bgz.tbi': f'{hail_sv_inputs}.tbi'},
            )['vcf.bgz']
            sv_vcf_arg = f'--labelled_sv {labelled_sv_vcf} '
        else:
            sv_vcf_arg = ''

        labelled_vcf = hail_batch.get_batch().read_input_group(
            **{
                'vcf.bgz': hail_inputs,
                'vcf.bgz.tbi': f'{hail_inputs}.tbi',
            },
        )['vcf.bgz']

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.ValidateMOI \\
                --labelled_vcf {labelled_vcf} \\
                --output {job.output} \\
                --panelapp {panelapp_data} \\
                --pedigree {pedigree} {sv_vcf_arg}
            """,
        )
        expected_out = self.expected_outputs(cohort)
        hail_batch.get_batch().write_output(job.output, expected_out)

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[MakeRuntimeConfig, ValidateVariantInheritance],
    analysis_type='aip-results',
    analysis_keys=['report'],
)
class HpoFlagging(stage.CohortStage):
    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'full_report.json'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        phenio_db = hail_batch.get_batch().read_input(config.config_retrieve(['HPOFlagging', 'phenio_db']))
        gene_to_phenotype = hail_batch.get_batch().read_input(
            config.config_retrieve(['HPOFlagging', 'gene_to_phenotype']),
        )

        job = set_up_job_with_resources(
            name=f'HPOFlagging: {cohort.id} ({cohort.dataset.name})',
            cpu=2,
            memory='highmem',
            storage='20Gi',
        )

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, key='config'))

        results_json = hail_batch.get_batch().read_input(
            inputs.as_path(cohort, ValidateVariantInheritance),
        )

        mane_json = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'mane_1.4', 'json']))

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.HPOFlagging \\
                --input {results_json} \\
                --mane_json {mane_json} \\
                --gen2phen {gene_to_phenotype} \\
                --phenio {phenio_db} \\
                --output {job.output!s}
            """,
        )

        hail_batch.get_batch().write_output(job.output, outputs)

        return self.make_outputs(target=cohort, jobs=job, data=outputs)


@stage.stage(
    required_stages=[HpoFlagging, UnifiedPanelAppParser, MakeRuntimeConfig],
    analysis_type='aip-report',
    analysis_keys=['dated'],
    tolerate_missing_output=True,
)
class CreateTalosHtml(stage.CohortStage):
    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        std_prefix = generate_dataset_prefix(
            dataset=cohort.dataset.name,
            stage_name=self.name,
            hash_value=get_date_string(),
        )
        web_prefix = generate_dataset_prefix(
            dataset=cohort.dataset.name,
            category='web',
            stage_name=self.name,
            hash_value=get_date_string(),
        )
        static_web_prefix = generate_dataset_prefix(
            dataset=cohort.dataset.name,
            category='web',
            stage_name=self.name,
            hash_value='talos_static',
        )
        return {
            'tar': std_prefix / 'reports.tar.gz',
            'id_map': std_prefix / 'int_ext_id_map.tsv',
            'dated': web_prefix / 'summary_output.html',
            'generic': static_web_prefix / 'summary_output.html',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        job = set_up_job_with_resources(name=f'CreateTalosHtml: {cohort.id} ({cohort.dataset.name})')

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))

        expected_out = self.expected_outputs(cohort)

        # generate an internal~external ID map for labelling the HTML file
        with expected_out['id_map'].open('w') as id_map_handle:
            for sg in cohort.get_sequencing_groups():
                id_map_handle.write(f'{sg.id}\t{sg.participant_id}\n')

        localised_ids = hail_batch.get_batch().read_input(expected_out['id_map'])
        results_json = hail_batch.get_batch().read_input(inputs.as_str(cohort, HpoFlagging))
        panelapp_data = hail_batch.get_batch().read_input(inputs.as_path(cohort, UnifiedPanelAppParser))

        # this will write output files directly to GCP
        job.command(f'export TALOS_CONFIG={runtime_config}')

        # create a new directory for the results
        job.command('mkdir html_outputs')
        job.command('cd html_outputs')
        job.command(
            f"""
            python -m talos.CreateTalosHTML \\
                --input {results_json} \\
                --panelapp {panelapp_data} \\
                --output summary_output.html \\
                --ext_ids {localised_ids}
            """,
        )

        # copy up to a date-specific run folder
        date_folder = str(expected_out['dated']).removesuffix('/summary_output.html')
        job.command(f'gcloud storage cp -r * {date_folder!s}')

        # copy to a dataset-generic folder
        generic_folder = str(expected_out['generic']).removesuffix('/summary_output.html')
        job.command(f'gcloud storage cp -r * {generic_folder!s}')

        # Create a tar'chive here, then use an image with GCloud to copy it up in a bit
        job.command(f'tar -czf {job.output} *')

        hail_batch.get_batch().write_output(job.output, expected_out['tar'])

        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[HpoFlagging, MakeRuntimeConfig],
    analysis_keys=['seqr_file', 'seqr_pheno_file'],
    analysis_type='custom',
    tolerate_missing_output=True,
)
class MinimiseOutputForSeqr(stage.CohortStage):
    """
    takes the results file from Seqr and produces a minimised form for Seqr ingestion
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        analysis_prefix = (
            generate_dataset_prefix(
                dataset=cohort.dataset.name,
                category='analysis',
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'seqr_files'
        )
        return {
            'seqr_file': analysis_prefix / f'{get_date_string()}_seqr.json',
            'seqr_pheno_file': analysis_prefix / f'{get_date_string()}_seqr_pheno.json',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        # pull out the config section relevant to this datatype & cohort
        # if it doesn't exist for this sequencing type, fail gracefully
        seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
        try:
            seqr_lookup = config.config_retrieve(['cohorts', cohort.dataset.name, seq_type, 'seqr_lookup'])
        except config.ConfigError:
            logger.warning(f'No Seqr lookup file for {cohort.id} ({cohort.dataset.name}) {seq_type}')
            return self.make_outputs(cohort, skipped=True)

        input_localised = hail_batch.get_batch().read_input(inputs.as_str(cohort, HpoFlagging))

        # create a job to run the minimisation script
        job = set_up_job_with_resources(
            name=f'MinimiseOutputForSeqr: {cohort.id} ({cohort.dataset.name})',
            memory='lowmem',
        )

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(cohort, MakeRuntimeConfig, 'config'))

        lookup_in_batch = hail_batch.get_batch().read_input(seqr_lookup)
        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            f"""
            python -m talos.cpg_internal_scripts.MinimiseOutputForSeqr \\
                --input {input_localised} \\
                --output {job.out_json} \\
                --pheno {job.pheno_json} \\
                --external_map {lookup_in_batch}
            """,
        )

        # write the results out
        expected_out = self.expected_outputs(cohort)
        hail_batch.get_batch().write_output(job.out_json, expected_out['seqr_file'])
        hail_batch.get_batch().write_output(job.pheno_json, expected_out['seqr_pheno_file'])
        return self.make_outputs(cohort, data=expected_out, jobs=job)


@stage.stage(required_stages=[CreateTalosHtml])
class UpdateIndexFile(stage.MultiCohortStage):
    """Update the index.html with all the latest reports."""

    def expected_outputs(self, cohort: targets.Cohort) -> Path:
        return (
            generate_dataset_prefix(
                dataset=config.config_retrieve(['workflow', 'dataset']),
                stage_name=self.name,
                hash_value=get_date_string(),
            )
            / 'index_created.txt'
        )

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        job = hail_batch.get_batch().new_bash_job('Create Index page')
        job.image(config.config_retrieve(['workflow', 'driver_image'])).memory('lowmem').cpu(1)

        # write out the index file
        expected_out = self.expected_outputs(cohort)

        dataset = config.config_retrieve(['workflow', 'dataset'])
        job.command(f'python -m talos.cpg_internal_scripts.BuildReportIndexPage --dataset {dataset}')

        return self.make_outputs(cohort, data=expected_out, jobs=job)
