"""
This is a central script for the Talos process, implemented at the CPG, using the cpg-flow workflow framework
"""

import toml
import datetime
import functools

from loguru import logger
from os.path import join
from random import randint
from typing import TYPE_CHECKING

from cpg_flow import stage, targets, workflow
from cpg_utils import Path, to_path, config, hail_batch

from talos.utils import get_granular_date
from talos.cpg_internal_scripts.prep_stages import SquashMtIntoTarballStage
from talos.cpg_internal_scripts.cpg_flow_utils import query_for_latest_analysis

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

    Returns:
        A job in the current Batch with all resources allocated
    """

    job = hail_batch.get_batch().new_job(name=name)
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


def tshirt_mt_sizing(sequencing_type: str, cohort_size: int) -> int:
    """
    Some way of taking the details we have (#SGs, sequencing type)
    and producing an estimate (with padding) of the MT size on disk
    used to determine VM provision

    This method is being copied here and modified - Talos is using minimised MTs, so the size is much smaller
    The current calculation method in cpg-workflows/flow is still useful in the Seqr pipeline

    Args:
        sequencing_type (str): exome or genome
        cohort_size (int): number of samples

    Returns:
        str, the value for job.storage(X)
    """

    if (sequencing_type == 'genome' and cohort_size < 400) or (sequencing_type == 'exome' and cohort_size < 3000):  # noqa: PLR2004
        return 50
    return 250


@functools.cache
def get_date_string() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses

    Returns:
        either an override in config, or the default (today, YYYY-MM-DD)
    """
    return config.config_retrieve(['workflow', 'date_folder_override'], get_granular_date())


@functools.cache
def get_date_folder() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses
    Returns:
        either an override in config, or the default "reanalysis/(today, YYYY-MM-DD)"
    """
    return join('reanalysis', get_date_string())


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
class GeneratePED(stage.DatasetStage):
    """
    revert to just using the metamist/CPG-flow Pedigree generation
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'pedigree.ped'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        expected_out = self.expected_outputs(dataset)
        pedigree = dataset.write_ped_file(out_path=expected_out)
        logger.info(f'PED file for {dataset.name} written to {pedigree}')

        return self.make_outputs(dataset, data=expected_out)


@stage.stage
class MakeRuntimeConfig(stage.DatasetStage):
    """
    create a config file for this run,
    this new config should include all elements specific to this Project and sequencing_type
    this new unambiguous config file should be used in all jobs
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'config.toml'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # start off with a fresh config dictionary, including generic content
        new_config: dict = {
            'categories': config.config_retrieve(['categories']),
            'GeneratePanelData': config.config_retrieve(['GeneratePanelData']),
            'RunHailFiltering': config.config_retrieve(['RunHailFiltering']),
            'ValidateMOI': config.config_retrieve(['ValidateMOI']),
            'HPOFlagging': config.config_retrieve(['HPOFlagging']),
            'CreateTalosHTML': {},
        }

        # pull the content relevant to this cohort + sequencing type (mandatory in CPG)
        seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
        dataset_conf = config.config_retrieve(['cohorts', dataset.name])
        seq_type_conf = dataset_conf.get(seq_type, {})

        # forbidden genes and forced panels
        new_config['GeneratePanelData'].update(
            {
                'forbidden_genes': dataset_conf.get('forbidden', []),
                'forced_panels': dataset_conf.get('forced_panels', []),
                'blacklist': dataset_conf.get('blacklist', None),
            },
        )

        # optionally, all SG IDs to remove from analysis
        new_config['ValidateMOI']['solved_cases'] = dataset_conf.get('solved_cases', [])

        # adapt to new hyperlink config structure
        if hyperlinks := seq_type_conf.get('hyperlinks'):
            new_config['CreateTalosHTML']['hyperlinks'] = hyperlinks

        if 'external_labels' in seq_type_conf:
            new_config['CreateTalosHTML']['external_labels'] = seq_type_conf['external_labels']

        # add a location for the run history files
        if 'result_history' in seq_type_conf:
            new_config['result_history'] = seq_type_conf['result_history']

        expected_outputs = self.expected_outputs(dataset)

        with expected_outputs.open('w') as write_handle:
            toml.dump(new_config, write_handle)

        return self.make_outputs(target=dataset, data=expected_outputs)


@stage.stage
class MakePhenopackets(stage.DatasetStage):
    """
    this calls the script which reads phenotype data from metamist
    and generates a phenopacket file (GA4GH compliant)
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'phenopackets.json'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """
        generate a pedigree from metamist
        script to generate an extended pedigree format - additional columns for Ext. ID and HPO terms
        """
        job = set_up_job_with_resources(name=f'MakePhenopackets: {dataset.name}', cpu=1)

        expected_out = self.expected_outputs(dataset)
        query_dataset = dataset.name
        if config.config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
            query_dataset += '-test'

        hpo_file = hail_batch.get_batch().read_input(config.config_retrieve(['GeneratePanelData', 'obo_file']))

        # mandatory argument
        seq_type = config.config_retrieve(['workflow', 'sequencing_type'])

        # insert a little stagger
        job.command(f'sleep {randint(0, 30)}')

        job.command(
            f'MakePhenopackets --dataset {query_dataset} --output {job.output} --type {seq_type} --hpo {hpo_file}',
        )
        hail_batch.get_batch().write_output(job.output, expected_out)
        logger.info(f'Phenopacket file for {dataset.name} going to {expected_out}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(required_stages=[DownloadPanelAppData, MakeRuntimeConfig, MakePhenopackets])
class UnifiedPanelAppParser(stage.DatasetStage):
    """
    Job to parse the PanelApp data, output specific to this Dataset
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'panelapp_data.json'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # create and resource a new job
        job = set_up_job_with_resources(name=f'UnifiedPanelAppParser: {dataset.name}', cpu=1)

        # read in the config made in MakeRuntimeConfig
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        # read in the output from DownloadPanelAppData
        panelapp_download = hail_batch.get_batch().read_input(
            inputs.as_path(target=workflow.get_multicohort(), stage=DownloadPanelAppData),
        )

        local_phenopackets = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=MakePhenopackets))

        hpo_file = hail_batch.get_batch().read_input(config.config_retrieve(['GeneratePanelData', 'obo_file']))

        expected_out = self.expected_outputs(dataset)

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(f"""
        UnifiedPanelAppParser \
            --input {panelapp_download} \
            --output {job.output} \
            --cohort {local_phenopackets} \
            --hpo {hpo_file}
        """)

        hail_batch.get_batch().write_output(job.output, expected_out)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[
        UnifiedPanelAppParser,
        GeneratePED,
        MakeRuntimeConfig,
        SquashMtIntoTarballStage,
    ],
)
class RunHailFiltering(stage.DatasetStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'hail_labelled.vcf.bgz'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # integrate this into the earlier workflow
        input_mt = inputs.as_str(dataset, SquashMtIntoTarballStage)

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        # Talos-prep workflow does a lot of compression, which is nice
        # MTs can vary from <10GB for a small exome, to 170GB for a larger one, Genomes ~200GB
        storage_estimate = tshirt_mt_sizing(
            sequencing_type=config.config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )

        storage = config.config_retrieve(['RunHailFiltering', 'storage', 'small_variants'], storage_estimate)
        cpu: int = config.config_retrieve(['RunHailFiltering', 'cores', 'small_variants'], 8)
        mem: str = config.config_retrieve(['RunHailFiltering', 'memory', 'small_variants'], 'highmem')

        job = set_up_job_with_resources(
            name=f'RunHailFiltering: {dataset.name}',
            storage=f'{storage * 2}Gi',
            cpu=cpu,
            memory=mem,
        )
        job.command('set -eux pipefail')

        panelapp_json = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=UnifiedPanelAppParser))

        # peds can't read cloud paths
        pedigree = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))
        expected_out = self.expected_outputs(dataset)

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

        # read in the massive MT, and unpack it
        localised_mt = hail_batch.get_batch().read_input(input_mt)

        job.command(f'tar -xf {localised_mt} -C $BATCH_TMPDIR && rm {localised_mt}')

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            'RunHailFiltering '
            f'--input "${{BATCH_TMPDIR}}/{dataset.name}.mt" '
            f'--panelapp {panelapp_json} '
            f'--pedigree {pedigree} '
            f'--output {job.output["vcf.bgz"]} '
            f'--clinvar "${{BATCH_TMPDIR}}/clinvarbitration_data/clinvar_decisions.ht" '
            f'--pm5 "${{BATCH_TMPDIR}}/clinvarbitration_data/clinvar_decisions.pm5.ht" '
            f'--checkpoint "${{BATCH_TMPDIR}}/checkpoint.mt" ',
        )
        hail_batch.get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(required_stages=[UnifiedPanelAppParser, GeneratePED, MakeRuntimeConfig])
class RunHailFilteringSV(stage.DatasetStage):
    """
    hail job to filter & label the SV MT
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        if (
            query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
            )
            is not None
        ):
            return dataset.prefix() / get_date_folder() / 'labelled_SVs.vcf.bgz'
        return {}

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # early skip if the stage has nothing to run on
        expected_out = self.expected_outputs(dataset)
        if (
            path_or_none := query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
            )
        ) is None:
            logger.info(f'No SV MT found for {dataset.name}, skipping')
            return self.make_outputs(dataset, data=expected_out)

        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))
        panelapp_json = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=UnifiedPanelAppParser))
        pedigree = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))

        cpu: int = config.config_retrieve(['RunHailFiltering', 'cores', 'sv'], 2)
        job = set_up_job_with_resources(
            name=f'RunHailFilteringSV: {dataset.name}, {path_or_none}',
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
            'RunHailFilteringSV '
            f'--input {annotated_vcf} '
            f'--panelapp {panelapp_json} '
            f'--pedigree {pedigree} '
            f'--mane_json {mane_json} '
            f'--output {job.output["vcf.bgz"]} ',
        )
        hail_batch.get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[
        GeneratePED,
        UnifiedPanelAppParser,
        RunHailFiltering,
        RunHailFilteringSV,
        MakeRuntimeConfig,
    ],
)
class ValidateMOI(stage.DatasetStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'summary_output.json'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # resource consumption here has dropped hugely
        job = set_up_job_with_resources(name=f'ValidateMOI: {dataset.name}', cpu=2.0)

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        panelapp_data = hail_batch.get_batch().read_input(inputs.as_path(dataset, UnifiedPanelAppParser))

        pedigree = hail_batch.get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))
        hail_inputs = inputs.as_path(dataset, RunHailFiltering)

        # either find a SV vcf, or None
        if (
            query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config.config_retrieve(['workflow', 'sequencing_type'])],
            )
            is not None
        ):
            hail_sv_inputs = inputs.as_path(dataset, RunHailFilteringSV)
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
            'ValidateMOI '
            f'--labelled_vcf {labelled_vcf} '
            f'--output {job.output} '
            f'--panelapp {panelapp_data} '
            f'--pedigree {pedigree} '
            f'{sv_vcf_arg}',
        )
        expected_out = self.expected_outputs(dataset)
        hail_batch.get_batch().write_output(job.output, expected_out)
        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[MakeRuntimeConfig, ValidateMOI],
    analysis_type='aip-results',
    analysis_keys=['pheno_annotated', 'pheno_filtered'],
)
class HPOFlagging(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        date_prefix = dataset.prefix() / get_date_folder()
        return {
            'pheno_annotated': date_prefix / 'pheno_annotated_report.json',
            'pheno_filtered': date_prefix / 'pheno_filtered_report.json',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(dataset)

        phenio_db = hail_batch.get_batch().read_input(config.config_retrieve(['HPOFlagging', 'phenio_db']))
        gene_to_phenotype = hail_batch.get_batch().read_input(
            config.config_retrieve(['HPOFlagging', 'gene_to_phenotype']),
        )

        job = set_up_job_with_resources(name=f'HPOFlagging: {dataset.name}', cpu=2, memory='highmem', storage='20Gi')

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        results_json = hail_batch.get_batch().read_input(inputs.as_path(dataset, ValidateMOI))

        mane_json = hail_batch.get_batch().read_input(config.config_retrieve(['references', 'mane_1.4', 'json']))

        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            'HPOFlagging '
            f'--input {results_json} '
            f'--mane_json {mane_json} '
            f'--gen2phen {gene_to_phenotype} '
            f'--phenio {phenio_db} '
            f'--output {job.output} '
            f'--phenout {job.phenout} ',
        )

        hail_batch.get_batch().write_output(job.output, outputs['pheno_annotated'])
        hail_batch.get_batch().write_output(job.phenout, outputs['pheno_filtered'])

        return self.make_outputs(target=dataset, jobs=job, data=outputs)


@stage.stage(required_stages=[HPOFlagging, UnifiedPanelAppParser, MakeRuntimeConfig])
class CreateTalosHtml(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix() / get_date_folder() / 'reports.tar.gz'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        job = set_up_job_with_resources(name=f'CreateTalosHtml: {dataset.name}')

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        results_json = hail_batch.get_batch().read_input(inputs.as_dict(dataset, HPOFlagging)['pheno_annotated'])
        panelapp_data = hail_batch.get_batch().read_input(inputs.as_path(dataset, UnifiedPanelAppParser))
        expected_out = self.expected_outputs(dataset)

        # this will write output files directly to GCP
        job.command(f'export TALOS_CONFIG={runtime_config}')

        # create a new directory for the results
        job.command('mkdir html_outputs')
        job.command('cd html_outputs')
        job.command(f'CreateTalosHTML --input {results_json} --panelapp {panelapp_data} --output summary_output.html')

        # Create a tar'chive here, then use an image with GCloud to copy it up in a bit
        job.command(f'tar -czf {job.output} *')

        hail_batch.get_batch().write_output(job.output, expected_out)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[CreateTalosHtml],
    analysis_type='aip-report',
    tolerate_missing_output=True,
)
class UploadTalosHtml(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix(category='web') / get_date_folder() / 'summary_output.html'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        job = hail_batch.get_batch().new_job(f'UploadTalosHtml: {dataset.name}')
        job.image(config.config_retrieve(['images', 'cpg_hail_gcloud']))
        job.memory('standard')
        job.cpu(1.0)

        input_tarchive = hail_batch.get_batch().read_input(inputs.as_path(dataset, CreateTalosHtml))

        expected_out = self.expected_outputs(dataset)

        hail_batch.authenticate_cloud_credentials_in_job(job)

        # create a new directory for the results
        output_folder = str(dataset.prefix(category='web') / get_date_folder())
        job.command('mkdir html_outputs')
        job.command(f'tar -xf {input_tarchive} -C html_outputs')
        job.command('cd html_outputs')
        job.command(f'gcloud storage cp -r * {output_folder}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage.stage(
    required_stages=[CreateTalosHtml],
    tolerate_missing_output=True,
    forced=True,
)
class UploadGenericTalosHtml(stage.DatasetStage):
    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        return dataset.prefix(category='web') / 'talos_static' / 'summary_output.html'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        job = hail_batch.get_batch().new_job(f'UploadGenericTalosHtml: {dataset.name}')
        job.image(config.config_retrieve(['images', 'cpg_hail_gcloud']))
        job.memory('standard')
        job.cpu(1.0)

        input_tarchive = hail_batch.get_batch().read_input(inputs.as_path(dataset, CreateTalosHtml))

        output = self.expected_outputs(dataset)

        hail_batch.authenticate_cloud_credentials_in_job(job)

        # unpack results into a local directory, then copy to GCP using gcloud rsync
        output_folder = dataset.prefix(category='web') / 'talos_static'
        job.command('mkdir html_outputs')
        job.command(f'tar -xf {input_tarchive} -C html_outputs')
        job.command(
            f'gcloud storage rsync --recursive html_outputs --delete-unmatched-destination-objects {output_folder!s}',
        )

        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(
    required_stages=[ValidateMOI, MakeRuntimeConfig],
    analysis_keys=['seqr_file', 'seqr_pheno_file'],
    analysis_type='custom',
    tolerate_missing_output=True,
)
class MinimiseOutputForSeqr(stage.DatasetStage):
    """
    takes the results file from Seqr and produces a minimised form for Seqr ingestion
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, Path]:
        analysis_folder_prefix = dataset.prefix(category='analysis') / 'seqr_files'
        return {
            'seqr_file': analysis_folder_prefix / f'{get_date_folder()}_seqr.json',
            'seqr_pheno_file': analysis_folder_prefix / f'{get_date_folder()}_seqr_pheno.json',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        # pull out the config section relevant to this datatype & cohort
        # if it doesn't exist for this sequencing type, fail gracefully
        seq_type = config.config_retrieve(['workflow', 'sequencing_type'])
        try:
            seqr_lookup = config.config_retrieve(['cohorts', dataset.name, seq_type, 'seqr_lookup'])
        except config.ConfigError:
            logger.warning(f'No Seqr lookup file for {dataset.name} {seq_type}')
            return self.make_outputs(dataset, skipped=True)

        input_localised = hail_batch.get_batch().read_input(inputs.as_path(dataset, ValidateMOI))

        # create a job to run the minimisation script
        job = set_up_job_with_resources(name=f'MinimiseOutputForSeqr: {dataset.name}', memory='lowmem')

        # use the new config file
        runtime_config = hail_batch.get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        lookup_in_batch = hail_batch.get_batch().read_input(seqr_lookup)
        job.command(f'export TALOS_CONFIG={runtime_config}')
        job.command(
            'MinimiseOutputForSeqr '
            f'--input {input_localised} '
            f'--output {job.out_json} '
            f'--pheno {job.pheno_json} '
            f'--external_map {lookup_in_batch}',
        )

        # write the results out
        expected_out = self.expected_outputs(dataset)
        hail_batch.get_batch().write_output(job.out_json, expected_out['seqr_file'])
        hail_batch.get_batch().write_output(job.pheno_json, expected_out['seqr_pheno_file'])
        return self.make_outputs(dataset, data=expected_out, jobs=job)
