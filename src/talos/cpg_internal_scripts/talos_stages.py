"""
This is a central script for the Talos process, implemented at the CPG, using the cpg-flow workflow framework
"""

import toml
from functools import cache, lru_cache
from os.path import join
from random import randint
from typing import TYPE_CHECKING

from cpg_flow.stage import DatasetStage, stage
from cpg_flow.utils import get_logger
from cpg_utils import Path
from cpg_utils.config import ConfigError, config_retrieve
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, get_batch

from metamist.graphql import gql, query
from talos.utils import get_granular_date

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets.dataset import Dataset
    from hailtop.batch.job import Job


METAMIST_ANALYSIS_QUERY = gql(
    """
    query MyQuery($dataset: String!, $type: String!) {
        project(name: $dataset) {
            analyses(active: {eq: true}, type: {eq: $type}, status: {eq: COMPLETED}) {
                output
                timestampCompleted
                meta
            }
        }
    }
""",
)

TALOS_PREP_TYPE = 'talos_prep'

# these are the values used in the metamist query
CLINVARBITRATION_PROJECT = 'fewgenomes'
CLINVARBITRATION_TYPE = 'clinvarbitration'

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
) -> 'Job':
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

    job = get_batch().new_job(name=name)
    if image:
        job.image(image)
    else:
        job.image(config_retrieve(['workflow', 'driver_image']))
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


@lru_cache(maxsize=1)
def get_date_string() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses

    Returns:
        either an override in config, or the default (today, YYYY-MM-DD)
    """
    return config_retrieve(['workflow', 'date_folder_override'], get_granular_date())


@lru_cache(1)
def get_date_folder() -> str:
    """
    allows override of the date folder to continue/re-run previous analyses
    Returns:
        either an override in config, or the default "reanalysis/(today, YYYY-MM-DD)"
    """
    return join('reanalysis', get_date_string())


@cache
def query_for_latest_analysis(
    dataset: str,
    analysis_type: str,
    sequencing_type: str = 'any',
) -> str | None:
    """
    Query for the latest analysis object of a given type in the requested project
    Analysis entries for Talos all have unique types, so we can use this generic query method

    Args:
        dataset (str):         project to query for
        analysis_type (str):   analysis type to query for - rd_combiner writes MTs to metamist as 'matrixtable',
                               seqr_loader used 'custom': using a config entry we can decide which type to use
        sequencing_type (str): optional, if set, only return entries with meta.sequencing_type == this
    Returns:
        str, the path to the latest object for the given type, or log a warning and return None
    """

    # swapping to a string we can freely modify
    query_dataset = dataset
    if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
        query_dataset += '-test'

    get_logger().info(f'Querying for {analysis_type} in {query_dataset}')

    result = query(METAMIST_ANALYSIS_QUERY, variables={'dataset': query_dataset, 'type': analysis_type})

    # get all the relevant entries, and bin by date
    analysis_by_date = {}
    for analysis in result['project']['analyses']:
        if analysis['output'] and (sequencing_type in {'all', analysis['meta'].get('sequencing_type')}):
            analysis_by_date[analysis['timestampCompleted']] = analysis['output']

    if not analysis_by_date:
        get_logger().warning(f'No Analysis Entries found for dataset {query_dataset}')
        return None

    # return the latest, determined by a sort on timestamp
    # 2023-10-10... > 2023-10-09..., so sort on strings
    return analysis_by_date[sorted(analysis_by_date)[-1]]


@stage
class GeneratePED(DatasetStage):
    """
    revert to just using the metamist/CPG-flow Pedigree generation
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'pedigree.ped'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        expected_out = self.expected_outputs(dataset)
        pedigree = dataset.write_ped_file(out_path=expected_out)
        get_logger().info(f'PED file for {dataset.name} written to {pedigree}')

        return self.make_outputs(dataset, data=expected_out)


@stage
class MakeRuntimeConfig(DatasetStage):
    """
    create a config file for this run,
    this new config should include all elements specific to this Project and sequencing_type
    this new unambiguous config file should be used in all jobs
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'config.toml'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # start off with a fresh config dictionary, including generic content
        new_config: dict = {
            'categories': config_retrieve(['categories']),
            'GeneratePanelData': config_retrieve(['GeneratePanelData']),
            'RunHailFiltering': config_retrieve(['RunHailFiltering']),
            'ValidateMOI': config_retrieve(['ValidateMOI']),
            'HPOFlagging': config_retrieve(['HPOFlagging']),
            'CreateTalosHTML': {},
        }

        # pull the content relevant to this cohort + sequencing type (mandatory in CPG)
        seq_type = config_retrieve(['workflow', 'sequencing_type'])
        dataset_conf = config_retrieve(['cohorts', dataset.name])
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


@stage
class MakePhenopackets(DatasetStage):
    """
    this calls the script which reads phenotype data from metamist
    and generates a phenopacket file (GA4GH compliant)
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'phenopackets.json'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        """
        generate a pedigree from metamist
        script to generate an extended pedigree format - additional columns for Ext. ID and HPO terms
        """
        job = set_up_job_with_resources(name=f'MakePhenopackets: {dataset.name}', cpu=1)

        expected_out = self.expected_outputs(dataset)
        query_dataset = dataset.name
        if config_retrieve(['workflow', 'access_level']) == 'test' and 'test' not in query_dataset:
            query_dataset += '-test'

        hpo_file = get_batch().read_input(config_retrieve(['GeneratePanelData', 'obo_file']))

        # mandatory argument
        seq_type = config_retrieve(['workflow', 'sequencing_type'])

        # insert a little stagger
        job.command(f'sleep {randint(0, 30)}')

        job.command(
            f'MakePhenopackets --dataset {query_dataset} --output {job.output} --type {seq_type} --hpo {hpo_file}',
        )
        get_batch().write_output(job.output, expected_out)
        get_logger().info(f'Phenopacket file for {dataset.name} going to {expected_out}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[MakePhenopackets, MakeRuntimeConfig])
class GeneratePanelData(DatasetStage):
    """
    PythonJob to find HPO-matched panels
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        """
        only one output, the panel data
        """
        return dataset.prefix() / get_date_folder() / 'hpo_panel_data.json'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # create and resource a new job
        job = set_up_job_with_resources(name=f'GeneratePanelData: {dataset.name}', cpu=1)

        # read in the config made in MakeRuntimeConfig
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        expected_out = self.expected_outputs(dataset)

        # read in a HPO file
        hpo_file = get_batch().read_input(config_retrieve(['GeneratePanelData', 'obo_file']))

        # read in the output from MakePhenopackets
        local_phenopacket = get_batch().read_input(inputs.as_path(target=dataset, stage=MakePhenopackets))

        job.command(f'export TALOS_CONFIG={runtime_config}')

        # insert a little stagger so PanelApp isn't completely swamped
        job.command(f'sleep {randint(0, 30)}')
        job.command(f'GeneratePanelData --input {local_phenopacket} --output {job.output} --hpo {hpo_file}')
        get_batch().write_output(job.output, expected_out)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[GeneratePanelData, MakeRuntimeConfig])
class QueryPanelapp(DatasetStage):
    """
    query PanelApp for up-to-date gene lists
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'panelapp_data.json'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        job = set_up_job_with_resources(name=f'QueryPanelApp: {dataset.name}', cpu=1)

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        hpo_panel_json = get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePanelData))
        expected_out = self.expected_outputs(dataset)

        job.command(f'export TALOS_CONFIG={runtime_config}')
        # insert a little stagger so PanelApp isn't completely swamped
        job.command(f'sleep {randint(20, 300)}')
        job.command(f'QueryPanelapp --input {hpo_panel_json} --output {job.output}')
        get_batch().write_output(job.output, expected_out)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED, MakeRuntimeConfig])
class RunHailFiltering(DatasetStage):
    """
    hail job to filter & label the MT
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'hail_labelled.vcf.bgz'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # if it's not set in config, try and get it through metamist
        if not (input_mt := config_retrieve(['workflow', 'matrix_table'], None)):
            input_mt = query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=TALOS_PREP_TYPE,
            )

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        # Talos-prep workflow does a lot of compression, which is nice
        # MTs can vary from <10GB for a small exome, to 170GB for a larger one, Genomes ~200GB
        storage_estimate = tshirt_mt_sizing(
            sequencing_type=config_retrieve(['workflow', 'sequencing_type']),
            cohort_size=len(dataset.get_sequencing_group_ids()),
        )

        storage = config_retrieve(['RunHailFiltering', 'storage', 'small_variants'], storage_estimate)
        cpu: int = config_retrieve(['RunHailFiltering', 'cores', 'small_variants'], 8)
        mem: str = config_retrieve(['RunHailFiltering', 'memory', 'small_variants'], 'highmem')

        job = set_up_job_with_resources(
            name=f'RunHailFiltering: {dataset.name}',
            storage=f'{storage * 2}Gi',
            cpu=cpu,
            memory=mem,
        )
        job.command('set -eux pipefail')

        panelapp_json = get_batch().read_input(inputs.as_path(target=dataset, stage=QueryPanelapp))

        # peds can't read cloud paths
        pedigree = get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))
        expected_out = self.expected_outputs(dataset)

        # copy vcf & index out manually
        job.declare_resource_group(
            output={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # find the clinvar table, localise, and expand
        clinvar_tar = query_for_latest_analysis(
            dataset=CLINVARBITRATION_PROJECT,
            analysis_type=CLINVARBITRATION_TYPE,
        )

        if clinvar_tar is None:
            raise ValueError('No ClinVar data found')

        job.command(f'tar -xzf {get_batch().read_input(clinvar_tar)} -C $BATCH_TMPDIR')

        # read in the massive MT, and unpack it
        localised_mt = get_batch().read_input(input_mt)

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
        get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(required_stages=[QueryPanelapp, GeneratePED, MakeRuntimeConfig])
class RunHailFilteringSV(DatasetStage):
    """
    hail job to filter & label the SV MT
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        if (
            query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config_retrieve(['workflow', 'sequencing_type'])],
            )
            is not None
        ):
            return dataset.prefix() / get_date_folder() / 'labelled_SVs.vcf.bgz'
        return {}

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # early skip if the stage has nothing to run on
        expected_out = self.expected_outputs(dataset)
        if (
            path_or_none := query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config_retrieve(['workflow', 'sequencing_type'])],
            )
        ) is None:
            get_logger().info(f'No SV MT found for {dataset.name}, skipping')
            return self.make_outputs(dataset, data=expected_out)

        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))
        panelapp_json = get_batch().read_input(inputs.as_path(target=dataset, stage=QueryPanelapp))
        pedigree = get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))

        cpu: int = config_retrieve(['RunHailFiltering', 'cores', 'sv'], 2)
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
        mane_json = get_batch().read_input(config_retrieve(['references', 'mane_1.4', 'json']))

        # copy the VCF in
        annotated_vcf = get_batch().read_input_group(
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
        get_batch().write_output(job.output, str(expected_out).removesuffix('.vcf.bgz'))

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[
        GeneratePED,
        GeneratePanelData,
        QueryPanelapp,
        RunHailFiltering,
        RunHailFilteringSV,
        MakeRuntimeConfig,
    ],
)
class ValidateMOI(DatasetStage):
    """
    run the labelled VCF -> results JSON stage
    """

    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'summary_output.json'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # resource consumption here has dropped hugely
        job = set_up_job_with_resources(name=f'ValidateMOI: {dataset.name}', cpu=2.0)

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        hpo_panels = get_batch().read_input(inputs.as_path(dataset, GeneratePanelData))
        pedigree = get_batch().read_input(inputs.as_path(target=dataset, stage=GeneratePED))
        hail_inputs = inputs.as_path(dataset, RunHailFiltering)

        # either find a SV vcf, or None
        if (
            query_for_latest_analysis(
                dataset=dataset.name,
                analysis_type=SV_ANALYSIS_TYPES[config_retrieve(['workflow', 'sequencing_type'])],
            )
            is not None
        ):
            hail_sv_inputs = inputs.as_path(dataset, RunHailFilteringSV)
            labelled_sv_vcf = get_batch().read_input_group(
                **{'vcf.bgz': hail_sv_inputs, 'vcf.bgz.tbi': f'{hail_sv_inputs}.tbi'},
            )['vcf.bgz']
            sv_vcf_arg = f'--labelled_sv {labelled_sv_vcf} '
        else:
            sv_vcf_arg = ''

        panel_input = get_batch().read_input(inputs.as_path(dataset, QueryPanelapp))
        labelled_vcf = get_batch().read_input_group(
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
            f'--panelapp {panel_input} '
            f'--pedigree {pedigree} '
            f'--participant_panels {hpo_panels} '
            f'{sv_vcf_arg}',
        )
        expected_out = self.expected_outputs(dataset)
        get_batch().write_output(job.output, expected_out)
        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[MakeRuntimeConfig, ValidateMOI],
    analysis_type='aip-results',
    analysis_keys=['pheno_annotated', 'pheno_filtered'],
)
class HPOFlagging(DatasetStage):
    def expected_outputs(self, dataset: 'Dataset') -> dict[str, Path]:
        date_prefix = dataset.prefix() / get_date_folder()
        return {
            'pheno_annotated': date_prefix / 'pheno_annotated_report.json',
            'pheno_filtered': date_prefix / 'pheno_filtered_report.json',
        }

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(dataset)

        phenio_db = get_batch().read_input(config_retrieve(['HPOFlagging', 'phenio_db']))
        gene_to_phenotype = get_batch().read_input(config_retrieve(['HPOFlagging', 'gene_to_phenotype']))

        job = set_up_job_with_resources(name=f'HPOFlagging: {dataset.name}', cpu=2, memory='highmem', storage='20Gi')

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        results_json = get_batch().read_input(inputs.as_path(dataset, ValidateMOI))

        mane_json = get_batch().read_input(config_retrieve(['references', 'mane_1.4', 'json']))

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

        get_batch().write_output(job.output, outputs['pheno_annotated'])
        get_batch().write_output(job.phenout, outputs['pheno_filtered'])

        return self.make_outputs(target=dataset, jobs=job, data=outputs)


@stage(required_stages=[HPOFlagging, QueryPanelapp, MakeRuntimeConfig])
class CreateTalosHtml(DatasetStage):
    def expected_outputs(self, dataset: 'Dataset') -> Path:
        return dataset.prefix() / get_date_folder() / 'reports.tar.gz'

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        job = set_up_job_with_resources(name=f'CreateTalosHtml: {dataset.name}')

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        results_json = get_batch().read_input(inputs.as_dict(dataset, HPOFlagging)['pheno_annotated'])
        panel_input = get_batch().read_input(inputs.as_path(dataset, QueryPanelapp))
        expected_out = self.expected_outputs(dataset)

        # this will write output files directly to GCP
        job.command(f'export TALOS_CONFIG={runtime_config}')

        # create a new directory for the results
        job.command('mkdir html_outputs')
        job.command('cd html_outputs')
        job.command(f'CreateTalosHTML --input {results_json} --panelapp {panel_input} --output summary_output.html')

        # Create a tar'chive here, then use an image with GCloud to copy it up in a bit
        job.command(f'tar -czf {job.output} *')

        get_batch().write_output(job.output, expected_out)

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[CreateTalosHtml],
    analysis_type='aip-report',
    analysis_keys=['results_html'],
    tolerate_missing_output=True,
)
class UploadTalosHtml(DatasetStage):
    def expected_outputs(self, dataset: 'Dataset') -> dict[str, str | Path]:
        date_folder_prefix = dataset.prefix(category='web') / get_date_folder()
        return {
            'results_html': date_folder_prefix / 'summary_output.html',
            'folder': str(date_folder_prefix),
        }

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        job = get_batch().new_job(f'UploadTalosHtml: {dataset.name}')

        # TODO this needs a basic gcloud image
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.memory('standard')
        job.cpu(1.0)

        input_tarchive = get_batch().read_input(inputs.as_path(dataset, CreateTalosHtml))

        expected_out = self.expected_outputs(dataset)

        authenticate_cloud_credentials_in_job(job)

        # create a new directory for the results
        job.command('mkdir html_outputs')
        job.command(f'tar -xf {input_tarchive} -C html_outputs')
        job.command('cd html_outputs')
        job.command(f'gcloud storage cp -r * {expected_out["folder"]}')

        return self.make_outputs(dataset, data=expected_out, jobs=job)


@stage(
    required_stages=[ValidateMOI, MakeRuntimeConfig],
    analysis_keys=['seqr_file', 'seqr_pheno_file'],
    analysis_type='custom',
    tolerate_missing_output=True,
)
class MinimiseOutputForSeqr(DatasetStage):
    """
    takes the results file from Seqr and produces a minimised form for Seqr ingestion
    """

    def expected_outputs(self, dataset: 'Dataset') -> dict[str, Path]:
        analysis_folder_prefix = dataset.prefix(category='analysis') / 'seqr_files'
        return {
            'seqr_file': analysis_folder_prefix / f'{get_date_folder()}_seqr.json',
            'seqr_pheno_file': analysis_folder_prefix / f'{get_date_folder()}_seqr_pheno.json',
        }

    def queue_jobs(self, dataset: 'Dataset', inputs: 'StageInput') -> 'StageOutput':
        # pull out the config section relevant to this datatype & cohort
        # if it doesn't exist for this sequencing type, fail gracefully
        seq_type = config_retrieve(['workflow', 'sequencing_type'])
        try:
            seqr_lookup = config_retrieve(['cohorts', dataset.name, seq_type, 'seqr_lookup'])
        except ConfigError:
            get_logger().warning(f'No Seqr lookup file for {dataset.name} {seq_type}')
            return self.make_outputs(dataset, skipped=True)

        input_localised = get_batch().read_input(inputs.as_path(dataset, ValidateMOI))

        # create a job to run the minimisation script
        job = set_up_job_with_resources(name=f'MinimiseOutputForSeqr: {dataset.name}', memory='lowmem')

        # use the new config file
        runtime_config = get_batch().read_input(inputs.as_path(dataset, MakeRuntimeConfig))

        lookup_in_batch = get_batch().read_input(seqr_lookup)
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
        get_batch().write_output(job.out_json, expected_out['seqr_file'])
        get_batch().write_output(job.pheno_json, expected_out['seqr_pheno_file'])
        return self.make_outputs(dataset, data=expected_out, jobs=job)
