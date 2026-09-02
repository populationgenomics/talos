import datetime
from os.path import join

from cpg_flow import stage, targets
from cpg_utils import Path, config, hail_batch, to_path

from talos.cpg_internal_scripts.cpg_flow_utils import set_up_job_with_resources

THIS_MONTH = datetime.datetime.now().strftime('%y-%m')  # noqa: DTZ005


@stage.stage(analysis_type='panelapp')
class DownloadPanelAppData(stage.MultiCohortStage):
    """
    runs a single instance of the stage which downloads the whole of PanelApp into a cached file
    """

    def expected_outputs(self, _multicohort: targets.MultiCohort) -> Path:
        return to_path(
            join(
                config.config_retrieve(['storage', 'common', 'analysis']),
                'panelapp_monthly',
                f'panelapp_data_{THIS_MONTH}.json',
            ),
        )

    def queue_jobs(
        self,
        multicohort: targets.MultiCohort,
        _inputs: stage.StageInput,
    ) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)
        job = set_up_job_with_resources(name='DownloadPanelAppData', cpu=1)
        job.command(f'python -m talos.download_panelapp --output {job.output}')

        hail_batch.get_batch().write_output(job.output, output)

        return self.make_outputs(target=multicohort, data=output, jobs=job)
