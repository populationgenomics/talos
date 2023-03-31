#!/usr/bin/env python3

"""
takes one or more files, checks they were created, register in metamist
"""
import logging
from os.path import join
from pathlib import Path
import sys

import click
from peddy import Ped

from cpg_utils import to_path
from cpg_utils.config import get_config

from sample_metadata.apis import AnalysisApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_model import AnalysisModel
from sample_metadata.model.analysis_status import AnalysisStatus
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.model.analysis_update_model import AnalysisUpdateModel


def inactivate_old_entry(name: str, metadata: dict[str, str], anal_type: AnalysisType):
    """
    for this file, metadata, and analysis type:
        - check for any pre-existing reanalysis entries in metamist
        - confirm an exact match:
            - same exome/genome status
            - same filename
            - same analysis type
        - update previous record to inactive

    Args:
        name (str): filename
        metadata (dict): dictionary of parameters
        anal_type (AnalysisType):  Custom/Web
    """

    # find any previous  AIP-specific AnalysisEntries... Update to active=False
    a_query_model = AnalysisQueryModel(
        projects=[get_config()['workflow']['dataset']], type=anal_type
    )
    for analysis in AnalysisApi().query_analyses(analysis_query_model=a_query_model):

        # only look for reanalysis entries
        if 'reanalysis' not in analysis['output']:
            continue

        # skip over reports that don't match this subtype
        for key, value in metadata.items():
            if analysis['meta'][key] != value:
                continue

        # check that the name is the same, check its active, then kill it
        output_path = Path(analysis['output'])
        if name.endswith(output_path.name) and analysis['active']:

            # update this entry to inactivate it (superseded)
            AnalysisApi().update_analysis_status(
                analysis_id=analysis['id'],
                analysis_update_model=AnalysisUpdateModel(
                    status=AnalysisStatus('completed'), active=False
                ),
            )


def register_html(file_path: str, samples: list[str]):
    """
    Takes the output HTML from this analysis and registers it in
    Metamist. Deprecates any existing HTML results

    Args:
        file_path (str): result file to register
        samples (list[str]): relevant samples
    """

    # Create object Meta - Exomes/genomes, Singletons/not, proxied html path
    report_meta = {
        'is_exome': bool('exome' in get_config()['workflow']['output_prefix']),
        'is_singleton': bool('singleton' in file_path),
    }

    analysis_type = (
        AnalysisType('web') if 'html' in file_path else AnalysisType('custom')
    )

    # find any previous versions of this AnalysisEntry, and deactivate
    inactivate_old_entry(name=file_path, metadata=report_meta, anal_type=analysis_type)

    # add HTML-specific elements
    if file_path.endswith('html'):

        web_template = get_config()['storage']['default']['web_url']
        file_name = Path(file_path).name
        display_url = join(
            web_template, get_config()['workflow']['output_prefix'], file_name
        )
        report_meta['display_url'] = display_url

    AnalysisApi().create_new_analysis(
        project=get_config()['workflow']['dataset'],
        analysis_model=AnalysisModel(
            sample_ids=samples,
            type=analysis_type,
            status=AnalysisStatus('completed'),
            output=file_path,
            meta=report_meta,
            active=True,
        ),
    )


def get_samples_from_pedigree(pedfile: str) -> list[str]:
    """
    parses the pedigree used in this analysis to find all sample IDs

    nb. prone to error - really we want an intersection between
        the pedigree and the actual VCF/MT
        e.g. this could result in the registration of results for
        samples who were in the pedigree and not the joint call

    Args:
        pedfile (str): path to the pedigree file

    Returns:
        list of all collected sample IDs
    """

    pedigree = Ped(pedfile)

    # yank out all the sample IDs used in this analysis
    samples = sorted(s.sample_id for s in pedigree.samples())

    return samples


@click.command()
@click.option('--pedigree', help='pedigree file to pull samples from')
@click.argument('files', nargs=-1)
def main(pedigree: str, files: list[str]):
    """
    Takes output file paths from this analysis and registers

    Args:
        files (list[str]): list of files to operate on
        pedigree (str): Pedigree file to pull samples from
    """

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )

    if get_config()['workflow']['access_level'] == 'test':
        logging.warning('Refusing to register files generated in test')
        # never update metamist in test mode - no permission
        sys.exit(0)

    samples = get_samples_from_pedigree(pedigree)

    logging.info(f'Metamist entries are lined to {len(samples)} samples')

    for file in files:
        if to_path(file).exists():
            register_html(file_path=file, samples=samples)
        else:
            logging.error(f'Could not see file {file}, will not be registered')
    logging.info('Completed file registration')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
