#!/usr/bin/env python3

"""
takes one or more files, checks they were created, register in metamist
"""
import logging
import sys
from os.path import join
from pathlib import Path

import click
from cpg_utils import to_path
from cpg_utils.config import get_config
from metamist.apis import AnalysisApi
from metamist.model.analysis import Analysis
from metamist.model.analysis_status import AnalysisStatus
from peddy import Ped


def register_html(file_path: str, samples: list[str]):
    """
    Takes the output HTML from this analysis and registers it in Metamist

    Args:
        file_path (str): result file to register
        samples (list[str]): relevant samples
    """

    # Create object Meta - Exomes/genomes, Singletons/not, proxied html path
    report_meta: dict[str, bool | str] = {
        'is_exome': bool('exome' in get_config()['workflow']['output_prefix']),
        'is_singleton': bool('singleton' in file_path),
    }

    # add HTML-specific elements
    if file_path.endswith('html'):
        analysis_type = 'aip-report'
        report_meta['display_url'] = join(
            get_config()['storage']['default']['web_url'],
            get_config()['workflow']['output_prefix'],
            Path(file_path).name,
        )
    else:
        analysis_type = 'aip-results'

    AnalysisApi().create_analysis(
        project=get_config()['workflow']['dataset'],
        analysis=Analysis(
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
    main()
