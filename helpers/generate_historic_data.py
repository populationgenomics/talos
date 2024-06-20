#! /usr/bin/env python3

"""
This is a script which fills a gap in the main pipeline

We are unable to easily generate data in the expected 'historical' format
from a run after that run's completion. If the 'historic' file was not
created, or was deleted for some reason, this would reduce the efficacy
of subsequent runs.
This script is intended to be run after the main pipeline has completed,
and will generate the historical prior file from the full results JSON.
"""

from argparse import ArgumentParser

from cpg_utils.config import get_config

from talos.models import ResultData
from talos.utils import generate_fresh_latest_results, read_json_from_path


def main(input_results: str):
    """
    read in the target full results JSON
    call the generate_fresh_latest_results function
    this generates new data and saves to a reasonable location

    Args:
        input_results ():
    """

    input_data: ResultData = read_json_from_path(input_results, return_model=ResultData)  # type: ignore
    generate_fresh_latest_results(input_data, prefix='', dataset=get_config()['workflow']['dataset'])


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', help='GCP path to full results JSON', required=True)
    args = parser.parse_args()
    main(input_results=args.i)
