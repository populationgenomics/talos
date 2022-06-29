"""
A method for constructing a dictionary of CPG sample ID to Seqr Family ID

Input:
use the developer console in a browser
load a study page in seqr with the network tab open
save the response of the `overview` API call as JSON

Seqr doesn't directly expose an API to access this data... yet
Plan A is contact Seqr directly, so these methods would still be useful
Plan B would be to read similar content from the Sample Metadata DB
Plan ... whatever, would be to do this
"""

import json
import os
from argparse import ArgumentParser

from reanalysis.utils import read_json_from_path


def get_seqr_details(seqr_meta: str) -> dict[str, dict[str, str]]:
    """
    this would be substituted for a GET call
    :return: dict of the semi-processed output
    """

    if not os.path.exists(seqr_meta):
        raise Exception(f'Input file {seqr_meta} inaccessible')

    details_dict = read_json_from_path(seqr_meta)

    # map CPG ID to individual GUID
    return {
        sample['sampleId']: sample['familyGuid']
        for seqr_sample_id, sample in details_dict['samplesByGuid'].items()
    }


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-i', type=str, help='input JSON file')
    parser.add_argument('-o', type=str, help='JSON path to write output on')
    args = parser.parse_args()

    json_digest = get_seqr_details(args.i)
    with open(args.o, 'w', encoding='utf-8') as handle:
        json.dump(json_digest, handle, indent=4, sort_keys=True)
