"""
compares a Seqr dump against results of a previous AIP run

See relevant documentation for a description of the algorithm used
"""

from collections import defaultdict
from enum import Enum
import logging
import re
import sys
from typing import Dict, List, Tuple

from argparse import ArgumentParser
import hail as hl

from reanalysis.interpretation_runner import read_json_dict_from_path


SAMPLE_NUM_RE = re.compile(r'sample_[0-9]+')
SAMPLE_ALT_TEMPLATE = 'num_alt_alleles_{}'


class Confidence(Enum):
    """
    enumeration for possible confidence values
    from AIP we get CERTAIN, as these are found
    from seqr dump, we use tags to determine confidence
        - where confidence means 'how confident are we that
            AIP should be capturing this result'
    """

    EXPECTED = 'AIP training: Expected'
    POSSIBLE = 'AIP training: Possible'
    UNLIKELY = 'AIP training: Unlikely'

    def __lt__(self, other):
        return self.value < other.value


class CommonFormatResult:
    """
    a common representation of variant details
    """

    def __init__(
        self, chrom: str, pos: int, ref: str, alt: str, confidence: List[Confidence]
    ):
        """

        :param chrom:
        :param pos:
        :param ref:
        :param alt:
        :param confidence:
        """
        self.chr: str = chrom
        self.pos: int = pos
        self.ref: str = ref
        self.alt: str = alt
        self.confidence: List[Confidence] = confidence

    def get_cyvcf2_pos(self):
        """
        get variant coordinates in string format
        :return:
        """
        return f'{self.chr}:{self.pos}-{self.pos}'

    def get_hail_pos(self) -> Tuple[hl.Locus, List[str]]:
        """
        get relevant values for finding a variant row within a MT
        :return:
        """
        return hl.Locus(contig=f'chr{self.chr}', position=self.pos), [
            self.ref,
            self.alt,
        ]

    def __repr__(self):
        return f'{self.chr}:{self.pos}_{self.ref}>{self.alt} - {self.confidence}'

    def __eq__(self, other):
        return (
            self.chr == other.chr
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
        )

    def __hash__(self):
        """
        hash function
        :return:
        """
        return hash(
            repr(self) + ''.join(list(map(lambda x: x.value, sorted(self.confidence))))
        )


CommonDict = Dict[str, List[CommonFormatResult]]


def common_format_from_results(aip_results: str) -> CommonDict:
    """
    Parses the JSON

    :param aip_results:
    :return:
    """
    results_dict = read_json_dict_from_path(aip_results)
    sample_dict: CommonDict = defaultdict(list)

    # collect all per-sample results into a separate index
    for sample, variants in results_dict.items():

        for var in variants.values():
            coords = var['var_data']['coords']
            sample_dict[sample].append(
                CommonFormatResult(
                    coords['chrom'],
                    int(coords['pos']),
                    coords['ref'],
                    coords['alt'],
                    [Confidence.EXPECTED],
                )
            )

    return sample_dict


def main(results: str):
    """

    :param results:
    :param seqr:
    :return:
    """

    # normalise data formats from AIP result file
    _result_dict = common_format_from_results(aip_results=results)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('--results')
    args = parser.parse_args()
    main(results=args.results)
