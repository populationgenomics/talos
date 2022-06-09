"""
compares a Seqr file against results of a previous AIP run

The Seqr file is a TSV of all variants flagged as interesting
This is designed to recognise flags in the format 'AIP training: Confidence'

See relevant documentation for a description of the algorithm used
"""

from collections import defaultdict
from csv import DictReader
from enum import Enum
import logging
import re
import sys
from typing import Any

from argparse import ArgumentParser
from cloudpathlib import AnyPath
import hail as hl
from peddy import Ped

from reanalysis.utils import read_json_from_path


SAMPLE_NUM_RE = re.compile(r'sample_[0-9]+')
SAMPLE_ALT_TEMPLATE = 'num_alt_alleles_{}'
VALID_VALUES = [
    'AIP training: Expected',
    'AIP training: Possible',
    'AIP training: Unlikely',
]


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
        self, chrom: str, pos: int, ref: str, alt: str, confidence: list[Confidence]
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
        self.confidence: list[Confidence] = confidence

    def get_cyvcf2_pos(self):
        """
        get variant coordinates in string format
        :return:
        """
        return f'{self.chr}:{self.pos}-{self.pos}'

    def get_hail_pos(self) -> tuple[hl.Locus, list[str]]:
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


CommonDict = dict[str, list[CommonFormatResult]]


def common_format_from_results(results_dict: dict[str, Any]) -> CommonDict:
    """
    Parses the JSON

    :param results_dict:
    :return:
    """
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


def common_format_from_seqr(seqr: str, affected: list[str]) -> CommonDict:
    """
    identify the most likely proband for each row, and create a variant object for them

    :param seqr:
    :param affected:
    :return:
    """

    assert seqr.endswith('.tsv'), f'process requires a TSV format export'
    sample_dict: CommonDict = defaultdict(list)

    with AnyPath(seqr).open() as handle:
        seqr_parser = DictReader(handle, delimiter='\t')

        # Each Sample ID is under a separate column heading, e.g. sample_1
        # this is instead of proband/mother/father; no mandatory family structure
        sample_cols = [
            sample_field
            for sample_field in seqr_parser.fieldnames
            if SAMPLE_NUM_RE.match(sample_field)
        ]

        for entry in seqr_parser:

            # get all valid tags
            tags = [
                Confidence(tag)
                for tag in entry['tags'].split('|')
                if tag in VALID_VALUES
            ]

            # no relevant tags, not interested...
            if len(tags) == 0:
                continue

            # create a variant object
            variant_obj = CommonFormatResult(
                entry['chrom'],
                int(entry['pos']),
                entry['ref'],
                entry['alt'],
                confidence=tags,
            )

            # seqr has no notion of 'proband', so add for each affected
            # see README for discussion
            for index, value in enumerate(sample_cols, 1):
                if (
                    entry[value] in affected
                    and int(entry[SAMPLE_ALT_TEMPLATE.format(index)]) > 0
                ):
                    sample_dict[entry[value]].append(variant_obj)

    logging.info(f'Variants from Seqr digest: {sample_dict}')

    return sample_dict


def find_missing(aip_results: CommonDict, seqr_results: CommonDict) -> CommonDict:
    """
    1 way comparison - find all results present in Seqr, and missing from AIP
    This is a check for False Negatives-
    :param aip_results:
    :param seqr_results:
    :return:
    """
    discrepancies: CommonDict = defaultdict(list)

    # get all common samples between the two data sets
    seqr_samples = set(seqr_results.keys())
    aip_samples = set(aip_results.keys())

    common_samples = aip_samples.intersection(seqr_samples)

    missing_samples = seqr_samples - common_samples
    if len(missing_samples) > 0:
        logging.error(
            f'Samples completely missing from AIP results: '
            f'{", ".join(missing_samples)}'
        )

        # for each of those missing samples, add all variants
        for miss_sample in missing_samples:
            discrepancies[miss_sample] = seqr_results[miss_sample]
            logging.error(
                f'Sample {miss_sample}: '
                f'{len(seqr_results[miss_sample])} missing variant(s)'
            )

    for sample in common_samples:

        # only finds discrepancies, not Matched results - revise
        sample_discrepancies = [
            variant
            for variant in seqr_results[sample]
            if variant not in aip_results[sample]
        ]

        # only populate the index if missing variants found
        if sample_discrepancies:
            discrepancies[sample] = sample_discrepancies

        # log the number of matches
        matched = len(
            [
                variant
                for variant in seqr_results[sample]
                if variant in aip_results[sample]
            ]
        )

        logging.info(f'Sample {sample} - {matched} matched variant(s)')
        logging.info(
            f'Sample {sample} - {len(sample_discrepancies)} missing variant(s)'
        )

    return discrepancies


def find_affected_samples(pedigree: Ped) -> list[str]:
    """
    finds all affected members of the provided pedigree

    :param pedigree:
    :return:
    """
    return [sam.sample_id for sam in pedigree.samples() if sam.affected]


def main(results: str, seqr: str, ped: str):
    """

    :param results:
    :param seqr:
    :param ped:
    :return:
    """

    # normalise data formats from AIP result file
    aip_json = read_json_from_path(results)
    result_dict = common_format_from_results(results_dict=aip_json)

    # Peddy can't read cloud paths, but we want the parsed Pedigree
    with open('i_am_a_temporary.ped', 'w', encoding='utf-8') as handle:
        handle.write(AnyPath(ped).read_text())
    pedigree_digest = Ped('i_am_a_temporary.ped')

    # Search for all affected sample IDs in the Pedigree
    affected = find_affected_samples(pedigree_digest)

    # parse the Seqr results table, specifically targeting variants in probands
    seqr_results = common_format_from_seqr(seqr=seqr, affected=affected)

    # compare the results of the two datasets
    _discrepancies = find_missing(seqr_results=seqr_results, aip_results=result_dict)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('--results')
    parser.add_argument('--seqr')
    parser.add_argument('--ped')
    args = parser.parse_args()
    main(results=args.results, seqr=args.seqr, ped=args.ped)
