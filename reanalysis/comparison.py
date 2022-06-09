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
from cyvcf2 import VCFReader
import hail as hl
from peddy import Ped

from cpg_utils.hail_batch import init_batch

from reanalysis.utils import read_json_from_path, canonical_contigs_from_vcf


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
        always normalise contig - upper case, with no CHR prefix
        :param chrom:
        :param pos:
        :param ref:
        :param alt:
        :param confidence:
        """
        self.chr: str = self.normalise_chrom(chrom)
        self.pos: int = pos
        self.ref: str = ref
        self.alt: str = alt
        self.confidence: list[Confidence] = confidence

    @staticmethod
    def normalise_chrom(chrom: str) -> str:
        """
        normalise the string name
        :param chrom:
        :return:
        """
        up_chrom = chrom.upper()
        chrom = up_chrom[up_chrom.startswith('CHR') and 3 :]
        return chrom

    def get_cyvcf2_pos(self, contigs: set[str]) -> tuple[str, str]:
        """
        get variant coordinates in string format
        returns the coordinates for a cyvcf2 query, AND
        the contig name, normalised to the VCF
        :param contigs: all the contigs present in the VCF
        :return:
        """
        if self.chr in contigs:
            string_chr = self.chr
        else:
            string_chr = f'chr{self.chr}'

        if string_chr not in contigs:
            raise Exception(
                f'Contigs in this VCF not compatible with provided locus: {self}'
            )
        return string_chr, f'{string_chr}:{self.pos}-{self.pos}'

    def get_hail_pos(self) -> tuple[hl.Locus, list[str]]:
        """
        get relevant values for finding a variant row within a MT
        :return:
        """
        return hl.Locus(
            contig=f'chr{self.chr}', position=self.pos, reference_genome='GRCh38'
        ), [
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
ReasonDict = dict[str, list[tuple[CommonFormatResult, list[str]]]]


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
    This is a check for False Negatives
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


def check_in_vcf(vcf_path: str, variants: CommonDict) -> tuple[CommonDict, CommonDict]:
    """
    check whether each variant is present within the labelled VCF
    split list of discrepancies into two collections;
    - in VCF (investigate via MOI)
    - not in VCF (investigate within MT)
    :param vcf_path:
    :param variants:
    :return:
    """

    in_vcf: CommonDict = defaultdict(list)
    not_in_vcf: CommonDict = defaultdict(list)

    # open the VCF for reading (VCF needs to be localised in hail)
    vcf_handler = VCFReader(vcf_path)
    vcf_contigs = canonical_contigs_from_vcf(vcf_handler)

    # iterate over all samples, and corresponding lists
    for sample, var_list in variants.items():
        for var in var_list:

            # assume missing until confirmed otherwise
            found = False
            normalised_chrom, coordinates = var.get_cyvcf2_pos(vcf_contigs)
            for vcf_var in vcf_handler(coordinates):

                # check position and alleles
                if (
                    vcf_var.CHROM == normalised_chrom
                    and vcf_var.POS == var.pos
                    and vcf_var.REF == var.ref
                    and vcf_var.ALT[0] == var.alt
                ):
                    found = True
                    in_vcf[sample].append(var)
                    break
            if not found:
                not_in_vcf[sample].append(var)

    return in_vcf, not_in_vcf


def find_variant_in_mt(
    matrix: hl.MatrixTable, var: CommonFormatResult
) -> hl.MatrixTable:
    """
    :param matrix:
    :param var:
    :return:
    """
    var_locus, var_alleles = var.get_hail_pos()
    logging.info(f'Querying for {var_locus}: {var_alleles}')
    return matrix.filter_rows(
        (matrix.locus == var_locus) & (matrix.alleles == var_alleles)
    )


def check_mt(
    matrix: hl.MatrixTable, variants: CommonDict
) -> tuple[CommonDict, ReasonDict]:
    """
    for all variants provided, check which are present in the matrix table
    export two lists; not in the MT, and in the MT
    :param matrix:
    :param variants:
    :return:
    """

    not_in_mt: CommonDict = defaultdict(list)
    untiered: ReasonDict = defaultdict(list)

    for sample, variant in [
        (sam, var) for (sam, var_list) in variants.items() for var in var_list
    ]:
        logging.info(f'running check for {sample}, variant {variant}')

        # filter the MT to the required locus
        var_mt = find_variant_in_mt(matrix=matrix, var=variant)

        # check if there are remaining variant rows
        if var_mt.count_rows() == 0:
            not_in_mt[sample].append(variant)

        # if it is found, but was not in the VCF, it was untiered
        else:
            untiered[sample].append((variant, []))

    return not_in_mt, untiered


def main(results: str, seqr: str, ped: str, vcf: str, mt: str):
    """

    :param results:
    :param seqr:
    :param ped:
    :param vcf:
    :param mt:
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
    discrepancies = find_missing(seqr_results=seqr_results, aip_results=result_dict)

    in_vcf, not_in_vcf = check_in_vcf(vcf_path=vcf, variants=discrepancies)

    # some logging content here
    for sample in not_in_vcf:
        logging.info(f'Sample: {sample}')
        for variant in not_in_vcf[sample]:
            logging.info(f'\tVariant {variant} missing from VCF')

    # some logging content here
    for sample in in_vcf:
        logging.info(f'Sample: {sample}')
        for variant in in_vcf[sample]:
            logging.info(f'\tVariant {variant} requires MOI checking')

    # if there were any variants missing from the VCF, attempt to find them in the MT
    if len(not_in_vcf) > 0:

        # if we need to check the MT, start Hail Query
        init_batch(driver_cores=8, driver_memory='highmem')

        # read in the MT
        matrix = hl.read_matrix_table(mt)

        _untiered, _not_present = check_mt(matrix, not_in_vcf)


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
    parser.add_argument('--vcf')
    parser.add_argument('--mt')
    args = parser.parse_args()
    main(results=args.results, seqr=args.seqr, ped=args.ped, vcf=args.vcf, mt=args.mt)
