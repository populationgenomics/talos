"""
compares a Seqr file against results of a previous AIP run

The Seqr file is a TSV of all variants flagged as interesting
This is designed to recognise flags in the format 'AIP training: Confidence'

See relevant documentation for a description of the algorithm used
"""
import json
import os
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

from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch

from reanalysis.hail_filter_and_label import (
    extract_annotations,
    filter_by_ab_ratio,
    filter_matrix_by_ac,
    filter_on_quality_flags,
    filter_to_population_rare,
    filter_to_well_normalised,
    green_and_new_from_panelapp,
    CONFLICTING,
    LOFTEE_HC,
    PATHOGENIC,
)

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
        return (
            f'{self.chr}:{self.pos}_{self.ref}>{self.alt} '
            f'- {", ".join(map(str, sorted(self.confidence)))}'
        )

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
        return hash(repr(self))


CommonDict = dict[str, list[CommonFormatResult]]
ReasonDict = dict[str, list[tuple[CommonFormatResult, list[str]]]]


def common_format_aip(results_dict: dict[str, Any]) -> CommonDict:
    """
    Parses the JSON

    :param results_dict:
    :return:
    """
    sample_dict: CommonDict = defaultdict(list)

    # collect all per-sample results into a separate index
    for sample, variants in results_dict.items():

        for var in variants:
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


def common_format_seqr(seqr: str, affected: list[str]) -> CommonDict:
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


def find_seqr_flags(
    aip_results: CommonDict, seqr_results: CommonDict
) -> dict[str, dict[str, list[str] | int]]:
    """
    check for exact matches to the Seqr flags, and report numbers
    returns a per-confidence dictionary of the variant details and counts
    :param aip_results:
    :param seqr_results:
    :return:
    """

    flag_matches = {
        key: {
            'matched': {'details': [], 'count': 0},
            'unmatched': {'details': [], 'count': 0},
        }
        for key in ['EXPECTED', 'UNLIKELY', 'POSSIBLE']
    }
    total_seqr_variants = 0

    for sample, variants in seqr_results.items():
        if sample not in aip_results:
            for v in variants:
                for conf in v.confidence:
                    flag_matches[conf.name]['unmatched']['details'].append(
                        f'{sample}::{repr(v)}'
                    )
                    flag_matches[conf.name]['unmatched']['count'] += 1
                total_seqr_variants += 1
            continue

        aip_variants = aip_results[sample]
        for v in variants:
            total_seqr_variants += 1
            match = 'matched' if v in aip_variants else 'unmatched'
            for conf in v.confidence:
                flag_matches[conf.name][match]['details'].append(f'{sample}::{repr(v)}')
                flag_matches[conf.name][match]['count'] += 1

    # print a summary into logging
    logging.info(f'Total Seqr Variants: {total_seqr_variants}')
    for confidence, match_types in flag_matches.items():
        logging.info(f'{confidence}')
        for match_type, match_dict in match_types.items():
            logging.info(f'\t{match_type} - {match_dict["count"]}')

    return flag_matches


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


def check_gene_is_green(
    matrix: hl.MatrixTable, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    :param matrix:
    :param green_genes:
    :return:
    """
    return matrix.filter_rows(green_genes.contains(matrix.geneIds))


def run_ac_check(matrix: hl.MatrixTable) -> list[str]:
    """
    if there are enough samples in the joint call, run the AC test
    :param matrix:
    :return:
    """
    if (
        filter_matrix_by_ac(matrix, get_config()['filter']['ac_threshold']).count_rows()
        == 0
    ):
        return ['QC: AC too high in joint call']
    return []


def run_quality_flag_check(matrix: hl.MatrixTable) -> list[str]:
    """
    check the variant wasn't flagged by the caller/VQSR
    :param matrix:
    :return:
    """
    if filter_on_quality_flags(matrix).count_rows() == 0:
        return ['QC: Variant has assigned quality flags']
    return []


def check_variant_was_normalised(matrix: hl.MatrixTable) -> list[str]:
    """
    checks the variant is normalised
    :param matrix:
    :return:
    """

    if filter_to_well_normalised(matrix).count_rows() == 0:
        return ['QC: Variant not well normalised']
    return []


def filter_sample_by_ab(matrix: hl.MatrixTable, sample_id: str) -> list[str]:
    """

    :param matrix:
    :param sample_id:
    :return:
    """

    # evaluating the AB test has to be sample ID specific
    ab_filt_mt = filter_by_ab_ratio(matrix)
    if len(ab_filt_mt.filter_cols(ab_filt_mt.s == sample_id).entries().collect()) == 0:
        return ['QC: Variant fails AB ratio']

    return []


def check_population_rare(matrix: hl.MatrixTable) -> tuple[hl.MatrixTable, list[str]]:
    """
    filter out all rare variants, return fail reason if everything is removed
    :param matrix:
    :return:
    """
    matrix = filter_to_population_rare(matrix)

    if matrix.count_rows() == 0:
        return matrix, ['AF: No variants remain after Rare filter']
    return matrix, []


def check_cat_1(matrix: hl.MatrixTable) -> list[str]:
    """
    test against all conditions of category 1
    pretty primitive approach - apply a filter, count the rows...
    :param matrix:
    :return:
    """
    reasons = []
    if matrix.filter_rows(matrix.info.clinvar_stars > 0).count_rows() == 0:
        reasons.append('C1: ClinVar stars 0 or missing')
    if (
        matrix.filter_rows(
            matrix.info.clinvar_sig.lower().contains(PATHOGENIC)
        ).count_rows()
        == 0
    ):
        reasons.append('C1: Not ClinVar Pathogenic')
    if (
        matrix.filter_rows(
            matrix.info.clinvar_sig.lower().contains(CONFLICTING)
        ).count_rows()
        > 0
    ):
        reasons.append('C1: ClinVar rating Conflicting')
    return reasons


def filter_csq_to_set(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    overwrite the transcript consequences by only retaining those which
    have a consequence term overlapping with the supplied set
    :param matrix:
    :return:
    """
    critical_consequences = hl.set(get_config()['filter']['critical_csq'])

    # consequence_terms is an array - cast as a set
    csq_mt = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: hl.len(
                    critical_consequences.intersection(hl.set(x.consequence_terms))
                )
                > 0
            )
        )
    )

    # return rows which have remaining consequences
    return csq_mt.filter_rows(hl.len(csq_mt.vep.transcript_consequences) > 0)


def check_cadd_revel(matrix: hl.MatrixTable) -> int:
    """
    :param matrix:
    :return:
    """

    return matrix.filter_rows(
        (matrix.info.cadd > get_config()['filter']['in_silico']['cadd'])
        | (matrix.info.revel > get_config()['filter']['in_silico']['revel'])
    ).count_rows()


def check_cat_2(matrix: hl.MatrixTable, new_genes: hl.SetExpression) -> list[str]:
    """
    test against all conditions of category 2
    :param matrix:
    :param new_genes:
    :return:
    """
    reasons = []
    if matrix.filter_rows(new_genes.contains(matrix.geneIds)).count_rows() == 0:
        reasons.append('C2: Gene was not New in PanelApp')

    # filter to high CSQ, so count_rows shows whether there were rows left
    csq_filtered = filter_csq_to_set(matrix=matrix)
    if not csq_filtered.count_rows():
        reasons.append('C2: No HIGH CSQs')

        # stop here, nothing left to search
        return reasons

    if (
        matrix.filter_rows(
            matrix.info.clinvar_sig.lower().contains(PATHOGENIC)
        ).count_rows()
        == 0
    ):
        reasons.append('C2: Not ClinVar Pathogenic')

    if check_cadd_revel(matrix) == 0:
        reasons.append('C2: CADD & REVEL not significant')

    return reasons


def check_cat_3(matrix: hl.MatrixTable) -> list[str]:
    """
    test against all conditions of category 3
    :param matrix:
    :return:
    """
    reasons: list[str] = []

    # row-dependent annotation check
    if (
        matrix.filter_rows(
            matrix.info.clinvar_sig.lower().contains(PATHOGENIC)
        ).count_rows()
        == 0
    ):
        reasons.append('C3: No ClinVar Pathogenic')

    csq_filtered = filter_csq_to_set(matrix=matrix)
    # if there are no high consequences, no more rows to test; continue
    if not csq_filtered.count_rows():
        reasons.append('C3: No HIGH CSQs')
        return reasons

    # further filter to rows which aren't ruled out by LOFTEE
    if (
        csq_filtered.filter_rows(
            csq_filtered.vep.transcript_consequences.any(
                lambda x: (x.lof == LOFTEE_HC) | (hl.is_missing(x.lof))
            )
        ).count_rows()
        == 0
    ):
        reasons.append('C3: No LOFTEE High Confidence')

    return reasons


def check_cat_support(matrix: hl.MatrixTable) -> list[str]:
    """
    test against all conditions of support category
    :param matrix:
    :return:
    """
    reasons: list[str] = []

    if check_cadd_revel(matrix) == 0:
        reasons.append('Support: CADD & REVEL not significant')
    if (
        matrix.filter_rows(
            matrix.vep.transcript_consequences.any(
                lambda x: x.sift_score <= get_config()['filter']['in_silico']['sift']
            )
        ).count_rows()
        == 0
    ):
        reasons.append('Support: SIFT not significant')

    if (
        matrix.filter_rows(
            matrix.vep.transcript_consequences.any(
                lambda x: x.polyphen_score
                >= get_config()['filter']['in_silico']['polyphen']
            )
        ).count_rows()
        == 0
    ):
        reasons.append('Support: PolyPhen not significant')
    if (
        matrix.filter_rows(
            (matrix.info.mutationtaster.contains('D'))
            | (matrix.info.mutationtaster == 'missing')
        ).count_rows()
        == 0
    ):
        reasons.append('Support: No significant MutationTaster')

    return reasons


def prepare_mt(matrix: hl.MatrixTable) -> hl.MatrixTable:
    """
    prepare the MT (splitting by gene, moving attributes, etc.)
    replicates split_rows_by_gene_and_filter_to_green
    :param matrix:
    :return:
    """

    # split out the different gene annotations onto separate rows
    # but don't remove non-green genes yet
    matrix = matrix.explode_rows(matrix.geneIds)

    # csq filter to single gene per row
    matrix = matrix.annotate_rows(
        vep=matrix.vep.annotate(
            transcript_consequences=matrix.vep.transcript_consequences.filter(
                lambda x: (matrix.geneIds == x.gene_id)
                & ((x.biotype == 'protein_coding') | (x.mane_select.contains('NM')))
            )
        )
    )

    # move annotations, use original method
    return extract_annotations(matrix)


def check_mt(
    matrix: hl.MatrixTable,
    variants: CommonDict,
    green_genes: hl.SetExpression,
    new_genes: hl.SetExpression,
) -> tuple[CommonDict, ReasonDict]:
    """
    for all variants provided, check which are present in the matrix table
    export two lists; not in the MT, and in the MT

    There's currently a fair bit of repetition here, simplify

    :param matrix:
    :param variants:
    :param green_genes:
    :param new_genes:
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
        if not var_mt.count_rows():
            not_in_mt[sample].append(variant)
            continue

        # shift attributes in the MT, and explode on gene ID
        var_mt = prepare_mt(var_mt)

        # check for any failure reasons in the QC tests
        reasons: list[str] = []

        # run all methods relating to the quality of the variant call
        reasons.extend(run_ac_check(var_mt))
        reasons.extend(run_quality_flag_check(var_mt))
        reasons.extend(check_variant_was_normalised(var_mt))
        # # not actually a reason for failing
        # reasons.extend(filter_sample_by_ab(var_mt, sample))
        var_mt = check_gene_is_green(matrix=var_mt, green_genes=green_genes)
        if not var_mt.count_rows():
            reasons.append('Gene is not GREEN in PanelApp')

        # break early if we find a QC/Green Gene failure
        if reasons:
            untiered[sample].append((variant, reasons))
            continue

        # remove common variants
        var_mt, af_reason = check_population_rare(var_mt)

        # break early if we find a CSQ failure?
        if af_reason:
            reasons.extend(af_reason)
            untiered[sample].append((variant, reasons))
            continue

        # pass through the classification methods
        reasons.extend(check_cat_1(matrix=var_mt))
        reasons.extend(check_cat_2(matrix=var_mt, new_genes=new_genes))
        reasons.extend(check_cat_3(matrix=var_mt))
        reasons.extend(check_cat_support(matrix=var_mt))

        # log all reasons, even if the list is empty
        untiered[sample].append((variant, reasons))

        if not reasons:
            logging.error(f'No Fail reasons for Sample {sample}, ' f'Variant {variant}')

    return not_in_mt, untiered


def main(results_folder: str, pedigree: str, seqr: str, vcf: str, mt: str, output: str):
    """
    runs a full match-seeking analysis of this AIP run against the
    expected variants (based on seqr training flags)

    :param results_folder:
    :param pedigree:
    :param seqr:
    :param vcf:
    :param mt:
    :param output:
    :return:
    """

    # normalise data formats from AIP result file
    aip_results = common_format_aip(
        results_dict=read_json_from_path(
            os.path.join(results_folder, 'summary_results.json')
        )
    )

    # Search for all affected sample IDs in the Peddy Pedigree
    affected = find_affected_samples(Ped(pedigree))

    # parse the Seqr results table, specifically targeting variants in probands
    seqr_results = common_format_seqr(seqr=seqr, affected=affected)

    # strict comparison
    flag_summary = find_seqr_flags(aip_results=aip_results, seqr_results=seqr_results)
    with AnyPath(f'{output}_match_summary.json').open('w') as handle:
        json.dump(flag_summary, handle, default=str, indent=4)

    # compare the results of the two datasets
    discrepancies = find_missing(seqr_results=seqr_results, aip_results=aip_results)
    if not discrepancies:
        logging.info('All variants resolved!')
        sys.exit(0)

    # load and digest panel data
    panel_dict = read_json_from_path(os.path.join(results_folder, 'panelapp_data.json'))
    green_genes, new_genes = green_and_new_from_panelapp(panel_dict)

    # if we had discrepancies, bin into classified and misc.
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
    if len(not_in_vcf) == 0:
        sys.exit(0)

    # if we need to check the MT, start Hail Query
    init_batch(driver_cores=8, driver_memory='highmem')

    # read in the MT
    matrix = hl.read_matrix_table(mt)

    not_present, untiered = check_mt(
        matrix=matrix,
        variants=not_in_vcf,
        green_genes=green_genes,
        new_genes=new_genes,
    )
    if untiered:
        with AnyPath(f'{output}_untiered.json').open('w') as handle:
            json.dump(untiered, handle, default=str, indent=4)

    if not_present:
        with AnyPath(f'{output}_missing.json').open('w') as handle:
            json.dump(not_present, handle, default=str, indent=4)


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(module)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        stream=sys.stderr,
    )
    parser = ArgumentParser()
    parser.add_argument('--results_folder')
    parser.add_argument('--pedigree')
    parser.add_argument('--seqr')
    parser.add_argument('--vcf')
    parser.add_argument('--mt')
    parser.add_argument('--output')
    args = parser.parse_args()
    main(
        results_folder=args.results_folder,
        pedigree=args.pedigree,
        seqr=args.seqr,
        vcf=args.vcf,
        mt=args.mt,
        output=args.output,
    )
