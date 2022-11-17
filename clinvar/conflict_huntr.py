"""
reads over the clinvar submissions
identifies consensus and disagreements

Requires two files from
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

submission_summary.txt
 - all individual submissions to ClinVar
relevant fields:
1.  VariationID: the identifier assigned by ClinVar
2.  ClinicalSignificance:
7.  ReviewStatus: the level of review for this submission, namely
10. Submitter

variant_summary.txt
 - links clinvar AlleleID, Variant ID, position and alleles
"""


import gzip
import json
import logging
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

import hail as hl
import pandas as pd

from cpg_utils import to_path, CloudPath
from cpg_utils.hail_batch import output_path, init_batch


BENIGN_SIGS = {'Benign', 'Likely benign', 'Benign/Likely benign', 'protective'}
CONFLICTING = 'conflicting data from submitters'
PATH_SIGS = {
    'Pathogenic',
    'Likely pathogenic',
    'Pathogenic, low penetrance',
    'Likely pathogenic, low penetrance',
    'Pathogenic/Likely pathogenic',
}
UNCERTAIN_SIGS = {'Uncertain significance', 'Uncertain risk allele'}
USELESS_RATINGS = {'no assertion criteria provided'}
MEGA_BLACKLIST = [
    'victorian clinical genetics services,murdoch childrens research institute'
]
MAJORITY_RATIO = 0.6
MINORITY_RATIO = 0.2
STRONG_REVIEWS = ['practice guideline', 'reviewed by expert panel']
ORDERED_ALLELES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y', 'M']]

# published Nov 2015, available pre-print since March 2015
# assumed to be influential since 2016
ACMG_THRESHOLD = datetime(year=2016, month=1, day=1)
VERY_OLD = datetime(year=1970, month=1, day=1)
LARGEST_COMPLEX_INDELS = 40


class Consequence(Enum):
    """
    csq enumeration
    """

    BENIGN = 'Benign'
    CONFLICTING = 'Conflicting'
    PATHOGENIC = 'Pathogenic'
    UNCERTAIN = 'VUS'
    UNKNOWN = 'Unknown'


# submitters which aren't trusted for a subset of consequences
QUALIFIED_BLACKLIST = [(Consequence.BENIGN, {'Illumina Laboratory Services; Illumina'})]


@dataclass
class Submission:
    """
    POPO to store details on each Submission
    """

    date: datetime
    submitter: str
    classification: Consequence
    review_status: str


def get_allele_locus_map(summary_file: str) -> tuple[dict, set]:
    """
    Process variant_summary.txt
     - links the allele ID, Locus/Alleles, and variant ID
    relevant fields:
    0 #AlleleID
    20 Chromosome
    30 VariationID
    31 Start
    32 ReferenceAllele
    33 AlternateAllele

    Args:
        summary_file (str): path to the gzipped text file

    Returns:

    """

    allele_dict = {}
    all_contigs = set()

    for line in lines_from_gzip(summary_file):
        if 'GRCh37' in line:
            continue

        # pull values from the line
        allele_id = line[0]
        chromosome = line[18] if 'chr' in line[18] else f'chr{line[18]}'
        var_id = line[30]
        pos = int(line[31])
        ref = line[32]
        alt = line[33]

        # skip chromosomal deletions and insertions, mito, or massive indels
        # this might break hail?
        if (
            ref == 'na'
            or alt == 'na'
            or 'm' in chromosome.lower()
            or (len(ref) + len(alt)) > LARGEST_COMPLEX_INDELS
        ):
            continue

        all_contigs.add(chromosome)

        allele_dict[var_id] = {
            'allele': allele_id,
            'chrom': chromosome,
            'pos': pos,
            'ref': ref,
            'alt': alt,
        }

    return allele_dict, all_contigs


def lines_from_gzip(filename: str) -> str:
    """
    generator for gzip reading
    copies file to local prior to reading

    Args:
        filename (): the gzipped input file

    Returns:
        generator - yields each line
    """

    if isinstance(to_path(filename), CloudPath):
        tempfile = 'file.txt.gz'
        to_path(filename).copy(tempfile)
        filename = tempfile

    with gzip.open(filename, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            yield line.rstrip().split('\t')


def county_county(subs: list[Submission]) -> Consequence:
    """
    count different consequence assignments
    """
    # pylint: disable=R0912

    # start with a default consequence
    decision = Consequence.UNCERTAIN

    # establish counts
    benign = 0
    pathogenic = 0
    uncertain = 0
    total = 0

    for each_sub in subs:

        # for 3/4-star ratings, don't look any further
        if each_sub.review_status in STRONG_REVIEWS:
            return each_sub.classification

        total += 1
        if each_sub.classification == Consequence.PATHOGENIC:
            pathogenic += 1
        elif each_sub.classification == Consequence.BENIGN:
            benign += 1
        elif each_sub.classification == Consequence.UNCERTAIN:
            uncertain += 1

    # no entries - no decision
    if total == 0:
        return decision

    if pathogenic and benign:
        if pathogenic == benign:
            decision = Consequence.CONFLICTING

        elif (max(pathogenic, benign) >= (total * MAJORITY_RATIO)) and (
            min(pathogenic, benign) <= (total * MINORITY_RATIO)
        ):
            decision = (
                Consequence.BENIGN if benign > pathogenic else Consequence.PATHOGENIC
            )
        else:
            decision = Consequence.CONFLICTING

    elif uncertain >= total / 2:
        decision = Consequence.UNCERTAIN

    elif pathogenic:
        decision = Consequence.PATHOGENIC

    elif benign:
        decision = Consequence.BENIGN

    return decision


def check_stars(subs: list[Submission]) -> int:
    """
    processes the submissions, and assigns a 'gold star' rating
    Args:
        subs (): list of all submissions at this allele

    Returns:
        integer, summarising the rating
    """
    minimum = 0
    for sub in subs:
        if sub.review_status == 'practice guideline':
            return 4
        if sub.review_status == 'reviewed by expert panel':
            return 3
        if 'criteria provided' in sub.review_status:
            minimum = 1

    return minimum


def process_line(data: str) -> tuple[str, Submission]:
    """
    takes a line array and strips out useful content
    :param data: the un-split TSV content
    """
    allele_id = data[0]
    if data[1] in PATH_SIGS:
        classification = Consequence.PATHOGENIC
    elif data[1] in BENIGN_SIGS:
        classification = Consequence.BENIGN
    elif data[1] in UNCERTAIN_SIGS:
        classification = Consequence.UNCERTAIN
    else:
        classification = Consequence.UNKNOWN
    date = datetime.strptime(data[2], '%b %d, %Y') if data[2] != '-' else VERY_OLD
    sub = data[9].lower()
    rev_status = data[6].lower()

    return allele_id, Submission(date, sub, classification, rev_status)


def dict_list_to_ht(list_of_dicts: list) -> hl.Table:
    """
    takes the per-allele results and aggregates into a hl.Table

    Args:
        list_of_dicts ():

    Returns:
        Hail table of the same content, indexed on locus & alleles

    """

    # convert list of dictionaries to a DataFrame
    pdf = pd.DataFrame(list_of_dicts)

    # convert DataFrame to a Table, keyed on Locus & Alleles
    return hl.Table.from_pandas(pdf, key=['locus', 'alleles'])


def get_all_decisions(
    submission_file: str, threshold_date: datetime, allele_ids: set
) -> dict[str, list[Submission]]:
    """
    obtains all submissions per-allele which pass basic criteria
        - not a blacklisted submitter
        - not a csq-specific blacklisted submitter
        - not after the user-specified date

    Args:
        submission_file ():
        threshold_date ():
        allele_ids ():

    Returns:

    """
    submission_dict = defaultdict(list)

    for line in lines_from_gzip(submission_file):

        a_id, line_sub = process_line(line)

        # skip any rows where the variantID isn't in this mapping
        # this saves a little effort on haplotypes, CNVs, and SVs
        if (
            (a_id not in allele_ids)
            or (line_sub.submitter in MEGA_BLACKLIST)
            or (line_sub.date > threshold_date)
            or (line_sub.review_status in USELESS_RATINGS)
            or (line_sub.classification == Consequence.UNKNOWN)
        ):
            continue

        # screen out some submitters per-consequence
        for consequence, submitters in QUALIFIED_BLACKLIST:
            if (
                line_sub.classification == consequence
                and line_sub.submitter in submitters
            ):
                continue

        submission_dict[a_id].append(line_sub)

    return submission_dict


def acmg_filter_submissions(subs: list[Submission]) -> list[Submission]:
    """
    filter submissions by dates
    if any submissions for this variant occur after the ACMG introduction
        - only return those
    if not
        - return all submissions

    Just to remove the possibility of removing expert curations, we won't
    date filter any expert/manual entries this way

    Args:
        subs (): list of submissions

    Returns:
        either all submissions if none are after the ACMG cut-off
        or, if newer entries exist, only those after cut-off
    """

    # apply the date threshold to all submissions
    date_filt_subs = [
        sub
        for sub in subs
        if sub.date >= ACMG_THRESHOLD or sub.review_status in STRONG_REVIEWS
    ]

    # if this contains results, return only those
    if date_filt_subs:
        return date_filt_subs

    # default to returning everything
    return subs


def sort_decisions(all_subs: list[dict]) -> list[dict]:
    """
    applies dual-layer sorting to the list of all decisions
    Args:
        all_subs ():

    Returns:
        a list of submissions, sorted hierarchically on chr & pos
    """
    return sorted(
        all_subs, key=lambda x: (ORDERED_ALLELES.index(x['contig']), x['position'])
    )


def main(
    submissions_file: str,
    summary: str,
    out_path: str,
    threshold_date: datetime,
):
    """

    Args:
        submissions_file (): file path to all submissions (gzipped)
        summary (): file path to variant summary (gzipped)
        out_path (): path to write JSON out to
        threshold_date (): date threshold to use for filtering submissions
    """

    logging.info('Getting all alleleID-VariantID-Loci from variant summary')
    allele_map, _all_contigs = get_allele_locus_map(summary)

    logging.info('Getting all decisions, indexed on clinvar AlleleID')
    decision_dict = get_all_decisions(
        submission_file=submissions_file,
        threshold_date=threshold_date,
        allele_ids=set(allele_map.keys()),
    )

    # placeholder to fill wth per-allele decisions
    all_decisions = []

    # now filter each set of decisions per allele
    for allele_id, submissions in decision_dict.items():
        # filter against ACMG date
        submissions = acmg_filter_submissions(submissions)

        # obtain an aggregate rating
        rating = county_county(submissions)

        # assess stars in remaining entries
        stars = check_stars(submissions)

        # for now, skip over variants which are not relevant to AIP
        if rating in [Consequence.UNCERTAIN, Consequence.UNKNOWN]:
            continue

        all_decisions.append(
            {
                'locus': hl.Locus(
                    contig=allele_map[allele_id]['chrom'],
                    position=allele_map[allele_id]['pos'],
                ),
                'alleles': [allele_map[allele_id]['ref'], allele_map[allele_id]['alt']],
                'contig': allele_map[allele_id]['chrom'],
                'position': allele_map[allele_id]['pos'],
                'id': allele_id,
                'rating': rating.value,
                'stars': stars,
                'allele_id': allele_map[allele_id]['allele'],
            }
        )

    # sort all collected decisions, trying to reduce overhead in HT later
    all_decisions = sort_decisions(all_decisions)

    with to_path(output_path(out_path).replace('.mt', '.json')).open('w') as handle:
        json.dump(all_decisions, handle)

    ht = dict_list_to_ht(all_decisions)

    # write out the aggregated table
    ht.write(output_path(out_path), overwrite=True)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    parser = ArgumentParser()
    parser.add_argument('-s', help='submission_summary.txt.gz from NCBI', required=True)
    parser.add_argument('-v', help='variant_summary.txt.gz from NCBI', required=True)
    parser.add_argument(
        '-d',
        help='date, format DD-MM-YYYY - submissions after this date '
        'will be removed. Un-dated submissions will pass this threshold',
        default=datetime.now(),
    )
    parser.add_argument('-o', help='output filename', required=True)
    args = parser.parse_args()

    processed_date = (
        datetime.strptime(args.d, '%d-%m-%Y') if isinstance(args.d, str) else args.d
    )

    init_batch(worker_memory='highmem', worker_cores=8)

    main(
        submissions_file=args.s,
        summary=args.v,
        threshold_date=processed_date,
        out_path=args.o,
    )
