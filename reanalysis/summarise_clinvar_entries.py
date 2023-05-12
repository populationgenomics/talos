"""
read clinvar submissions; identify consensus and disagreement

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
import logging
import json
import re
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

import hail as hl
import pandas as pd

from cpg_utils import to_path, CloudPath
from cpg_utils.hail_batch import output_path, init_batch
from cpg_utils.config import ConfigError

from reanalysis.utils import get_cohort_config

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

MAJORITY_RATIO = 0.6
MINORITY_RATIO = 0.2
STRONG_REVIEWS = ['practice guideline', 'reviewed by expert panel']
ORDERED_ALLELES = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']

# published Nov 2015, available pre-print since March 2015
# assumed to be influential since 2016
ACMG_THRESHOLD = datetime(year=2016, month=1, day=1)
VERY_OLD = datetime(year=1970, month=1, day=1)
LARGEST_COMPLEX_INDELS = 40
BASES = re.compile(r'[ACGTN]+')


class Consequence(Enum):
    """
    csq enumeration
    """

    BENIGN = 'Benign'
    CONFLICTING = 'Conflicting'
    PATHOGENIC = 'Pathogenic'
    UNCERTAIN = 'VUS'
    UNKNOWN = 'Unknown'


# submitters not trusted for a subset of consequences - after Consequence is defined
try:
    cohort_conf = get_cohort_config()
    assert 'clinvar' in cohort_conf
    QUALIFIED_BLACKLIST = [
        (Consequence.BENIGN, cohort_conf['clinvar'].get('filter_benign', []))
    ]
except AssertionError:
    QUALIFIED_BLACKLIST = [
        (Consequence.BENIGN, ['illumina laboratory services; illumina'])
    ]


@dataclass
class Submission:
    """
    POPO to store details on each Submission
    """

    date: datetime
    submitter: str
    classification: Consequence
    review_status: str


def get_allele_locus_map(summary_file: str) -> dict:
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
        dictionary of each variant ID to the positional details
    """

    allele_dict = {}

    for line in lines_from_gzip(summary_file):
        if 'GRCh37' in line:
            continue

        # pull values from the line
        allele_id = int(line[0])
        chromosome = line[18] if 'chr' in line[18] else f'chr{line[18]}'
        var_id = int(line[30])
        pos = int(line[31])
        ref = line[32]
        alt = line[33]

        # skip chromosomal deletions and insertions, mito, or massive indels
        if (
            ref == 'na'
            or alt == 'na'
            or ref == alt
            or 'm' in chromosome.lower()
            or (len(ref) + len(alt)) > LARGEST_COMPLEX_INDELS
        ):
            continue

        # don't include any of the trash bases in ClinVar
        if BASES.match(ref) and BASES.match(alt):
            allele_dict[var_id] = {
                'allele': allele_id,
                'chrom': chromosome,
                'pos': pos,
                'ref': ref,
                'alt': alt,
            }

    return allele_dict


def lines_from_gzip(filename: str) -> str:
    """
    generator for gzip reading, copies file locally before reading

    Args:
        filename (str): the gzipped input file

    Returns:
        generator; yields each line
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


def consequence_decision(subs: list[Submission]) -> Consequence:
    """
    determine overall consequence assignment based on submissions

    Args:
        subs (): a list of submission objects for this allele

    Returns:
        a single Consequence object
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


def process_line(data: list[str]) -> tuple[int, Submission]:
    """
    takes a line, strips out useful content as a 'Submission'

    Args:
        data (): the array of line content

    Returns:
        the allele ID and corresponding Submission details
    """
    allele_id = int(data[0])
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
        submission_file (): file containing submission-per-line
        threshold_date (): ignore submissions after this date
        allele_ids (): only process alleleIDs we have pos data for

    Returns:
        dictionary of alleles and their corresponding submissions
    """

    submission_dict = defaultdict(list)

    # remove all entries from these providers
    # for now require a mandatory dataset, may amend later
    try:
        cohort_config = get_cohort_config()
        blacklist = cohort_config.get('clinvar_filter', [])
    except (AssertionError, KeyError):
        blacklist = []

    for line in lines_from_gzip(submission_file):
        a_id, line_sub = process_line(line)

        # skip rows where the variantID isn't in this mapping
        # this saves a little effort on haplotypes, CNVs, and SVs
        if (
            (a_id not in allele_ids)
            or (line_sub.submitter in blacklist)
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
    # default to returning everything
    return date_filt_subs or subs


def sort_decisions(all_subs: list[dict]) -> list[dict]:
    """
    applies dual-layer sorting to the list of all decisions

    Args:
        all_subs (): list of all submissions

    Returns:
        a list of submissions, sorted hierarchically on chr & pos
    """

    return sorted(
        all_subs, key=lambda x: (ORDERED_ALLELES.index(x['contig']), x['position'])
    )


def parse_into_table(json_path: str, out_path: str) -> hl.Table:
    """
    takes the file of one clinvar variant per line
    processes that line into a table based on the schema

    Args:
        json_path (): path to the JSON file (temp)
        out_path (): where to write the Hail table

    Returns:
        the Hail Table object created
    """

    # start a hail runtime
    init_batch()

    # define the schema for each written line
    schema = hl.dtype(
        'struct{'
        'alleles:array<str>,'
        'contig:str,'
        'position:int32,'
        'id:int32,'
        'clinical_significance:str,'
        'gold_stars:int32,'
        'allele_id:int32'
        '}'
    )

    # import the table, and transmute to top-level attributes
    ht = hl.import_table(json_path, no_header=True, types={'f0': schema})
    ht = ht.transmute(**ht.f0)

    # create a locus and key
    ht = ht.transmute(locus=hl.locus(ht.contig, ht.position))
    ht = ht.key_by(ht.locus, ht.alleles)

    # write out
    ht.write(out_path, overwrite=True)
    return ht


def snv_missense_filter(clinvar_table: hl.Table, vcf_path: str):
    """
    takes a clinvar table and a filters to SNV & Pathogenic
    Writes results to a VCF file

    Args:
        clinvar_table (hl.Table): the table of re-summarised clinvar loci
        vcf_path (str): Path to write resulting VCF to
    """

    # filter to Pathogenic SNVs
    # there is at least one ClinVar submission which is Pathogenic
    # without being a changed base?
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/1705890/
    clinvar_table = clinvar_table.filter(
        (hl.len(clinvar_table.alleles[0]) == 1)
        & (hl.len(clinvar_table.alleles[1]) == 1)
        & (clinvar_table.clinical_significance == 'Pathogenic')
    )

    # persist the clinvar annotations in VCF
    clinvar_table = clinvar_table.annotate(
        info=hl.struct(
            allele_id=clinvar_table.allele_id, gold_stars=clinvar_table.gold_stars
        )
    )
    hl.export_vcf(clinvar_table, vcf_path, tabix=True)
    logging.info(f'Wrote SNV VCF to {vcf_path}')


def main(
    subs: str, date: datetime, variants: str, out: str, path_snv: str | None = None
):
    """
    Redefines what it is to be a clinvar summary

    Args:
        subs (str): file path to all submissions (gzipped)
        variants (str): file path to variant summary (gzipped)
        out (str): path to write JSON out to
        date (str): date threshold to use for filtering submissions
        path_snv (str): if defined, path to write SNV VCF file
    """

    logging.info('Getting alleleID-VariantID-Loci from variant summary')
    allele_map = get_allele_locus_map(variants)

    logging.info('Getting all decisions, indexed on clinvar AlleleID')
    decision_dict = get_all_decisions(
        submission_file=subs, threshold_date=date, allele_ids=set(allele_map.keys())
    )

    # placeholder to fill wth per-allele decisions
    all_decisions = []

    # now filter each set of decisions per allele
    for allele_id, submissions in decision_dict.items():
        # filter against ACMG date, if appropriate
        submissions = acmg_filter_submissions(submissions)

        # obtain an aggregate rating
        rating = consequence_decision(submissions)

        # assess stars in remaining entries
        stars = check_stars(submissions)

        # for now, skip over variants which are not relevant to AIP
        if rating in [Consequence.UNCERTAIN, Consequence.UNKNOWN]:
            continue

        all_decisions.append(
            {
                'alleles': [allele_map[allele_id]['ref'], allele_map[allele_id]['alt']],
                'contig': allele_map[allele_id]['chrom'],
                'position': allele_map[allele_id]['pos'],
                'id': allele_id,
                'clinical_significance': rating.value,
                'gold_stars': stars,
                'allele_id': allele_map[allele_id]['allele'],
            }
        )

    # sort all collected decisions, trying to reduce overhead in HT later
    all_decisions = sort_decisions(all_decisions)

    # if there's no defined config, write a local file
    try:
        temp_output = output_path(
            f'{date.strftime("%Y-%m-%d")}_clinvar_table.json', category='tmp'
        )
    except (ConfigError, KeyError):
        temp_output = f'{date.strftime("%Y-%m-%d")}_clinvar_table.json'

    logging.info(f'temp JSON location: {temp_output}')

    # open this temp path and write the json contents, line by line
    with to_path(temp_output).open('w') as handle:
        for each_dict in all_decisions:
            handle.write(f'{json.dumps(each_dict)}\n')

    ht = parse_into_table(json_path=temp_output, out_path=out)

    if path_snv:
        logging.info('Writing out SNV VCF')
        snv_missense_filter(ht, path_snv)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    parser = ArgumentParser()
    parser.add_argument('-s', help='submission_summary.txt.gz from NCBI', required=True)
    parser.add_argument('-v', help='variant_summary.txt.gz from NCBI', required=True)
    parser.add_argument('-o', help='output table name', required=True)
    parser.add_argument(
        '-d',
        help=(
            'date, format YYYY-MM-DD. Individual submissions after this date are '
            'removed. Un-dated submissions will pass this threshold.'
        ),
        default=datetime.now(),
    )
    parser.add_argument('--path_snv', help='Output VCF, sites-only, Pathogenic SNVs')
    args = parser.parse_args()

    processed_date = (
        datetime.strptime(args.d, '%Y-%m-%d') if isinstance(args.d, str) else args.d
    )

    logging.info(f'Date threshold: {processed_date}')

    main(
        subs=args.s,
        variants=args.v,
        out=args.o,
        date=processed_date,
        path_snv=args.path_snv,
    )
