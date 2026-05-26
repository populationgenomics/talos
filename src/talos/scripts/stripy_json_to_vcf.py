#!/usr/bin/env python3

"""
Borrowing _heavily_ from https://github.com/broadinstitute/gatk-sv/blob/7eb2af1feea9f81a390caff211f6832c6e7f7ac5/src/stripy/stripy_to_vcf.py
Has been extended to take multiple single-sample STRipy reports as input, and to characterise the GT field as:
- '0/0' 2 called alleles, no pathogenic
- '1/0' a1 pathogenic, a2 non-pathogenic
- '1/.' a1 pathogenic, a2 missing

Uses fixed contigs (GRCh38) for the header, instead of relying on a reference input

This script creates VCFs in the format expected by Talos. Specific fields used:
  - INFO/GENE              - used to detect the appropriate MOI for this STR
  - INFO/SYMBOL            - used in variant display
  - FORMAT/REPCN           - counts detected at each allele
  - FORMAT/STR_FILTER      - rule out quality-failed STR calls
  - FORMAT/DISEASE_DETAILS - dense formatted String, presenting each disease and corresponding thresholds
"""

import argparse
import json
import logging
import re
from pathlib import Path

import pysam

from cpg_utils import config, to_path


CONTIG_ORDER = [f'chr{x}' for x in list(range(1, 23))] + ['chrX', 'chrY', 'chrM']
HEADER_LINES = [
    ('chr1', 248956422),
    ('chr2', 242193529),
    ('chr3', 198295559),
    ('chr4', 190214555),
    ('chr5', 181538259),
    ('chr6', 170805979),
    ('chr7', 159345973),
    ('chr8', 145138636),
    ('chr9', 138394717),
    ('chr10', 133797422),
    ('chr11', 135086622),
    ('chr12', 133275309),
    ('chr13', 114364328),
    ('chr14', 107043718),
    ('chr15', 101991189),
    ('chr16', 90338345),
    ('chr17', 83257441),
    ('chr18', 80373285),
    ('chr19', 58617616),
    ('chr20', 64444167),
    ('chr21', 46709983),
    ('chr22', 50818468),
    ('chrX', 156040895),
    ('chrY', 57227415),
    ('chrM', 16569),
]

VCF_HEADER = {
    'INFO': [
        {'ID': 'END', 'Number': '1', 'Type': 'Integer', 'Description': 'Stop position of the interval'},
        {'ID': 'SVTYPE', 'Number': '1', 'Type': 'String', 'Description': 'Type of structural variant'},
        {'ID': 'RU', 'Number': '1', 'Type': 'String', 'Description': 'Repeat unit in the reference orientation'},
        {'ID': 'PERIOD', 'Number': '1', 'Type': 'Integer', 'Description': 'Length of the repeat unit'},
        {
            'ID': 'DISEASES',
            'Number': '.',
            'Type': 'String',
            'Description': 'Associated disease symbols for this locus (| separated)',
        },
        {'ID': 'LOCUS', 'Number': '1', 'Type': 'String', 'Description': 'Gene/locus identifier from STRipy'},
        {'ID': 'SYMBOL', 'Number': '1', 'Type': 'String', 'Description': 'Gene inferred from STRipy locus identifier'},
        {'ID': 'GENE', 'Number': '1', 'Type': 'String', 'Description': 'ENSG ID inferred from Gene Symbol'},
    ],
    'FORMAT': [
        {'ID': 'GT', 'Number': '1', 'Type': 'String', 'Description': 'Unphased genotype'},
        {'ID': 'REPCN', 'Number': '2', 'Type': 'Float', 'Description': 'Number of repeat units spanned by each allele'},
        {
            'ID': 'DISEASE_DETAILS',
            'Number': '1',
            'Type': 'String',
            'Description': '|-delimited details for each disease, in the form diseaseSymbol__normal__intermediate__pathogenic',  # noqa: E501
        },
        {
            'ID': 'REPCI1',
            'Number': '2',
            'Type': 'Integer',
            'Description': '95% CI min,max on repeat counts of first allele',
        },
        {
            'ID': 'REPCI2',
            'Number': '2',
            'Type': 'Integer',
            'Description': '95% CI min,max on repeat counts of second allele',
        },
        {
            'ID': 'OUTLIER',
            'Number': '2',
            'Type': 'Integer',
            'Description': 'Allelic population outlier flags (0/1) assigned by STRipy',
        },
        {
            'ID': 'ZSCORE',
            'Number': '2',
            'Type': 'Float',
            'Description': 'Allelic population Z-scores assigned by STRipy',
        },
        {'ID': 'DP', 'Number': '1', 'Type': 'Integer', 'Description': 'Total Depth'},
        {'ID': 'STR_FILTER', 'Number': '.', 'Type': 'String', 'Description': 'Filter status assigned by STRipy'},
    ],
}


def get_limited_locus_list() -> list[str] | None:
    """Try to get the finite list of loci from config, if available."""
    try:
        return config.config_retrieve(['stripy', 'loci_lists', 'default_with_exclusions'])
    except config.ConfigError:
        print('No CPG config detected, using all loci')
        return None


def parse_coords(coord: str) -> tuple[str, int, int] | None:
    m = re.match(r'(chr)?(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)', str(coord))
    if not m:
        return None
    chrom = m.group('chrom')
    if not chrom.startswith('chr'):
        chrom = 'chr' + chrom
    return chrom, int(m.group('start')), int(m.group('end'))


def _to_float(x: str | float) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def _to_outlier(x: str | int) -> int | None:
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _ci_tuple(allele: dict[str, dict[str, int]]) -> tuple[int, int]:
    return allele['CI'].get('Min', 0), allele['CI'].get('Max', 0)


def parse_disease_ranges(content: dict[str, str | dict[str, int]]) -> str:
    """
    Read the CorrespondingDisease content of the locus' JSON dictionary.
    This creates a hybrid String of Gene, MOI, normal and intermediate ranges, and pathogenic threshold.
    These occur per-disease, and there can be multiple per allele.

    example results: "HD__AD__min9max26__min27max35__36" or "DMD__XLR__min11max33__.__59"

    Args:
        content: the dictionary from STRipy's CorrespondingDisease data block

    Returns:
        A single aggregated String for all the details
    """
    # get the symbol for this disease
    disease = content['DiseaseSymbol']

    # the inheritance pattern
    inheritance = content['Inheritance']

    # parse the normal range, or '.'
    normal_range: dict[str, int] | str = content['NormalRange']
    normal = '.'
    if isinstance(normal_range, dict):
        normal = f'min{normal_range["Min"]}max{normal_range["Max"]}'

    # parse the intermediate range, or '.'
    inter_range: dict[str, int] | str = content['IntermediateRange']
    intermediate = '.'
    if isinstance(inter_range, dict):
        intermediate = f'min{inter_range["Min"]}max{inter_range["Max"]}'

    # add the pathogenic cutoff used for this disease, pull it all together
    return f'{disease}__{inheritance}__{normal}__{intermediate}__{content["PathogenicCutoff"]}'


def load_sample(json_path: str) -> tuple[str, dict]:
    """Load a STRipy JSON and return (sample_name, loci_dict keyed by locus_id)."""
    with to_path(json_path).open() as f:
        data = json.load(f)

    # hard-wired to pull the SampleID from the input filename - consider reverting
    input_file = data.get('JobDetails', {}).get('InputFile', '')
    if input_file:
        stem = Path(input_file).stem
        sample_name = stem.split('__')[0]
    else:
        sample_name = Path(json_path).stem.split('.')[0]

    loci = {}
    for entry in data.get('GenotypingResults', []):
        for locus_name, locus in entry.items():
            tl = locus['TargetedLocus']

            # create a |-delimited string for all ranges for all relevant disease for this locus
            disease_deets = '|'.join(parse_disease_ranges(corr) for corr in tl['CorrespondingDisease'].values())

            # add flexibility if the Alleles aren't populated at all
            if 'Alleles' not in locus:
                explain_string = f'Parsing {sample_name}, locus {tl["LocusID"]}, no Alleles present - skipping. '
                if 'Filter' in locus:
                    explain_string += f'Filter: {locus["Filter"]} - '
                logging.debug(explain_string)
                continue

            parsed = parse_coords(tl['Coordinates'])
            if not parsed:
                continue
            chrom, start, end = parsed
            alleles = locus['Alleles']
            a1 = alleles[0]
            a2 = alleles[1] if len(alleles) > 1 else None
            motif = tl['Motif']
            loci[locus_name] = {
                'chrom': chrom,
                'pos': start,
                'end': end,
                'id': locus_name,
                'motif': motif,
                'period': len(motif),
                'a1_rep': _to_float(a1['Repeats']),
                'a2_rep': _to_float(a2['Repeats']) if a2 is not None else None,
                'a1_ci': _ci_tuple(a1),
                'a2_ci': _ci_tuple(a2) if a2 is not None else (None, None),
                'a1_range': a1['Range'],
                'a2_range': a2['Range'] if a2 else None,
                'a1_out': _to_outlier(a1['IsPopulationOutlier']),
                'a2_out': _to_outlier(a2['IsPopulationOutlier']) if a2 is not None else None,
                'a1_z': _to_float(a1['PopulationZscore']),
                'a2_z': _to_float(a2['PopulationZscore']) if a2 is not None else None,
                'coverage': locus['Metadata']['Coverage'],
                'filter': locus['Filter'],
                'diseases': '|'.join(sorted({meta['DiseaseSymbol'] for meta in tl['CorrespondingDisease'].values()})),
                'disease_details': disease_deets,
            }

    return sample_name, loci


def get_header(sample_names: list[str]) -> pysam.VariantHeader:
    """Fill the header."""
    header = pysam.VariantHeader()
    header.add_meta('source', value='STRipy2VCF')
    header.add_meta('ALT', items=[('ID', 'STR'), ('Description', 'Short tandem repeat')])
    for info in VCF_HEADER['INFO']:
        header.add_meta('INFO', items=list(info.items()))
    for fmt in VCF_HEADER['FORMAT']:
        header.add_meta('FORMAT', items=list(fmt.items()))

    # generate the contig header lines
    for contig, length in HEADER_LINES:
        header.contigs.add(contig, length=length)

    for sample_name in sample_names:
        header.add_sample(sample_name)

    return header


def convert_range_to_gt(range_val: str | None) -> int | None:
    """Converts the per-allele range to an integer, for embedding in GT."""
    if range_val is None:
        return None
    if range_val.lower() == 'pathogenic':
        return 1
    return 0


def get_gene_lookup(map_path: str | None = None) -> dict[str, str]:
    """Pull out the dictionary reading logic."""
    lookup: dict[str, str] = {}
    if map_path is not None:
        with to_path(map_path).open() as handle:
            lookup = json.load(handle)
    return lookup


def write_multisample_vcf(  # noqa: PLR0915
    samples: list[tuple[str, dict[str, dict]]],
    out_path: str,
    gene_map: str | None = None,
) -> None:
    """
    Integrate the per-sample data into a union VCF

    Args:
        samples: list of (sample_name, loci_dict) tuples
        out_path: vcf path to write output to. PySam will write compressed if ending is gz/bgz
        gene_map: optional, file path to a JSON dict, mapping gene symbols to gene IDs
    """

    gene_lookup = get_gene_lookup(gene_map)

    header = get_header(sample_names=[sam_bit[0] for sam_bit in samples])

    only_loci = get_limited_locus_list()

    canonical: dict[str, dict] = {}

    for _, loci in samples:
        for locus_id, loc in loci.items():
            # if we got a list of specific loci, only use that finite list - ignore everything else
            if (only_loci is not None) and (locus_id not in only_loci):
                continue
            canonical.setdefault(locus_id, loc)

    # double sort contigs by both chrom and position
    sorted_loci = sorted(canonical.values(), key=lambda x: (CONTIG_ORDER.index(x['chrom']), x['pos']))

    with pysam.VariantFile(out_path, mode='w', header=header) as vf:
        for loc in sorted_loci:
            rec = header.new_record()
            rec.contig = loc['chrom']
            rec.start = loc['pos'] - 1
            rec.stop = loc['end']
            rec.id = str(loc['id'])
            rec.ref = 'N'
            rec.alts = ('<STR>',)
            rec.filter.add('PASS')
            rec.info['SVTYPE'] = 'STR'
            if loc['motif']:
                rec.info['RU'] = loc['motif']
            if loc['period']:
                rec.info['PERIOD'] = int(loc['period'])
            rec.info['DISEASES'] = loc['diseases']
            rec.info['LOCUS'] = loc['id']

            # extract the gene symbol from STRipy annotations
            symbol = loc['id'].split('_')[0]
            rec.info['SYMBOL'] = symbol
            rec.info['GENE'] = gene_lookup.get(symbol, symbol)

            for sample_name, loci in samples:
                s_loc = loci.get(loc['id'])
                s = rec.samples[sample_name]
                if s_loc is None:
                    s['GT'] = (None, None)
                    s['REPCN'] = (None, None)
                    s['REPCI1'] = (None, None)
                    s['REPCI2'] = (None, None)
                    s['OUTLIER'] = (None, None)
                    s['ZSCORE'] = (None, None)
                    s['DP'] = None
                    s['STR_FILTER'] = ['.']
                    s['DISEASE_DETAILS'] = '.'
                else:
                    s['GT'] = (convert_range_to_gt(s_loc['a1_range']), convert_range_to_gt(s_loc['a2_range']))
                    s['REPCN'] = (s_loc['a1_rep'], s_loc['a2_rep'])
                    s['REPCI1'] = s_loc['a1_ci']
                    s['REPCI2'] = s_loc['a2_ci']
                    s['OUTLIER'] = (s_loc['a1_out'], s_loc['a2_out'])
                    s['ZSCORE'] = (s_loc['a1_z'], s_loc['a2_z'])
                    s['DP'] = int(s_loc['coverage'])
                    s['STR_FILTER'] = [str(s_loc['filter'])] if s_loc['filter'] else ['PASS']
                    s['DISEASE_DETAILS'] = s_loc['disease_details']

            vf.write(rec)


def main() -> None:
    ap = argparse.ArgumentParser(description='Convert STRipy JSON output(s) into a multi-sample VCF.')
    ap.add_argument('--json', required=True, nargs='+', help='STRipy JSON report(s); one per sample')
    ap.add_argument('--output', required=True, help='Output VCF')
    ap.add_argument('--mapping', default=None, help='A JSON dict of Gene Symbol:Gene ID to update annotation')
    args = ap.parse_args()

    samples = [load_sample(p) for p in args.json]

    write_multisample_vcf(samples=samples, out_path=args.output, gene_map=args.mapping)


if __name__ == '__main__':
    main()
