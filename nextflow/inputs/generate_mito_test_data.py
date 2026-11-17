"""
script used to create an artificial Mito VCF for testing
"""

from pathlib import Path

from generate_test_data import HEADER, SAMPLES, Variant, bgzip_and_tabix, mat_inherited

OUTPUT_DIR = Path(__file__).parent

CONTIGS = {
    'chrM': 16569,
}

VARIANTS = [
    # cat 3, no canonical consequences, here to check normalisation
    Variant('chrM', 3243, 'A', 'G', mat_inherited),
    Variant('chrM', 8528, 'T', 'C', mat_inherited),
    Variant('chrM', 12278, 'T', 'C', mat_inherited),
]


def render(samples: list[str], variants: list[Variant]) -> str:
    """build the full VCF text for the given sample column order"""
    contigs = '\n'.join(f'##contig=<ID={contig},length={length},assembly=GRCh38>' for contig, length in CONTIGS.items())
    columns = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples])
    rows = '\n'.join(variant.as_row(samples) for variant in variants)
    return f'{HEADER.format(contigs=contigs)}{columns}\n{rows}\n'


def main():
    """write the joint Mito VCF, all bgzipped and tabixed"""

    bgzip_and_tabix(render(SAMPLES, VARIANTS), OUTPUT_DIR / 'joint')

    print(f'wrote {len(VARIANTS)} variants to joint_mito.vcf.bgz and {len(SAMPLES)} single-sample VCFs')


if __name__ == '__main__':
    main()
