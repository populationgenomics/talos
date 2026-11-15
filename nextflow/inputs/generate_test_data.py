"""
script used to create an artificial VCF for testing

Hand-writes the joint-called trio VCF and its three single-sample derivatives directly, in the
same style as generate_sv_test_data.py. Hail and talos.data_model have both been deleted, but the
annotation payload they used to attach (VEP consequences, gnomAD, ClinVar, SpliceAI) never actually
reached the VCF anyway - the old hl.export_vcf() call only had a slot for `info` and per-sample
FORMAT fields, so every other field on the Hail row was silently dropped on export. The only things
that have to round-trip here are locus, alleles, INFO AC/AF/AN, and the per-sample
GT/AD/DP/GQ/PL/PS values run_small_filtering.py and friends read - annotation gets attached
downstream by the real annotation.nf workflow (bcftools csq, echtvar), not by this script.

Genotype/PL defaults mirror the deleted data_model.Entry: AD defaults to [15, 15] (DP=30), GQ
defaults to 60, and PL is chosen by how many ALT alleles the GT carries (homref/het/homalt). The
magic ad=[-1] sentinel from Entry (used below for the "tricky" de novo calls) meant "AD/DP unknown"
in Hail terms - it is rendered here as a plain missing '.' field.

Run from the repository root:
    python nextflow/inputs/generate_test_data.py
"""

import subprocess
from dataclasses import dataclass
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent
INDIVIDUAL_DIR = OUTPUT_DIR / 'individual_vcfs'

# matches the sample columns already checked in to nextflow/inputs/joint.vcf.bgz
SAMPLES = ['proband', 'mother', 'father']

# only the contigs actually used below - lengths taken from the real GRCh38 reference, matching
# ##contig lines emitted by e.g. AnnotateSvWithGatk's GTF (see generate_sv_test_data.py)
CONTIGS = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr6': 170805979,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr16': 90338345,
    'chrX': 156040895,
}

HEADER = """\
##fileformat=VCFv4.2
{contigs}
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in the joint call">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele frequency in the joint call">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number in the joint call">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the observed alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=.,Type=Integer,Description="Phred-scaled genotype likelihoods">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
"""

# every variant here uses the same notional joint-call frequency - none of the original
# BaseFields() calls in the deleted script overrode it
INFO = 'AC=1;AF=0.001;AN=1000'


@dataclass
class Call:
    """one sample's genotype call, mirroring the deleted data_model.Entry defaults"""

    gt: str
    ad: list[int] | None = None
    gq: int = 60

    def format(self) -> str:
        """render as the GT:AD:DP:GQ:PL:PS FORMAT string"""
        gt = './.' if self.gt == '.' else self.gt
        # Entry's magic ad=[-1] sentinel meant "no read-depth data" - render as missing
        missing_depth = self.ad == [-1]
        ad = '.' if missing_depth else ','.join(str(x) for x in (self.ad or [15, 15]))
        dp = '.' if missing_depth else str(sum(self.ad or [15, 15]))
        # PL is picked by how many ALT alleles the GT carries: homref, het, homalt
        pl = [[0, self.gq, 1000], [self.gq, 0, 1000], [1000, self.gq, 0]][self.gt.count('1')]
        return f'{gt}:{ad}:{dp}:{self.gq}:{",".join(str(x) for x in pl)}:.'


HET_CALL = Call('0/1', ad=[15, 15])
HOM_CALL = Call('1/1', ad=[0, 30])
HOMREF_CALL = Call('0/0', ad=[30, 0])
MISSING_CALL = Call('.', ad=[-1])
# AD is Number=R - one value per allele regardless of ploidy, so haploid calls still carry a REF and an ALT depth
HEMI_VAR = Call('1', ad=[0, 30])
HEMI_REF = Call('0', ad=[30, 0])

# --- shared genotype patterns from the trio, matching the deleted script's fixtures ---
comp_het = {'proband': HOM_CALL, 'mother': HET_CALL, 'father': HET_CALL}
de_novo = {'proband': HET_CALL, 'mother': HOMREF_CALL, 'father': HOMREF_CALL}
tricky_de_novo = {'proband': HET_CALL, 'mother': MISSING_CALL, 'father': MISSING_CALL}
mat_inherited = {'proband': HET_CALL, 'mother': HET_CALL, 'father': HOMREF_CALL}
pat_inherited = {'proband': HET_CALL, 'mother': HOMREF_CALL, 'father': HET_CALL}
# this is for a parsing check that we can accurately process hemizygous alleles
mat_inherited_hemi = {'proband': HEMI_VAR, 'mother': HET_CALL, 'father': HEMI_REF}


@dataclass
class Variant:
    """one hand-authored small variant: its locus, alleles, and each trio member's call"""

    contig: str
    pos: int
    ref: str
    alt: str
    genotypes: dict[str, Call]

    def as_row(self, samples: list[str]) -> str:
        """render as a VCF data line, restricted to the given sample columns"""
        columns = [
            self.contig,
            self.pos,
            '.',
            self.ref,
            self.alt,
            '60',
            'PASS',
            INFO,
            'GT:AD:DP:GQ:PL:PS',
            *(self.genotypes[sample].format() for sample in samples),
        ]
        return '\t'.join(str(column) for column in columns)


VARIANTS = [
    # cat 3, no canonical consequences, here to check normalisation
    Variant('chr1', 21706892, 'GA', 'GAA', comp_het),
    # cat 1, recessive, HFE - not retained by the filtering logic when HFE is a phenotype-match only gene
    Variant('chr6', 26090951, 'C', 'G', comp_het),
    # PKHD1 (AR) | v4a pm5 | v4b cat1 | both should be reported
    Variant('chr6', 52043699, 'T', 'A', de_novo),
    Variant('chr6', 52043102, 'C', 'G', mat_inherited),
    # DARS (AR) | v3a cat3 | v3b cat6 | v3a should be reported, v3b supporting only
    Variant('chr2', 135920591, 'G', 'C', pat_inherited),
    Variant('chr2', 135912503, 'G', 'A', mat_inherited),
    # cat 6, dominant, DAAM2 - will be filtered out
    Variant('chr6', 39887558, 'C', 'T', de_novo),
    # Cat 1, 3, dominant, WT1
    Variant('chr11', 32392032, 'G', 'A', de_novo),
    Variant('chr11', 62841775, 'T', 'C', de_novo),
    # cats 1 5 3, POC1B, AR
    Variant('chr12', 89470359, 'A', 'C', comp_het),
    # IL2RG (hemi/bi in females) | cat 1 | should be reported
    Variant('chrX', 71109321, 'G', 'A', mat_inherited_hemi),
    # rnu4-2 (AD/AR) | cat 1, de novo | should be reported
    Variant('chr12', 120291834, 'A', 'G', de_novo),
    # HSPA8, not in any panels
    Variant('chr11', 123057736, 'A', 'AATC', tricky_de_novo),
    # HSPA8, not in any panels
    Variant('chr16', 89279566, 'CCTTCGGGG', 'C', tricky_de_novo),
]


def render(samples: list[str], variants: list[Variant]) -> str:
    """build the full VCF text for the given sample column order"""
    contigs = '\n'.join(f'##contig=<ID={contig},length={length},assembly=GRCh38>' for contig, length in CONTIGS.items())
    columns = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *samples])
    rows = '\n'.join(variant.as_row(samples) for variant in variants)
    return f'{HEADER.format(contigs=contigs)}{columns}\n{rows}\n'


def bgzip_and_tabix(text: str, stem: Path) -> None:
    """write `text` to `{stem}.vcf`, then bgzip and tabix it in place"""
    plain = stem.with_suffix('.vcf')
    plain.write_text(text)

    bgzipped = plain.with_suffix('.vcf.bgz')
    with bgzipped.open('wb') as handle:
        subprocess.run(['bgzip', '-c', str(plain)], stdout=handle, check=True)  # noqa: S603, S607
    subprocess.run(['tabix', '-p', 'vcf', str(bgzipped)], check=True)  # noqa: S603, S607
    plain.unlink()


def main():
    """write the joint VCF and one single-sample VCF per trio member, all bgzipped and tabixed"""
    INDIVIDUAL_DIR.mkdir(exist_ok=True)

    # VARIANTS is ordered by intent, not coordinate - tabix requires each contig's rows contiguous
    # and sorted by position, so group by contig in the order CONTIGS declares them
    contig_order = list(CONTIGS)
    ordered = sorted(VARIANTS, key=lambda variant: (contig_order.index(variant.contig), variant.pos))

    bgzip_and_tabix(render(SAMPLES, ordered), OUTPUT_DIR / 'joint')
    for sample in SAMPLES:
        bgzip_and_tabix(render([sample], ordered), INDIVIDUAL_DIR / sample)

    print(f'wrote {len(VARIANTS)} variants to joint.vcf.bgz and {len(SAMPLES)} single-sample VCFs')


if __name__ == '__main__':
    main()
