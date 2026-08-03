"""
script used to create an artificial joint-called SV VCF for testing the SV_ANNOTATION workflow

Unlike generate_test_data.py and generate_mito_test_data.py this does not go via Hail. The variants here are
symbolic SV records carrying INFO fields Hail's VCF export would not round-trip, so the VCF is written out
directly and compressed with bgzip.

The variants are deliberately adversarial rather than nominal - see the Testing guidance section of
docs/SvAnnotationPlan.md. Coordinates are taken from the real reference data:

- gene spans come from MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz, the same GTF AnnotateSvWithGatk uses
- gnomAD records come from SVAFotate_SV_popAFs.GRCh38.v4.1.bed.gz, the same BED AnnotateSvWithSvafotate uses.
  BED START is 0-based and VCF POS is 1-based, so a query matching a BED record exactly is POS = START + 1
  and END = END

Nothing in the annotation chain writes AC, AF, AN, N_HET, N_HOMALT, MALE_AF or FEMALE_AF - they have to be
present in the joint-called input already, or run_hail_filtering_sv raises at rearrange_annotations(). They
are written here for that reason, with MALE_AF/FEMALE_AF array-typed because filter_matrix_by_ac subscripts
them.

Run from the repository root:
    python nextflow/inputs/generate_sv_test_data.py
"""

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

OUTPUT = Path(__file__).parent / 'joint_sv.vcf'

# matches nextflow/inputs/pedigree.ped - proband and father are male, mother is female
SAMPLES = ['proband', 'father', 'mother']

HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Breakend">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome of the second breakend">
##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of the second breakend">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Callers supporting this variant">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in the joint call">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency in the joint call">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number in the joint call">
##INFO=<ID=N_HET,Number=1,Type=Integer,Description="Number of heterozygous samples">
##INFO=<ID=N_HOMALT,Number=1,Type=Integer,Description="Number of homozygous-alt samples">
##INFO=<ID=MALE_AF,Number=A,Type=Float,Description="Allele frequency in male samples">
##INFO=<ID=FEMALE_AF,Number=A,Type=Float,Description="Allele frequency in female samples">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""


# the notional joint call these frequency fields come from. The VCF itself holds one trio, because that is
# what the test pedigree describes, but AC/AF/AN and the sex-specific AFs describe the whole joint call - a
# real cohort SV VCF is a subset of a much larger callset. Deriving them from three samples instead would put
# every variant above RunHailFiltering.callset_af_sv_recessive (0.03) and filter_matrix_by_ac would empty the
# callset
JOINT_CALL_AN = 2000


@dataclass
class SV:
    """one symbolic SV record: its coordinates, the trio's genotypes, and its joint-call frequency"""

    identifier: str
    pos: int
    end: int
    svtype: str
    svlen: int
    # proband, father, mother
    genotypes: list[str]
    # allele frequency in the notional joint call, NOT in the three samples below
    callset_af: float
    # only set for breakends
    alt: str | None = None
    extra_info: dict[str, str] = field(default_factory=dict)

    def info(self) -> str:
        """build the INFO column, expanding the joint-call frequency into the fields Talos reads"""
        allele_count = round(self.callset_af * JOINT_CALL_AN)

        fields = {
            'SVTYPE': self.svtype,
            'SVLEN': self.svlen,
            'END': self.end,
            'ALGORITHMS': 'depth',
            'AC': allele_count,
            'AF': self.callset_af,
            'AN': JOINT_CALL_AN,
            # every carrier in the joint call taken as heterozygous, which is enough for these fields to be
            # present and self-consistent - nothing in Talos cross-checks them against the genotypes
            'N_HET': allele_count,
            'N_HOMALT': 0,
            'MALE_AF': self.callset_af,
            'FEMALE_AF': self.callset_af,
            **self.extra_info,
        }
        return ';'.join(f'{key}={value}' for key, value in fields.items())

    def as_row(self) -> str:
        """render as a VCF data line"""
        alt = self.alt or f'<{self.svtype}>'
        columns = ['chr1', self.pos, self.identifier, 'N', alt, '.', 'PASS', self.info(), 'GT', *self.genotypes]
        return '\t'.join(str(column) for column in columns)


VARIANTS = [
    # THE positive control. A deletion spanning the whole of PADI6 (chr1:17372196-17401699), which is a
    # PanelApp gene. Every gnomAD DEL in that window is either far smaller or far larger, so nothing clears
    # -f 0.5 - and the highest-ranked candidate, gnomAD-SV_v3_DEL_chr1_7b67fb38 at OFP 0.13, is exactly the
    # spurious Best_gnomAD_ID the rename has to withhold. Rare in gnomAD and rare in the callset, so this is
    # the one variant expected to survive to the report as CategoryBooleanSV1
    # --- this was correctly detected as Cat. SV1, but the gene is Biallelic.
    SV('lof_del_padi6', 17372196, 17401699, 'DEL', 29503, ['1/1', '0/0', '0/0'], callset_af=0.001),
    # matches gnomAD-SV_v3_DEL_chr1_143f6519 (BED 1:509967-690000) exactly, AF 0.2238. It carries
    # PREDICTED_LOF=OR4F16, so it survives the LoF filter and has to be dropped on gnomAD frequency by
    # filter_matrix_by_af - it is rare in the callset, so nothing else can catch it
    SV('common_del', 509968, 690000, 'DEL', 180032, ['0/1', '0/1', '0/0'], callset_af=0.001),
    # the mirror image: matches gnomAD-SV_v3_INV_chr1_b34a7097 (BED 1:1171037-1205428) at AF 8e-06, so it is
    # vanishingly rare in gnomAD, but common in the joint call. It carries PREDICTED_LOF=TNFRSF18, another
    # PanelApp gene, so only filter_matrix_by_ac can drop it
    SV('inv_common_in_callset', 1171038, 1205428, 'INV', 34390, ['0/1', '0/0', '0/0'], callset_af=0.2),
    # a 200bp deletion buried inside that same 180kb gnomAD DEL. Reciprocal overlap is ~0.001, so it must NOT
    # inherit the AF - and its Best_gnomAD_ID must be withheld by the rename's gnomAD_Count gate. This is the
    # case a mis-set -f would break
    SV('small_del_in_large', 600001, 600200, 'DEL', 200, ['0/1', '0/0', '0/0'], callset_af=0.001),
    # overlaps four gnomAD duplications with very different frequencies, the extremes being
    # gnomAD-SV_v3_DUP_chr1_eb9a525a (AF 0.704) and gnomAD-SV_v3_DUP_chr1_34304cb4 (AF 8e-06).
    # Max_AF must take the highest, which is what makes it the conservative choice
    SV('dup_two_matches', 21489967, 21494856, 'DUP', 4889, ['0/1', '0/0', '0/1'], callset_af=0.001),
    # wholly inside GJA9 (chr1:38874069-38881587) and covering its first exon, so it annotates as
    # PREDICTED_INTRAGENIC_EXON_DUP and never acquires a PREDICTED_LOF
    SV('intragenic_dup_gja9', 38874500, 38876500, 'DUP', 2001, ['0/1', '0/1', '0/0'], callset_af=0.001),
    # an insertion where END == POS, which is how callers emit them and which --ins exists to handle. It
    # matches gnomAD-SV_v3_INS_chr1_3f94b1dc, proving the expansion works
    SV('ins_end_eq_pos', 66340, 66340, 'INS', 161, ['0/1', '0/0', '0/1'], callset_af=0.001),
    # a breakend. SVAnnotate never assigns PREDICTED_LOF to a BND - this one lands intronic in PADI6 and
    # exonic in GJA9 - so it is dropped by run_hail_filtering_sv before any frequency filter is consulted.
    # See the Risks section of docs/SvAnnotationPlan.md
    SV(
        'bnd_1',
        17390000,
        17390001,
        'BND',
        -1,
        ['0/1', '0/0', '0/1'],
        callset_af=0.001,
        alt='N]chr1:38875000]',
        extra_info={'CHR2': 'chr1', 'END2': '38875000'},
    ),
]


def main():
    """write the VCF, then bgzip and tabix it in place"""
    # VARIANTS is ordered by intent, not coordinate - tabix requires position order
    rows = '\n'.join(variant.as_row() for variant in sorted(VARIANTS, key=lambda variant: variant.pos))
    columns = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *SAMPLES])
    OUTPUT.write_text(f'{HEADER}{columns}\n{rows}\n')

    bgzipped = OUTPUT.with_suffix('.vcf.bgz')
    with bgzipped.open('wb') as handle:
        subprocess.run(['bgzip', '-c', str(OUTPUT)], stdout=handle, check=True)  # noqa: S603, S607
    subprocess.run(['tabix', '-p', 'vcf', str(bgzipped)], check=True)  # noqa: S603, S607
    OUTPUT.unlink()

    print(f'wrote {len(VARIANTS)} variants to {bgzipped}')


if __name__ == '__main__':
    main()
