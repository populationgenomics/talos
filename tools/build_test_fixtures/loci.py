"""Single source of truth for the test fixture loci.

Every fixture (mini GFF, mini FASTA, mini AlphaMissense, echtvar zips, panelapp JSON,
clinvarbitration HTs, MANE subset) is derived from this list. Update this module if
nextflow/inputs/generate_test_data.py or generate_mito_test_data.py changes.

Each entry mirrors the variants emitted by those scripts. The `note` field documents
the test scenario the variant exercises. The `expected_gene` field records the gene
symbol the prod GFF assigns at that position - it is *not* used to drive subsetting
(coordinates do that), but documents intent so that reviewers can sanity-check that
the mini GFF actually covers each test scenario.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Locus:
    chrom: str
    pos: int
    ref: str
    alt: str
    expected_gene: str | None
    note: str


# Nuclear loci from nextflow/inputs/generate_test_data.py
NUCLEAR_LOCI: list[Locus] = [
    Locus('chr1', 21706892, 'GA', 'GAA', None, 'normalisation check, no canonical csq'),
    Locus('chr6', 26090951, 'C', 'G', 'HFE', 'cat 1 recessive, dropped when HFE is phenotype-match-only'),
    Locus('chr6', 52043699, 'T', 'A', 'PKHD1', 'PKHD1 (AR), part-a, pm5 (de novo)'),
    Locus('chr6', 52043102, 'C', 'G', 'PKHD1', 'PKHD1 (AR), part-b, cat1 (mat inherited)'),
    Locus('chr2', 135920591, 'G', 'C', 'DARS1', 'DARS (AR), part-a cat3 (pat inherited)'),
    Locus('chr2', 135912503, 'G', 'A', 'DARS1', 'DARS (AR), part-b cat6 supporting (mat inherited)'),
    Locus('chr6', 39887558, 'C', 'T', 'DAAM2', 'cat 6 dominant, will be filtered out'),
    Locus('chr11', 32392032, 'G', 'A', 'WT1', 'cat 1+3 dominant, de novo'),
    Locus('chr11', 62841775, 'T', 'C', 'RNU2-2P', 'rnu22 dominant de novo'),
    Locus('chr12', 89470359, 'A', 'C', 'POC1B', 'POC1B AR, cat 1+3+5 comp het'),
    Locus('chrX', 71109321, 'G', 'A', 'IL2RG', 'IL2RG hemi/bi female, cat 1, mat inherited hemi'),
    Locus('chr12', 120291834, 'A', 'G', 'RNU4-2', 'RNU4-2 AD/AR cat 1 de novo'),
    Locus('chr11', 123057736, 'A', 'AATC', None, 'tricky de novo, gene off-panel (HSPA8 region)'),
    Locus('chr16', 89279566, 'CCTTCGGGG', 'C', None, 'tricky de novo deletion, gene off-panel'),
]

# Mitochondrial loci from nextflow/inputs/generate_mito_test_data.py
MITO_LOCI: list[Locus] = [
    Locus('chrM', 3243, 'A', 'G', 'MT-TL1', 'MELAS m.3243A>G, mat inherited'),
    Locus('chrM', 8528, 'T', 'C', 'MT-ATP8', 'NARP-region m.8528T>C, mat inherited'),
    Locus('chrM', 12278, 'T', 'C', 'MT-TL2', 'tRNA-Leu(CUN) m.12278T>C, mat inherited'),
]

ALL_LOCI: list[Locus] = NUCLEAR_LOCI + MITO_LOCI

# Pedigree HPO terms from nextflow/inputs/pedigree.ped
PEDIGREE_HPO_TERMS: list[str] = [
    'HP:0002779',
    'HP:0004322',
    'HP:0009145',
]

# Flanking window (bp) of real reference bases to retain around each locus
# in the mini-FASTA, so bcftools csq can resolve codons / splice regions.
FASTA_FLANK_BP: int = 2000


def expected_gene_symbols() -> set[str]:
    """Gene symbols documented as the intended hit for each test variant.

    Used as a cross-check, not as the subset key. The actual subset is driven by
    coordinate-overlap against the prod GFF.
    """
    return {locus.expected_gene for locus in ALL_LOCI if locus.expected_gene}
