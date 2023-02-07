"""
tests relevant to the home-brewed methods file
"""

import hail as hl

from reanalysis.homebrewed import custom_de_novo


def test_de_novo_calling(trio_ped, de_novo_matrix):
    """
    test that we correctly identify the de novo variant(s)
    """
    hail_pedigree = hl.Pedigree.read(str(trio_ped))
    table = custom_de_novo(de_novo_matrix, hail_pedigree)
    table = table.key_by(table.locus, table.alleles)
    table = table.select(table.id).collect_by_key()
    table = table.annotate(values=hl.delimit(hl.map(lambda x: x.id, table.values), ','))
    # chr1 11111111 het in 1, wt in 2
    # chr2 22222222 wt in 1, het in 2
    # chr3 33333333 hom in 1, wt in 2
    # chr20 63406931 het in both
    # chrX 55555555 hom in male, het in female

    for locus, expected in [
        (hl.Locus('chr1', 11111111, reference_genome='GRCh38'), ['PROBAND1']),
        (hl.Locus('chr2', 22222222, reference_genome='GRCh38'), ['PROBAND2']),
        (
            hl.Locus('chr20', 63406931, reference_genome='GRCh38'),
            ['PROBAND1,PROBAND2'],
        ),
        (
            hl.Locus('chrX', 55555555, reference_genome='GRCh38'),
            ['PROBAND1,PROBAND2'],
        ),
    ]:
        tiny_table = table.filter(table.locus == locus)
        assert tiny_table.values.collect() == expected

    assert (
        table.filter(
            table.locus == hl.Locus('chr3', 33333333, reference_genome='GRCh38')
        ).count()
        == 0
    )
