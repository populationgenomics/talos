import pytest

import hail as hl

from talos.pedigree_parser import PedigreeParser
from talos.RunHailFiltering import annotate_category_4, pad_homref_ad


@pytest.mark.parametrize(
    'clinvar_talos,consequence,biotype,result',
    [
        (0, 'frameshift', 'protein_coding', ['male']),
        (1, 'synonymous', 'protein_coding', ['male']),
        (0, 'synonymous', 'protein_coding', ['missing']),
    ],
)
def test_dn_working(clinvar_talos, consequence, biotype, result, make_a_de_novo_mt, pedigree_path):
    """check that the de novo annotation works"""
    dn_matrix = make_a_de_novo_mt.annotate_rows(
        info=make_a_de_novo_mt.info.annotate(clinvar_talos=clinvar_talos),
        transcript_consequences=hl.array(
            [
                hl.Struct(consequence=consequence, biotype=biotype),
            ],
        ),
    )
    dn_matrix = annotate_category_4(dn_matrix, pedigree_data=PedigreeParser(pedigree_path))

    assert dn_matrix.info.categorysample4.collect() == result


def test_dn_bch_one(make_a_bch_de_novo_mt, pedigree_path):
    """check that the de novo annotation works"""
    dn_matrix = make_a_bch_de_novo_mt.annotate_rows(
        info=make_a_bch_de_novo_mt.info.annotate(clinvar_talos=1),
        transcript_consequences=hl.array(
            [
                hl.Struct(consequence='frameshift', biotype='protein_coding'),
            ],
        ),
    )

    with pytest.raises(hl.utils.java.HailUserError):
        _dn_matrix = annotate_category_4(dn_matrix, pedigree_data=PedigreeParser(pedigree_path))


def test_dn_bch_working(make_a_bch_de_novo_mt, pedigree_path):
    """check that the de novo annotation works"""
    dn_matrix = make_a_bch_de_novo_mt.annotate_rows(
        info=make_a_bch_de_novo_mt.info.annotate(clinvar_talos=1),
        transcript_consequences=hl.array(
            [
                hl.Struct(consequence='frameshift', biotype='protein_coding'),
            ],
        ),
    )

    pad_homref_ad(dn_matrix)
    dn_matrix = annotate_category_4(dn_matrix, pedigree_data=PedigreeParser(pedigree_path), strict_ad=True)
    assert dn_matrix.info.categorysample4.collect() == ['missing']
