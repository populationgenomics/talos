import pytest
from mendelbrot.pedigree_parser import PedigreeParser

import hail as hl

from talos.RunHailFiltering import annotate_category_de_novo


@pytest.mark.parametrize(
    'clinvar_talos,consequence,biotype,result',
    [
        (0, 'frameshift', 'protein_coding', ['male']),
        (1, 'synonymous', 'protein_coding', ['male']),
        (0, 'synonymous', 'protein_coding', ['missing']),
        (0, 'inframe_deletion', 'protein_coding', ['male']),
        (0, 'inframe_insertion', 'protein_coding', ['male']),
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
    dn_matrix = annotate_category_de_novo(dn_matrix, pedigree_data=PedigreeParser(pedigree_path))

    assert dn_matrix.info.categorysampledenovo.collect() == result


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
        _dn_matrix = annotate_category_de_novo(dn_matrix, pedigree_data=PedigreeParser(pedigree_path))
