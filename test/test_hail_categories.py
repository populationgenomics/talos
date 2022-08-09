"""
unit testing collection for the hail MT methods
"""

import pytest

import hail as hl
import pandas as pd

from reanalysis.hail_filter_and_label import (
    annotate_category_1,
    annotate_category_2,
    annotate_category_3,
    annotate_category_5,
    annotate_category_support,
    green_and_new_from_panelapp,
    filter_to_population_rare,
    split_rows_by_gene_and_filter_to_green,
    filter_to_categorised,
)


category_1_keys = ['locus', 'clinvar_sig', 'clinvar_stars']
category_2_keys = [
    'locus',
    'clinvar_sig',
    'cadd',
    'revel',
    'geneIds',
    'consequence_terms',
]
category_3_keys = ['locus', 'clinvar_sig', 'lof', 'consequence_terms']
support_category_keys = [
    'locus',
    'cadd',
    'revel',
    'mutationtaster',
    'gerp',
    'eigen',
]

category_conf = {
    'critical_csq': ['frameshift_variant'],
    'in_silico': {
        'cadd': 28.0,
        'revel': 0.4,
        'polyphen': 0.85,
        'sift': 0.05,
        'gerp': 1.0,
        'eigen': 0.25,
    },
}

hl_locus = hl.Locus(contig='chr1', position=1, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'pathogenic', 0], 0),
        ([hl_locus, 'pathogenic', 1], 1),
        ([hl_locus, 'conflicting_interpretations_of_pathogenicity', 1], 0),
        ([hl_locus, 'benign', 1], 0),
        ([hl_locus, 'pathogenic&something&else', 2], 1),
        ([hl_locus, 'pathogenic&sbenignng&else', 2], 0),
    ],
)
def test_class_1_assignment(values, classified, hail_matrix):
    """
    use some fake annotations, apply to the single fake variant
    check that the classification process works as expected based
    on the provided annotations
    :param values:
    :param classified:
    :param hail_matrix:
    """
    # cast the input as a dictionary
    row_dict = dict(zip(category_1_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            clinvar_stars=hail_table[hail_matrix.locus].clinvar_stars,
        )
    )

    anno_matrix = annotate_category_1(anno_matrix)
    assert anno_matrix.info.categoryboolean1.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'pathogenic', 0.0, 0.0, 'GREEN', 'missense'], 1),
        ([hl_locus, 'benign', 0.0, 0.0, 'GREEN', 'missense'], 0),
        ([hl_locus, 'pathogenic', 99.0, 1.0, 'RED', 'frameshift_variant'], 0),
        ([hl_locus, 'pathogenic', 99.0, 1.0, 'GREEN', 'frameshift_variant'], 1),
        ([hl_locus, 'benign', 30.0, 0.0, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 28.01, 0.0, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 0, 0.41, 'GREEN', 'synonymous'], 1),
        ([hl_locus, 'meh', 0, 0.0, 'GREEN', 'synonymous'], 0),
    ],
)
def test_class_2_assignment(values, classified, hail_matrix):
    """
    the fields in the input are, respectively:
    :param values: locus clinvar_sig cadd revel geneIds consequence_terms
    :param classified:
    :param hail_matrix:
    """

    csq = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(category_2_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=hail_table[hail_matrix.locus].geneIds,
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            cadd=hail_table[hail_matrix.locus].cadd,
            revel=hail_table[hail_matrix.locus].revel,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [hl.Struct(consequence_terms=hl.set([csq]))]
            ),
        ),
    )

    anno_matrix = annotate_category_2(
        anno_matrix, config=category_conf, new_genes=hl.set(['GREEN'])
    )
    assert anno_matrix.info.categoryboolean2.collect() == [classified]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'benign', 'hc', 'frameshift_variant'], 0),
        ([hl_locus, 'benign', 'HC', 'frameshift_variant'], 1),
        ([hl_locus, 'pathogenic', 'lc', 'frameshift_variant'], 1),
        ([hl_locus, 'pathogenic', hl.missing(hl.tstr), 'frameshift_variant'], 1),
    ],
)
def test_class_3_assignment(values, classified, hail_matrix):
    """
    the fields in the input are, respectively:
    :param values: locus clinvar_sig loftee consequence_terms
    :param classified:
    :param hail_matrix:
    """

    csq = values.pop()
    lof = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(category_3_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        consequence_terms=hl.set([csq]),
                        lof=lof,
                    )
                ]
            ),
        ),
    )

    anno_matrix = annotate_category_3(anno_matrix, config=category_conf)
    assert anno_matrix.info.categoryboolean3.collect() == [classified]


@pytest.mark.parametrize(
    'spliceai_score,flag',
    [(0.1, 0), (0.11, 0), (0.3, 0), (0.49, 0), (0.5, 1), (0.69, 1), (0.9, 1)],
)
def test_category_5_assignment(spliceai_score: float, flag: int, hail_matrix):
    """
    :param hail_matrix:
    :return:
    """

    matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(splice_ai_delta=spliceai_score)
    )
    conf = {'spliceai_full': 0.5}
    matrix = annotate_category_5(matrix, conf)
    assert matrix.info.categoryboolean5.collect() == [flag]


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 0.0, 0.0, 'n', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.9, 'n', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 29.5, 0.9, 'n', 0.0, 0.0, 1.0, 0.0], 1),
        ([hl_locus, 0.0, 0.0, 'D', 0.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.0, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 1.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 0.0, 0.0], 0),
        ([hl_locus, 0.0, 0.0, 'D', 10.0, 0.5, 0.0, 0.9], 1),
    ],
)
def test_support_assignment(values, classified, hail_matrix):
    """
    :param values: value order in the class4_keys list
    :param classified: expected classification
    :param hail_matrix: test fixture
    """

    polyphen = values.pop()
    sift = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(support_category_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            cadd=hail_table[hail_matrix.locus].cadd,
            eigen_phred=hail_table[hail_matrix.locus].eigen,
            gerp_rs=hail_table[hail_matrix.locus].gerp,
            mutationtaster=hail_table[hail_matrix.locus].mutationtaster,
            revel=hail_table[hail_matrix.locus].revel,
        ),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        polyphen_score=polyphen,
                        sift_score=sift,
                    )
                ]
            ),
        ),
    )

    anno_matrix = annotate_category_support(anno_matrix, config=category_conf)
    assert anno_matrix.info.categorysupport.collect() == [classified]


def test_green_and_new_from_panelapp(panel_changes):
    """
    check that the set expressions from panelapp data are correct
    this is collection of ENSG names from panelapp
    2 set expressions, one for all genes, one for new genes only
    :param panel_changes:
    """
    green_expression, new_expression = green_and_new_from_panelapp(panel_changes)

    # check types
    assert isinstance(green_expression, hl.SetExpression)
    assert isinstance(new_expression, hl.SetExpression)

    # check content by collecting
    assert sorted(list(green_expression.collect()[0])) == [
        'ENSG00ABCD',
        'ENSG00EFGH',
        'ENSG00IJKL',
    ]
    assert list(new_expression.collect()[0]) == ['ENSG00EFGH']


@pytest.mark.parametrize(
    'exomes,genomes,length',
    [
        (0, 0, 1),
        (1.0, 0, 0),
        (0.04, 0.04, 1),
    ],
)
def test_filter_rows_for_rare(exomes, genomes, length, hail_matrix):
    """
    :param hail_matrix:
    :return:
    """
    conf = {'af_semi_rare': 0.05}
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            gnomad_ex_af=exomes,
            gnomad_af=genomes,
        )
    )
    matrix = filter_to_population_rare(anno_matrix, conf)
    assert matrix.count_rows() == length


@pytest.mark.parametrize(
    'gene_ids,length',
    [
        ({'not_green'}, 0),
        ({'green'}, 1),
        ({'gene'}, 1),
        ({'gene', 'not_green'}, 1),
        ({'green', 'gene'}, 2),
        ({hl.missing(t=hl.tstr)}, 0),
    ],
)
def test_filter_to_green_genes_and_split(gene_ids, length, hail_matrix):
    """
    :param hail_matrix:x
    :return:
    """
    green_genes = hl.literal({'green', 'gene'})
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=hl.literal(gene_ids),
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [hl.Struct(gene_id='gene', biotype='protein_coding', mane_select='')]
            ),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == length


def test_filter_to_green_genes_and_split__consequence(hail_matrix):
    """
    :param hail_matrix:
    :return:
    """

    green_genes = hl.literal({'green'})
    anno_matrix = hail_matrix.annotate_rows(
        geneIds=green_genes,
        vep=hl.Struct(
            transcript_consequences=hl.array(
                [
                    hl.Struct(
                        gene_id='green', biotype='protein_coding', mane_select=''
                    ),
                    hl.Struct(gene_id='green', biotype='batman', mane_select='NM_Bane'),
                    hl.Struct(gene_id='green', biotype='non_coding', mane_select=''),
                    hl.Struct(
                        gene_id='NOT_GREEN', biotype='protein_coding', mane_select=''
                    ),
                ]
            ),
        ),
    )
    matrix = split_rows_by_gene_and_filter_to_green(anno_matrix, green_genes)
    assert matrix.count_rows() == 1
    matrix = matrix.filter_rows(hl.len(matrix.vep.transcript_consequences) == 2)
    assert matrix.count_rows() == 1


@pytest.mark.parametrize(
    'one,two,three,four,five,support,length',
    [
        (0, 0, 0, 'missing', 0, 0, 0),
        (0, 1, 0, 'missing', 0, 0, 1),
        (0, 0, 1, 'missing', 0, 0, 1),
        (0, 0, 0, 'missing', 0, 1, 1),
        (0, 0, 0, 'not_blank', 0, 0, 1),
        (0, 0, 0, 'missing', 1, 0, 1),
        (0, 1, 1, 'missing', 0, 1, 1),
        (1, 0, 0, 'missing', 0, 1, 1),
    ],
)
def test_filter_to_classified(
    one, two, three, four, five, support, length, hail_matrix
):
    """
    :param hail_matrix:
    """
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            categoryboolean1=one,
            categoryboolean2=two,
            categoryboolean3=three,
            categorysample4=four,
            categoryboolean5=five,
            categorysupport=support,
        )
    )
    matrix = filter_to_categorised(anno_matrix)
    assert matrix.count_rows() == length
