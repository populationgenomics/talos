"""
unit testing collection for the hail MT methods

aim - reconstruct some fields into a struct
simulate the
"""

import os
import hail as hl

from hail.utils.java import FatalError
import pytest
import pandas as pd

from reanalysis.hail_filter_and_classify import (
    annotate_class_1,
    annotate_class_2,
    annotate_class_3,
    annotate_class_4,
)


PWD = os.path.dirname(__file__)
INPUT = os.path.join(PWD, 'input')

# contains a single variant at chr1:1, with minimal info
HAIL_VCF = os.path.join(INPUT, 'single_hail.vcf.bgz')

class1_keys = ['locus', 'clinvar_sig', 'clinvar_stars']
class2_keys = ['locus', 'clinvar_sig', 'cadd', 'revel', 'geneIds', 'consequence_terms']
class3_keys = ['locus', 'clinvar_sig', 'lof', 'consequence_terms']
class4_keys = [
    'locus',
    'cadd',
    'revel',
    'mutationtaster',
    'gerp',
    'eigen',
]
class_conf = {
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


@pytest.fixture(name='hail_matrix')
def fixture_hail_matrix():
    """
    loads the single variant as a matrix table
    :return:
    """
    try:
        hl.init(default_reference='GRCh38')
    except FatalError:
        print('failure - hail already initiated')
    return hl.import_vcf(HAIL_VCF, reference_genome='GRCh38')


@pytest.mark.parametrize(
    'values,classified',
    [
        ([hl_locus, 'pathogenic', 0], 0),
        ([hl_locus, 'pathogenic', 1], 1),
        ([hl_locus, 'conflicting_interpretations_of_pathogenicity', 1], 0),
        ([hl_locus, 'benign', 1], 0),
        ([hl_locus, 'pathogenic&something&else', 2], 1),
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
    row_dict = dict(zip(class1_keys, values))

    # create a single row dataframe using the input
    dataframe = pd.DataFrame(row_dict, index=[0])

    hail_table = hl.Table.from_pandas(dataframe, key='locus')
    anno_matrix = hail_matrix.annotate_rows(
        info=hail_matrix.info.annotate(
            clinvar_sig=hail_table[hail_matrix.locus].clinvar_sig,
            clinvar_stars=hail_table[hail_matrix.locus].clinvar_stars,
        )
    )

    anno_matrix = annotate_class_1(anno_matrix)
    assert anno_matrix.info.Class1.collect() == [classified]


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
    row_dict = dict(zip(class2_keys, values))

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

    anno_matrix = annotate_class_2(
        anno_matrix, config=class_conf, new_genes=hl.set(['GREEN'])
    )
    assert anno_matrix.info.Class2.collect() == [classified]


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
    row_dict = dict(zip(class3_keys, values))

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

    anno_matrix = annotate_class_3(anno_matrix, config=class_conf)
    assert anno_matrix.info.Class3.collect() == [classified]


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
def test_class_4_assignment(values, classified, hail_matrix):
    """
    :param values: value order in the class4_keys list
    :param classified: expected classification
    :param hail_matrix: test fixture
    """

    polyphen = values.pop()
    sift = values.pop()
    # cast the input as a dictionary
    row_dict = dict(zip(class4_keys, values))

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

    anno_matrix = annotate_class_4(anno_matrix, config=class_conf.get('in_silico'))
    assert anno_matrix.info.Class4.collect() == [classified]
