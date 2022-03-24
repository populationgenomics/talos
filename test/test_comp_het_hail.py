"""
specific test class for the compound het calculation

input is
2-sample VCF containing 7 dummy variants, all het,
3 in one gene (both samples het for all)
2 in a second gene (only one sample het)
2 in a third gene, but both flagged as c4 only (only one sample het)

| Gene  |  Locus   | C4_only | Sample1 | Sample2 |
| Gene1 | chr1:1   |    0    |    X    |    X    |
| Gene1 | chr1:21  |    0    |    X    |    X    |
| Gene1 | chr1:41  |    X    |    X    |    0    |
| Gene2 | chr1:91  |    0    |    X    |    0    |
| Gene2 | chr1:921 |    0    |    X    |    0    |
| Gene3 | chr1:111 |    X    |    X    |    0    |
| Gene3 | chr1:112 |    X    |    X    |    0    |

all variants are GC>G
"""

import hail as hl
import pandas as pd

from reanalysis.hail_filter_and_categorise import extract_comp_het_details


def test_comp_hets(hail_comp_het):
    """
    applies multiple annotations to a multi-variant VCF
    pulls out all compound-het combinations from the VCF
    CH-check is conditional on no variant pair being 2 Class4s

    :param hail_comp_het:
    :return:
    """

    # define annotations as a list of lists
    annotation_cols = ['locus', 'gene_id', 'class_4_only']
    annotations = [
        [hl.Locus(contig='chr1', position=1, reference_genome='GRCh38'), 'Gene1', 0],
        [hl.Locus(contig='chr1', position=21, reference_genome='GRCh38'), 'Gene1', 0],
        [hl.Locus(contig='chr1', position=41, reference_genome='GRCh38'), 'Gene1', 1],
        [hl.Locus(contig='chr1', position=91, reference_genome='GRCh38'), 'Gene2', 0],
        [hl.Locus(contig='chr1', position=921, reference_genome='GRCh38'), 'Gene2', 1],
        [hl.Locus(contig='chr1', position=1111, reference_genome='GRCh38'), 'Gene3', 1],
        [hl.Locus(contig='chr1', position=1112, reference_genome='GRCh38'), 'Gene3', 1],
    ]

    # transform to a DataFrame, then a hl.Table
    anno_table = hl.Table.from_pandas(
        pd.DataFrame(
            [dict(zip(annotation_cols, each_anno)) for each_anno in annotations]
        ),
        key='locus',
    )

    # apply those annotations to the VCF content
    hail_comp_het = hail_comp_het.annotate_rows(
        class_4_only=anno_table[hail_comp_het.locus].class_4_only,
        info=hail_comp_het.info.annotate(
            gene_id=anno_table[hail_comp_het.locus].gene_id
        ),
    )

    # do the compound het check
    ch_results = extract_comp_het_details(hail_comp_het)

    # may need to go through and sort/set the lists for comparison
    expected = {
        'SAMPLE': {
            'Gene1': [
                ['1-1-GC-G', '1-21-GC-G'],
                ['1-1-GC-G', '1-41-GC-G'],
                ['1-21-GC-G', '1-1-GC-G'],
                ['1-21-GC-G', '1-41-GC-G'],
                ['1-41-GC-G', '1-1-GC-G'],
                ['1-41-GC-G', '1-21-GC-G'],
            ],
            'Gene2': [['1-91-GC-G', '1-921-GC-G'], ['1-921-GC-G', '1-91-GC-G']],
        },
        'SAMPLE2': {'Gene1': [['1-1-GC-G', '1-21-GC-G'], ['1-21-GC-G', '1-1-GC-G']]},
    }
    assert ch_results == expected
