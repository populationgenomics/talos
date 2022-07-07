"""
CPG-home-brewed implementations of Hail built-ins
"""


import hail as hl


def custom_de_novo(
    matrix: hl.MatrixTable,
    pedigree: hl.Pedigree,
    ab_cut_off: float = 0.05,
    depth_minimum: int = 10,
) -> hl.Table:
    """
    implements a home-brewed version of the de novo identifier
    here we accept the variants as called by the variant caller
    :param matrix:
    :param pedigree:
    :param ab_cut_off:
    :param depth_minimum:
    """

    # parse family structure into the trio matrix format
    trio_matrix = hl.trio_matrix(matrix, pedigree=pedigree, complete_trios=True)

    dad = trio_matrix.father_entry
    mum = trio_matrix.mother_entry

    # conditional value expressions
    parents_ref = dad.GT.is_hom_ref() & mum.GT.is_hom_ref()
    proband_depth_fail = hl.sum(trio_matrix.proband_entry.AD) < depth_minimum
    parents_ab_fail = (dad.AD[1] / hl.sum(dad.AD) > ab_cut_off) | (
        mum.AD[1] / hl.sum(mum.AD) > ab_cut_off
    )

    # regional expressions
    autosomal = trio_matrix.locus.in_autosome()
    chr_x = trio_matrix.locus.in_x_nonpar()
    is_female = trio_matrix.is_female

    trio_matrix = trio_matrix.annotate_entries(
        de_novo=(
            hl.case()
            .when((proband_depth_fail) | (parents_ab_fail) | (~parents_ref), 0)
            .when(
                (autosomal | (chr_x & is_female))
                & (trio_matrix.proband_entry.GT.is_het()),
                1,
            )
            .when(
                (chr_x) & (~is_female) & (trio_matrix.proband_entry.GT.is_hom_var()),
                1,
            )
            .default(0)
        )
    )

    trio_matrix = trio_matrix.filter_entries(trio_matrix.de_novo == 1)
    return trio_matrix.entries()
