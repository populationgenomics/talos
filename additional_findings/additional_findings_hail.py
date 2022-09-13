"""
a hail classifier specific to the additional findings process
"""

from argparse import ArgumentParser

import hail as hl

from reanalysis.utils import read_json_from_path
from reanalysis.hail_methods import (
    checkpoint_and_repartition,
    filter_to_population_rare,
    filter_on_quality_flags,
    split_rows_by_gene_and_filter_to_green,
)
from reanalysis.hail_filter_and_label import (
    PATHOGENIC,
    CONFLICTING,
    ONE_INT,
    MISSING_INT,
)


# this will be a dictionary mapping the exact consequences with the VEP CSQs
# currently the only category that exists is 'truncating'
VARIANT_CONSEQUENCES = {
    'truncating': hl.set(
        [
            'frameshift_variant',
            'splice_acceptor_variant',
            'splice_donor_variant',
            'stop_gained',
        ]
    )
}


# not sure how best to structure this whole thing...
# either we take the MT, annotate it, and return it
# which seems increasingly wasteful as the search goals get hyper-specific
# or process the MT, get the loci, and annotate the flags back against those loci....


def add_top_level_clinvar_flag(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    bop bop
    Parameters
    ----------
    mt :

    Returns
    -------
    same MT with two fields = clinvar yes/no, and a keep==no
    """
    return mt.annotate_rows(
        clinvar=hl.if_else(
            (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
            & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING)),
            ONE_INT,
            MISSING_INT,
        ),
        keep=MISSING_INT,
    )


def read_and_filter_mt(
    mt_path: str, config: dict, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    import the MT from file, and run some standard filters
    recycling as much as possible, inc splitting CSQ by gene

    Parameters
    ----------
    mt_path : path to the MT to use for input
    config : all settings for this run
    green_genes : ENSG values to filter to

    Returns
    -------
    the filtered MT

    TODO add a field audit call...
    TODO work out what would go in that audit...
    """

    mt = hl.read_matrix_table(mt_path)
    mt = filter_on_quality_flags(mt)
    mt = filter_to_population_rare(mt, config['filter_af'])
    return split_rows_by_gene_and_filter_to_green(mt, green_genes=green_genes)


def find_standard_interesting(mt: hl.MatrixTable, acmg_data: dict[str, dict]):
    """
    Parameters
    ----------
    mt :
    acmg_data : ACMG data

    Returns
    -------
    same MT, with clinvar P/LP highlighted
    """

    plain_ids = hl.set(
        [
            key
            for key, value in acmg_data.items()
            if not any([value.get('specific_type'), value.get('specific_variant')])
        ]
    )

    mt = mt.annotate_rows(
        info=mt.info.annotate(
            categorybooleangeneral=hl.if_else(
                (plain_ids.contains(mt.geneIds)) & (mt.clinvar == ONE_INT),
                ONE_INT,
                MISSING_INT,
            )
        )
    )

    # either update to 1, or keep original value
    return mt.annotate_rows(
        keep=hl.if_else(mt.info['categorybooleangeneral'], ONE_INT, mt.keep)
    )


def find_specific_types(
    mt: hl.MatrixTable, acmg_data: dict[str, dict]
) -> hl.MatrixTable:
    """
    Parameters
    ----------
    mt :
    acmg_data : ACMG data

    Returns
    -------
    same MT, with categories flagged where the variant types are appropriate
    """

    specific_types = {
        key: value for key, value in acmg_data.items() if value.get('specific_type')
    }
    types = {gene['specific_type'] for gene in specific_types.values()}
    for var_type in types:
        assert (
            var_type in VARIANT_CONSEQUENCES
        ), f'this script doesn\'t tolerate {var_type}'  # so fix it
        consequences = VARIANT_CONSEQUENCES[var_type]
        genes_for_this_consequence = hl.set(
            [key for key, value in specific_types if var_type in value['specific_type']]
        )
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **{
                    f'categoryboolean{var_type}': hl.if_else(
                        (genes_for_this_consequence.contains(mt.geneIds))
                        & (
                            hl.len(
                                mt.vep.transcript_consequences.filter(
                                    lambda x: hl.len(
                                        consequences.intersection(
                                            hl.set(x.consequence_terms)
                                        )
                                    )
                                    > 0
                                )
                            )
                            > 0
                        )
                        & (mt.clinvar == ONE_INT),
                        ONE_INT,
                        MISSING_INT,
                    )
                }
            )
        )
        # either update to 1, or keep original value
        mt = mt.annotate_rows(
            keep=hl.if_else(mt.info[f'categoryboolean{var_type}'], ONE_INT, mt.keep)
        )

    return mt


def find_exact_variants(
    mt: hl.MatrixTable, acmg_data: dict[str, dict]
) -> hl.MatrixTable:
    """
    Parameters
    ----------
    mt : the whole matrix table
    acmg_data : ACMG data

    Returns
    -------
    same MT, with categories flagged where the exact variants were found
    """

    specific_variants = {
        key: value['specific_variant']
        for key, value in acmg_data.items()
        if value.get('specific_variant')
    }

    # this approach is pretty miserable... won't scale
    # some sexy _case_ work could solve
    for gene, variants in specific_variants.items():
        variant_set = hl.set(variants)
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **{
                    f'categoryboolean{gene}specific': hl.if_else(
                        (mt.geneIds == gene)
                        & (
                            hl.len(
                                mt.vep.transcript_consequences.filter(
                                    lambda x: variant_set.contains(x.hgvsp)
                                )
                            )
                            > 0
                        ),
                        ONE_INT,
                        MISSING_INT,
                    )
                }
            )
        )

        # either update to 1, or keep original value
        mt = mt.annotate_rows(
            keep=hl.if_else(mt.info[f'categoryboolean{gene}specific'], ONE_INT, mt.keep)
        )

    return mt


def filter_mt_to_keep(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    reduces the MT to relevant rows only
    Parameters
    ----------
    mt : all variants in the MT

    Returns
    -------
    the same MT, filtered to all rows where 'keep' is set to 1
    keep is set to 1 when any individual category is approved
    """
    return mt.filter_rows(mt.keep == ONE_INT)


def main(mt_path: str, acmg_path: str, output: str):
    """
    main method
    """

    # get the relevant acmg data in full
    acmg_data = read_json_from_path(acmg_path)

    # read the MT & do filtering
    mt = read_and_filter_mt(
        mt_path=mt_path,
        config={'filter_af': 0.05},
        green_genes=hl.set(acmg_data.keys()),
    )

    # assign a clinvar boolean flag once to simplify logic
    mt = add_top_level_clinvar_flag(mt)
    mt = checkpoint_and_repartition(
        mt,
        checkpoint_root='',
        extra_logging='after applying filters and splitting genes',
    )

    # NOW FOR THE FUN PART!
    mt = find_standard_interesting(mt=mt, acmg_data=acmg_data)
    mt = find_specific_types(mt=mt, acmg_data=acmg_data)
    mt = find_exact_variants(mt=mt, acmg_data=acmg_data)
    mt = filter_mt_to_keep(mt)

    print(mt.count())
    print(output)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--mt', help='the MT to use as input')
    parser.add_argument('--acmg', help='the parsed ACMG data')
    parser.add_argument('-o', help='VCF output path/prefix')
    args = parser.parse_args()
    print(args)
