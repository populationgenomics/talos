"""
A streaming filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category
This expects data annotated by GATK's SVAnnotate tool

This is a cyvcf2 process - no Hail, no Spark - despite the module name, which is retained
so that existing invocation paths keep working until the Hail-era names are retired.

CategoryBooleanSV1:
- rare
- predicted LoF in a listed gene
"""

from argparse import ArgumentParser
from typing import Any

from cyvcf2 import VCF, Writer
from loguru import logger

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path
from talos.vcf_streaming import (
    first_value,
    header_has_field,
    parse_pedigree,
    subset_reader_to_pedigree,
    variant_is_pass,
)

GNOMAD_POP = config_retrieve(['RunHailFilteringSv', 'gnomad_population'], 'gnomad_v4.1')

# INFO fields written to the labelled output, in order, as (ID, Number, Type)
OUTPUT_INFO_FIELDS = [
    ('AC', 'A', 'Integer'),
    ('AF', 'A', 'Float'),
    ('AN', '1', 'Integer'),
    ('algorithms', '.', 'String'),
    ('gnomad_sv_ID', '1', 'String'),
    ('gnomad_sv_AF', '1', 'Float'),
    ('lof', '.', 'String'),
    ('n_het', '1', 'Integer'),
    ('n_homalt', '1', 'Integer'),
    ('svlen', '1', 'Integer'),
    ('svtype', '1', 'String'),
    ('status', '1', 'String'),
    ('end', '1', 'Integer'),
    ('chr2', '1', 'String'),
    ('end2', '1', 'Integer'),
    ('male_af', 'A', 'Float'),
    ('female_af', 'A', 'Float'),
    ('lof_ensg', '.', 'String'),
    ('gene_id', '1', 'String'),
    ('categorybooleansv1', '1', 'Integer'),
]


def rearrange_annotations(
    info: dict[str, Any],
    gene_mapping: dict[str, str],
    uses_af_male_spelling: bool,
) -> dict[str, Any]:
    """
    Rearrange the INFO annotations into the curated set Talos reads downstream.

    Args:
        info (dict): the raw INFO fields for one variant
        gene_mapping (dict): gene symbol to gene ID mapping
        uses_af_male_spelling (bool): whether this callset spells the per-sex
            frequencies AF_MALE/AF_FEMALE instead of MALE_AF/FEMALE_AF

    Returns:
        a new dict of the curated annotations
    """

    predicted_lof = info['PREDICTED_LOF']
    lof = set(predicted_lof.split(',')) if isinstance(predicted_lof, str) else set(predicted_lof)

    # trying to sustain backwards compatibility with a changing GATK-SV/gCNV pipeline combination
    if uses_af_male_spelling:
        male_af, female_af = info.get('AF_MALE'), info.get('AF_FEMALE')
    else:
        male_af, female_af = info.get('MALE_AF'), info.get('FEMALE_AF')

    return {
        'AC': info['AC'],
        'AF': info['AF'],
        'AN': info['AN'],
        # accept, but don't force, this GATK-SV field
        'algorithms': info.get('ALGORITHMS', 'gCNV'),
        'gnomad_sv_ID': info[f'{GNOMAD_POP}_sv_SVID'],
        'gnomad_sv_AF': info[f'{GNOMAD_POP}_sv_AF'],
        'lof': lof,
        'n_het': info['N_HET'],
        'n_homalt': info['N_HOMALT'],
        'svlen': info['SVLEN'],
        'svtype': info['SVTYPE'],
        'status': info.get('STATUS'),
        'end': info.get('END'),
        'chr2': info.get('CHR2'),
        'end2': info.get('END2'),
        'male_af': male_af,
        'female_af': female_af,
        # match the symbols to gene IDs; unmapped symbols pass through unchanged
        'lof_ensg': {gene_mapping.get(gene, gene) for gene in lof},
        'gene_id': {gene_mapping.get(gene, gene) for gene in lof},
    }


def passes_af_filter(info: dict[str, Any], af_threshold: float = 0.03) -> bool:
    """gnomAD AF below threshold; a missing AF is treated as 0 and survives."""
    return (first_value(info['gnomad_sv_AF']) or 0) < af_threshold


def passes_callset_af_filter(info: dict[str, Any], ac_threshold: float = 0.03) -> bool:
    """
    Both per-sex callset frequencies at or below threshold.
    We don't need to worry about minimum cohort size due to the minimum group size of GATK-SV.
    A missing frequency fails, matching the Hail filter this replaces.
    """
    male_af = first_value(info['male_af'])
    female_af = first_value(info['female_af'])
    return male_af is not None and female_af is not None and male_af <= ac_threshold and female_af <= ac_threshold


def diploidise_genotypes(genotypes: list[list]) -> tuple[list[list], bool]:
    """
    Recast haploid calls to a biallelic representation, going with Hom-Alt/Hom-Ref.

    if GT == 1, recast as [1, 1]
    if GT == 0, recast as [0, 0]
    missing haploid calls are left untouched

    Args:
        genotypes: cyvcf2-style genotypes, each [allele, ..., phased]

    Returns:
        the (possibly modified) genotypes, and whether any modification was made
    """

    modified = False
    result = []
    for gt in genotypes:
        alleles, phased = gt[:-1], gt[-1]
        if len(alleles) == 1 and alleles[0] >= 0:
            result.append([alleles[0], alleles[0], phased])
            modified = True
        else:
            result.append(gt)
    return result, modified


def prepare_output_header(reader: VCF) -> bool:
    """
    Add the curated INFO field definitions to the reader's header, ahead of Writer creation.

    Returns:
        whether this callset spells the per-sex frequencies AF_MALE/AF_FEMALE
        (a header-level property) instead of MALE_AF/FEMALE_AF
    """
    for field_id, number, field_type in OUTPUT_INFO_FIELDS:
        if not header_has_field(reader, field_id):
            reader.add_info_to_header(
                {'ID': field_id, 'Number': number, 'Type': field_type, 'Description': 'Talos SV annotation'},
            )
    return header_has_field(reader, 'AF_MALE')


def cli_main():
    """
    main method wrapper for console script execution
    """
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        required=True,
        help='path to input VCF',
    )
    parser.add_argument(
        '--panelapp',
        required=True,
        help='GeneratePanelData JSON',
    )
    parser.add_argument(
        '--pedigree',
        required=True,
        help='Cohort Pedigree',
    )
    parser.add_argument(
        '--output',
        help='Where to write the VCF',
        required=True,
    )
    args = parser.parse_args()
    main(
        vcf_path=args.input,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.output,
    )


def main(vcf_path: str, panelapp_path: str, pedigree: str, vcf_out: str):
    """
    Stream the annotated SV VCF, filter, apply category annotation, write as a VCF.

    Args:
        vcf_path (str): where to find vcf output
        panelapp_path ():
        pedigree ():
        vcf_out (str): where to write VCF out
    """
    logger.info(
        r"""Welcome To
 ███████████   █████████   █████          ███████     █████████
 █   ███   █  ███     ███   ███         ███     ███  ███     ███
     ███      ███     ███   ███        ███       ███ ███
     ███      ███████████   ███        ███       ███  █████████
     ███      ███     ███   ███        ███       ███         ███
     ███      ███     ███   ███      █  ███     ███  ███     ███
    █████    █████   █████ ███████████    ███████     █████████
        """,
    )

    # read the parsed panelapp data
    logger.info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path, return_model=PanelApp)
    if not isinstance(panelapp, PanelApp):
        raise TypeError(f'PanelApp was not a PanelApp object: {panelapp}')

    gene_id_mapping = get_symbol_to_ensg_mapping(panelapp)

    # pull green genes from the panelapp data - all genes in the model are green
    green_genes = set(panelapp.genes)
    logger.info(f'Extracted {len(green_genes)} green genes')

    pedigree_data = parse_pedigree(pedigree)

    reader = VCF(vcf_path)

    # subset to currently considered samples
    if not subset_reader_to_pedigree(reader, pedigree_data):
        raise ValueError('No samples shared between pedigree and VCF')

    uses_af_male_spelling = prepare_output_header(reader)

    writer = Writer(vcf_out, reader, mode='wz')

    ac_threshold = config_retrieve(['RunHailFiltering', 'callset_af_sv_recessive'])

    kept = 0
    for variant in reader:
        # remove filtered variants
        if not variant_is_pass(variant):
            continue

        # drop rows with no LOF consequences
        if not variant.INFO.get('PREDICTED_LOF'):
            continue

        raw_info = dict(variant.INFO)

        # rearrange the annotations
        info = rearrange_annotations(raw_info, gene_id_mapping, uses_af_male_spelling)

        # apply blanket filters
        if not (
            passes_callset_af_filter(info, ac_threshold=ac_threshold)
            and passes_af_filter(info, af_threshold=ac_threshold)
        ):
            continue

        # hard filter to PanelApp green genes; everything left is `SV1`
        green_gene_ids = info['gene_id'] & green_genes
        if not green_gene_ids:
            continue

        # replace the raw INFO fields with the curated set
        for key in raw_info:
            del variant.INFO[key]
        for field_id, _number, _type in OUTPUT_INFO_FIELDS:
            value = info.get(field_id)
            if value is None:
                continue
            if isinstance(value, set):
                value = ','.join(sorted(value))
            variant.INFO[field_id] = value
        variant.INFO['categorybooleansv1'] = 1

        # VCF export doesn't handle hemizygous calls - recast haploid GTs as biallelic
        new_genotypes, modified = diploidise_genotypes(variant.genotypes)
        if modified:
            variant.genotypes = new_genotypes

        # one output row per green gene this SV ablates
        for gene_id in sorted(green_gene_ids):
            variant.INFO['gene_id'] = gene_id
            writer.write_record(variant)
            kept += 1

    writer.close()
    reader.close()
    logger.info(f'Wrote {kept} labelled SV rows to {vcf_out}')


if __name__ == '__main__':
    cli_main()
