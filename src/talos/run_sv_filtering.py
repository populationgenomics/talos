"""
A streaming filtering process for labelling analysis-relevant SVs
Initially this will only contain a single category
This expects data annotated by GATK's SVAnnotate tool

CategoryBooleanSV1:
- rare
- predicted LoF in a listed gene

CategoryBooleanRegulatoryCNV:
- rare
- CNV of the same type (loss/gain) as a PanelApp region, with feature overlap of at least the required %,
  or fully contained in the region
"""

from argparse import ArgumentParser
from typing import Any

from cyvcf2 import VCF, Writer
from loguru import logger

from talos.config import config_retrieve
from talos.models import PanelApp, RegionDetail
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path
from talos.vcf_streaming import (
    first_value,
    header_has_field,
    parse_pedigree,
    subset_reader_to_pedigree,
    variant_is_pass,
)

GNOMAD_POP = config_retrieve(['RunSvFiltering', 'gnomad_population'], 'gnomad_v4.1')

# PanelApp region variant types, and the SVTYPEs which can match them
REGION_SVTYPES = {'cnv_loss': {'DEL'}, 'cnv_gain': {'DUP'}, 'cnv_both': {'DEL', 'DUP'}}

# curated INFO keys the callset AF filter reads, preferred first
SEX_AF_FIELDS = ('male_af', 'female_af')

OUTPUT_FILTER_FIELDS = [
    {'ID': 'UNRESOLVED', 'Description': 'Variant is unresolved'},
]

# INFO fields written to the labelled output, in order, as (ID, Number, Type, Description)
# these are the only INFO fields the labelled output declares; a field already in the input header keeps its definition
OUTPUT_INFO_FIELDS = [
    ('AC', 'A', 'Integer', 'Callset allele count'),
    ('AF', 'A', 'Float', 'Callset allele frequency'),
    ('AN', '1', 'Integer', 'Callset allele number'),
    ('algorithms', '.', 'String', 'Calling algorithms, from ALGORITHMS (gCNV if absent)'),
    ('gnomad_sv_ID', '1', 'String', 'gnomAD SV ID of the best qualifying match, from SVAFotate Best_gnomAD_ID'),
    ('gnomad_sv_AF', '1', 'Float', 'gnomAD SV AF, the maximum across qualifying matches, from SVAFotate Max_AF'),
    ('bnd', '.', 'String', 'Gene symbols with an exonic breakend, from PREDICTED_BREAKEND_EXONIC'),
    ('lof', '.', 'String', 'Gene symbols with a predicted loss of function, from PREDICTED_LOF'),
    ('n_het', '1', 'Integer', 'Callset heterozygous sample count, from N_HET'),
    ('n_homalt', '1', 'Integer', 'Callset hom-alt sample count, from N_HOMALT'),
    ('svlen', '1', 'Integer', 'SV length, from SVLEN'),
    ('svtype', '1', 'String', 'SV type, from SVTYPE'),
    ('status', '1', 'String', 'Variant status, from STATUS'),
    ('end', '1', 'Integer', 'SV end position, from END'),
    ('chr2', '1', 'String', 'Second contig of a BND/CTX, from CHR2'),
    ('end2', '1', 'Integer', 'Position on the second contig, from END2'),
    ('male_af', 'A', 'Float', 'Callset male AF, from MALE_AF or AF_MALE'),
    ('female_af', 'A', 'Float', 'Callset female AF, from FEMALE_AF or AF_FEMALE'),
    ('bnd_ensg', '.', 'String', 'Gene IDs with an exonic breakend'),
    ('lof_ensg', '.', 'String', 'Gene IDs with a predicted loss of function'),
    ('gene_id', '1', 'String', 'Green gene ID or PanelApp region name this row is labelled against'),
    ('categorybooleansv1', '1', 'Integer', 'Predicted loss of function in a green gene'),
    ('categorybooleanbnd', '1', 'Integer', 'Exonic breakend in a green gene'),
    ('categorybooleanregulatorycnv', '1', 'Integer', 'CNV matching a PanelApp region'),
]
OUTPUT_INFO_IDS = frozenset(field[0] for field in OUTPUT_INFO_FIELDS)


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

    predicted_lof = info.get('PREDICTED_LOF', [])
    predicted_bnd = info.get('PREDICTED_BREAKEND_EXONIC', [])
    lof = set(predicted_lof.split(',')) if isinstance(predicted_lof, str) else set(predicted_lof)
    bnd = set(predicted_bnd.split(',')) if isinstance(predicted_bnd, str) else set(predicted_bnd)

    # trying to sustain backwards compatibility with a changing GATK-SV/gCNV pipeline combination
    if uses_af_male_spelling:
        male_af, female_af = info.get('AF_MALE'), info.get('AF_FEMALE')
    else:
        male_af, female_af = info.get('MALE_AF'), info.get('FEMALE_AF')

    # AC/AF/AN, N_HET/N_HOMALT and SVLEN aren't written by every caller; absent values are omitted from the output
    return {
        'AC': info.get('AC'),
        'AF': info.get('AF'),
        'AN': info.get('AN'),
        # accept, but don't force, this GATK-SV field
        'algorithms': info.get('ALGORITHMS', 'gCNV'),
        'gnomad_sv_ID': info.get(f'{GNOMAD_POP}_sv_SVID'),
        'gnomad_sv_AF': info.get(f'{GNOMAD_POP}_sv_AF', 0.0),
        'lof': lof,
        'bnd': bnd,
        'n_het': info.get('N_HET'),
        'n_homalt': info.get('N_HOMALT'),
        'svlen': info.get('SVLEN'),
        'svtype': info['SVTYPE'],
        'status': info.get('STATUS'),
        'end': info.get('END'),
        'chr2': info.get('CHR2'),
        'end2': info.get('END2'),
        'male_af': male_af,
        'female_af': female_af,
        # match the symbols to gene IDs; unmapped symbols pass through unchanged
        'lof_ensg': {gene_mapping.get(gene, gene) for gene in lof},
        'bnd_ensg': {gene_mapping.get(gene, gene) for gene in bnd},
        'gene_id': {gene_mapping.get(gene, gene) for gene in lof | bnd},
    }


def matching_regions(
    chrom: str,
    start: int,
    end: int | None,
    svtype: str | None,
    regions: list[RegionDetail],
) -> set[str]:
    """
    Find PanelApp regions this CNV matches - same type, and either region overlap of at least the region's
    required %, or fully contained in the region (subject to a minimum percentage of region overlap)

    Args:
        chrom (str): CNV contig, with or without a chr prefix
        start (int): CNV start (POS)
        end (int | None): CNV END
        svtype (str | None): CNV SVTYPE
        regions (list[RegionDetail]): PanelApp regions

    Returns:
        names of the matched regions
    """
    if end is None:
        return set()
    # a CNV fully enclosed by a region must still cover at least this % of the region to be considered a match
    min_region_percentage = config_retrieve(['RunSvFiltering', 'regulatory_region_percentage'], 10)
    matched = set()
    for region in regions:
        if chrom.removeprefix('chr') != region.chrom or svtype not in REGION_SVTYPES.get(region.variant_type, set()):
            continue
        # the overlap must be at least the required % of the region, integer maths avoids rounding
        overlap = (min(end, region.end) - max(start, region.start)) * 100
        contained = (
            start >= region.start
            and end <= region.end
            and (end - start) * 100 >= min_region_percentage * (region.end - region.start)
        )
        if contained or overlap >= region.required_overlap * (region.end - region.start):
            matched.add(region.name)
    return matched


def passes_af_filter(info: dict[str, Any], af_threshold: float = 0.03) -> bool:
    """gnomAD AF below threshold; a missing AF is treated as 0 and survives."""
    return (first_value(info['gnomad_sv_AF']) or 0) < af_threshold


def passes_callset_af_filter(
    info: dict[str, Any],
    ac_threshold: float = 0.03,
    af_fields: tuple[str, ...] = SEX_AF_FIELDS,
) -> bool:
    """
    Every callset frequency in af_fields at or below threshold - per-sex by default, overall AF as a fallback,
    or none at all (always passes) when the callset has neither.
    We don't need to worry about minimum cohort size due to the minimum group size of GATK-SV.
    A missing frequency fails, matching the Hail filter this replaces.
    """
    return all((af := first_value(info[field])) is not None and af <= ac_threshold for field in af_fields)


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


def prepare_output_header(reader: VCF) -> tuple[bool, tuple[str, ...]]:
    """
    Add the curated INFO field definitions to the reader's header, ahead of Writer creation.
    Callset-frequency fields are optional; each one absent is warned about, as its filter is skipped.

    Returns:
        whether this callset spells the per-sex frequencies AF_MALE/AF_FEMALE
        (a header-level property) instead of MALE_AF/FEMALE_AF, and the curated fields the callset AF filter
        should read: per-sex if present, else overall AF, else none
    """
    # checked before the output fields are added, as those include AC/AF/AN
    if missing := [field for field in ('AC', 'AF', 'AN', 'N_HET', 'N_HOMALT') if not header_has_field(reader, field)]:
        logger.warning(
            f'SV VCF has no {", ".join(missing)} - omitted from the output; '
            'if AC/AF are among them, the callset AC/AF filter in ValidateMOI is skipped',
        )
    af_fields: tuple[str, ...] = SEX_AF_FIELDS
    if not (header_has_field(reader, 'MALE_AF') or header_has_field(reader, 'AF_MALE')):
        af_fields = ('AF',) if header_has_field(reader, 'AF') else ()
        logger.warning(
            'SV VCF has no MALE_AF/FEMALE_AF (or AF_MALE/AF_FEMALE) - callset AF filter '
            + ('falls back to overall AF' if af_fields else 'is skipped'),
        )

    for filter_dict in OUTPUT_FILTER_FIELDS:
        reader.add_filter_to_header(filter_dict)
    for field_id, number, field_type, description in OUTPUT_INFO_FIELDS:
        if not header_has_field(reader, field_id):
            reader.add_info_to_header(
                {'ID': field_id, 'Number': number, 'Type': field_type, 'Description': description},
            )
    return header_has_field(reader, 'AF_MALE'), af_fields


def curated_header(reader: VCF) -> str:
    """
    The reader's header with INFO definitions trimmed to OUTPUT_INFO_FIELDS, as every other INFO field is
    deleted from the records before they're written. Every non-INFO line is kept.

    The reader's own header can't be trimmed in place: it still has to parse the raw INFO fields of each record.
    """
    lines = []
    for line in reader.raw_header.splitlines():
        if line.startswith('##INFO=<ID=') and line.removeprefix('##INFO=<ID=').split(',', 1)[0] not in OUTPUT_INFO_IDS:
            continue
        lines.append(line)
    return '\n'.join(lines) + '\n'


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
    regions = list(panelapp.regions.values())
    logger.info(f'Extracted {len(green_genes)} green genes and {len(regions)} PanelApp regions')

    pedigree_data = parse_pedigree(pedigree)

    reader = VCF(vcf_path)

    # subset to currently considered samples
    if not subset_reader_to_pedigree(reader, pedigree_data):
        raise ValueError('No samples shared between pedigree and VCF')

    uses_af_male_spelling, callset_af_fields = prepare_output_header(reader)

    # the writer has its own trimmed header, so each record is re-parsed against it before writing
    writer = Writer.from_string(vcf_out, curated_header(reader), mode='wz')

    ac_threshold = config_retrieve(['RunSmallFiltering', 'callset_af_sv_recessive'])

    kept = 0
    for variant in reader:
        # remove filtered variants
        if not variant_is_pass(variant):
            continue

        region_matches = matching_regions(
            variant.CHROM, variant.POS, variant.INFO.get('END'), variant.INFO.get('SVTYPE'), regions
        )

        # drop rows with no LOF consequences or region matches
        if not (variant.INFO.get('PREDICTED_LOF') or variant.INFO.get('PREDICTED_BREAKEND_EXONIC') or region_matches):
            continue

        raw_info = dict(variant.INFO)

        # rearrange the annotations
        info = rearrange_annotations(raw_info, gene_id_mapping, uses_af_male_spelling)

        # apply blanket filters
        if not (
            passes_callset_af_filter(info, ac_threshold=ac_threshold, af_fields=callset_af_fields)
            and passes_af_filter(info, af_threshold=ac_threshold)
        ):
            continue

        # hard filter to PanelApp green genes or matched PanelApp regions
        green_gene_ids = info['gene_id'] & green_genes
        if not (green_gene_ids or region_matches):
            continue

        # replace the raw INFO fields with the curated set
        for key in raw_info:
            del variant.INFO[key]

        for field_id, *_ in OUTPUT_INFO_FIELDS:
            value = info.get(field_id)
            if value is None:
                continue
            if isinstance(value, set):
                value = ','.join(sorted(value))
            variant.INFO[field_id] = value

        # VCF export doesn't handle hemizygous calls - recast haploid GTs as biallelic
        new_genotypes, modified = diploidise_genotypes(variant.genotypes)
        if modified:
            variant.genotypes = new_genotypes

        # one output row per green gene this SV ablates, and one per PanelApp region it matches
        # region rows only carry the region category, so they're not paired with gene variants as comp-hets
        for gene_id in sorted(green_gene_ids | region_matches):
            is_region = gene_id in region_matches
            variant.INFO['gene_id'] = gene_id
            # tried and tested - SVs predicted to cause LOF by SvAnnotate
            variant.INFO['categorybooleansv1'] = int(not is_region and bool(info['lof']))
            # second category introduced - exonic breakends
            variant.INFO['categorybooleanbnd'] = int(not is_region and bool(info['bnd']))
            variant.INFO['categorybooleanregulatorycnv'] = int(is_region)
            writer.write_record(writer.variant_from_string(str(variant).rstrip('\n')))
            kept += 1

    writer.close()
    reader.close()
    logger.info(f'Wrote {kept} labelled SV rows to {vcf_out}')


if __name__ == '__main__':
    cli_main()
