"""
A streaming filtering and labelling process for small variants - the Stage E replacement for
run_hail_filtering.py and annotation_scripts/annotated_vcf_into_matrixtable.py combined.

This is a cyvcf2 process - no Hail, no Spark. One variant is held in memory at a time.

Input contract (per shard):
- a normalised, BCSQ-annotated VCF (bcftools csq + echtvar, from the annotation workflow)
- fresh ClinvArbitration decisions already applied to the VCF by an `echtvar annotate` stream,
  as the INFO fields clinvar_significance / clinvar_stars / clinvar_allele
- PM5 data as a JSON map of {"TRANSCRIPT::CODON": "alleleID::stars+alleleID::stars..."}
- MANE data as the JSON written by parse_mane_into_json.py
- AC/AN/AF must be present (joint-called data carries them; otherwise run
  `bcftools +fill-tags -- -t AC,AN,AF` upstream). AF alone is derived from AC/AN here.

Output: a labelled VCF shard, one row per (variant, green gene), gathered downstream
with `bcftools concat` and consumed by ValidateMOI unchanged.
"""

import gzip
import json
from argparse import ArgumentParser
from dataclasses import dataclass, field
from typing import Any

from cyvcf2 import VCF, Variant, Writer
from loguru import logger
from mendelbrot.pedigree_parser import PedigreeParser

from talos.config import config_retrieve
from talos.models import PanelApp
from talos.utils import get_symbol_to_ensg_mapping, read_json_from_path
from talos.vcf_streaming import (
    MISSING_STRING,
    PATHOGENIC,
    consequences_to_csq_string,
    first_value,
    header_has_field,
    normalise_chrom,
    parse_bcsq_entries,
    parse_pedigree,
    split_csq_header,
    subset_reader_to_pedigree,
    variant_is_pass,
)

ADDITIONAL_CSQ_DEFAULT = ['missense', 'inframe_deletion', 'inframe_insertion']
CRITICAL_CSQ_DEFAULT = [
    'frameshift',
    'splice_acceptor',
    'splice_donor',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'transcript_ablation',
]

# cyvcf2 gt_types values (default gts012=False)
HOM_REF = 0
HET = 1
UNKNOWN = 2
HOM_ALT = 3

# regions used by the de novo configuration test
AUTOSOME_OR_PAR = 'autosome_or_par'
X_NONPAR = 'x_nonpar'
Y_NONPAR = 'y_nonpar'
MITO = 'mito'

AUTOSOMES = {str(i) for i in range(1, 23)}

# GRCh38 pseudo-autosomal regions, as Hail defines them: start inclusive, end exclusive
PAR_REGIONS = {
    'X': ((10001, 2781479), (155701383, 156030895)),
    'Y': ((10001, 2781479), (56887903, 57217415)),
}

# raw INFO names written by the gnomAD echtvar zip, and the names Talos reads downstream
GNOMAD_SOURCE_FIELDS = {
    'gnomad_AC': 'gnomad_AC_joint',
    'gnomad_AF': 'gnomad_AF_joint',
    'gnomad_AC_XY': 'gnomad_AC_joint_XY',
    'gnomad_HomAlt': 'gnomad_HomAlt_joint',
}

# per-transcript AlphaMissense INFO fields, folded into the csq string then dropped
AM_INFO_FIELDS = ('am_transcript', 'am_class', 'am_score')

# INFO fields added to the labelled output
NEW_INFO_HEADERS = [
    {'ID': 'clinvar_significance', 'Number': '1', 'Type': 'String', 'Description': 'ClinvArbitration significance'},
    {'ID': 'clinvar_stars', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration gold stars'},
    {'ID': 'clinvar_allele', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration allele ID'},
    {'ID': 'clinvar_talos', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinvArbitration Pathogenic flag'},
    {'ID': 'categorybooleanclinvarplp', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinVar P/LP, 1+ stars'},
    {'ID': 'categorybooleanclinvar0star', 'Number': '1', 'Type': 'Integer', 'Description': 'ClinVar P/LP, 0 stars'},
    {
        'ID': 'categorybooleanclinvar0starnewgene',
        'Number': '1',
        'Type': 'Integer',
        'Description': 'ClinVar P/LP, 0 stars, in a PanelApp new gene',
    },
    {'ID': 'categorybooleanalphamissense', 'Number': '1', 'Type': 'Integer', 'Description': 'AlphaMissense P/LP'},
    {'ID': 'categorybooleanhighimpact', 'Number': '1', 'Type': 'Integer', 'Description': 'Critical consequence'},
    {'ID': 'categorybooleanspliceai', 'Number': '1', 'Type': 'Integer', 'Description': 'Significant SpliceAI score'},
    {'ID': 'categorybooleanavi', 'Number': '1', 'Type': 'Integer', 'Description': 'Significant AVI score'},
    {'ID': 'categorysampledenovo', 'Number': '1', 'Type': 'String', 'Description': 'Samples with a de novo call'},
    {'ID': 'categorydetailspm5', 'Number': '1', 'Type': 'String', 'Description': 'ClinVar alleles at this codon'},
    {'ID': 'gnomad_AC', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD joint AC'},
    {'ID': 'gnomad_AF', 'Number': '1', 'Type': 'Float', 'Description': 'gnomAD joint AF'},
    {'ID': 'gnomad_AC_XY', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD joint AC XY'},
    {'ID': 'gnomad_HomAlt', 'Number': '1', 'Type': 'Integer', 'Description': 'gnomAD joint Hom-Alt count'},
    {'ID': 'gene_id', 'Number': '1', 'Type': 'String', 'Description': 'Green gene this row is labelled against'},
]


@dataclass
class DeNovoConfig:
    """Thresholds for the naive de novo configuration test, from RunHailFiltering.de_novo."""

    min_depth: int = 5
    max_depth: int = 1000
    min_proband_gq: int = 25
    min_all_sample_gq: int = 19
    apply_min_all_sample_gq: bool = True

    @classmethod
    def from_config(cls) -> 'DeNovoConfig':
        conf = config_retrieve(['RunHailFiltering', 'de_novo'], {})
        return cls(
            min_depth=conf.get('min_depth', 5),
            max_depth=conf.get('max_depth', 1000),
            min_proband_gq=conf.get('min_proband_gq', 25),
            min_all_sample_gq=conf.get('min_all_sample_gq', 19),
            apply_min_all_sample_gq=conf.get('apply_min_all_sample_gq', True),
        )


@dataclass
class Trio:
    """A complete trio - all three members present in both pedigree and VCF."""

    child_id: str
    child_idx: int
    father_idx: int
    mother_idx: int
    # True/False from the pedigree, None for unknown sex
    child_is_female: bool | None
    child_affected: bool
    father_affected: bool
    mother_affected: bool


@dataclass
class StreamingContext:
    """Everything label_variant needs, resolved once before the streaming loop."""

    writer: Writer
    csq_fields: list[str]
    mane: dict[str, dict[str, str]]
    symbol_to_ensg: dict[str, str]
    green_genes: set[str]
    new_genes: set[str]
    pm5: dict[str, str] | None
    # None when de novo testing is disabled for this run
    trios: list[Trio] | None
    dn_conf: DeNovoConfig
    af_semi_rare: float
    am_threshold: float
    spliceai_threshold: float
    avi_threshold: float | None
    critical_csqs: set[str]
    dn_relevant_csqs: set[str]
    derive_af: bool
    has_spliceai: bool
    has_avi: bool
    kept: int = field(default=0)


def load_json_dict(path: str) -> dict:
    """Read a JSON file (optionally gzipped) into a dict."""
    open_method = gzip.open if path.endswith('.gz') else open
    with open_method(path, 'rt') as handle:
        return json.load(handle)


def load_pm5_json(path: str | None) -> dict[str, str] | None:
    """
    Read the PM5 JSON: {"TRANSCRIPT::CODON": "alleleID::stars+alleleID::stars"}.
    List values are tolerated and joined with '+', matching the ClinvArbitration TSV/HT content.
    """
    if path is None:
        logger.info('PM5 not required or requested, skipping annotation')
        return None
    data = load_json_dict(path)
    logger.info(f'Read {len(data)} PM5 residue entries from {path}')
    return {key: '+'.join(sorted(value)) if isinstance(value, list) else value for key, value in data.items()}


def check_frequency_fields(reader: VCF) -> bool:
    """
    Port of populate_callset_frequencies: require AC/AN/AF, deriving AF from AC & AN if absent.

    Whole-callset frequency generation (the Hail variant_qc fallback) is not possible in a
    per-shard stream - that job belongs to `bcftools +fill-tags` upstream.

    Returns:
        whether AF needs deriving per-variant
    """
    present = {fid for fid in ('AC', 'AN', 'AF') if header_has_field(reader, fid)}
    if present == {'AC', 'AN', 'AF'}:
        logger.info('AC, AN, AF already present, skipping annotation')
        return False
    if {'AC', 'AN'} <= present:
        logger.info('AC, AN present, deriving AF from existing annotations')
        reader.add_info_to_header(
            {'ID': 'AF', 'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency, derived from AC/AN'},
        )
        return True
    raise ValueError(
        'AC/AN INFO fields are missing from this VCF. '
        'Populate them upstream with: bcftools +fill-tags in.vcf.gz -- -t AC,AN,AF',
    )


def require_gnomad_fields(reader: VCF):
    """The gnomAD echtvar annotations are a hard requirement, matching the Hail schema access."""
    if missing := [src for src in GNOMAD_SOURCE_FIELDS.values() if not header_has_field(reader, src)]:
        raise ValueError(f'Required gnomAD INFO fields are missing from this VCF: {missing}')


def structure_consequences(
    variant: Variant,
    csq_fields: list[str],
    mane: dict[str, dict[str, str]],
    symbol_to_ensg: dict[str, str],
) -> list[dict[str, Any]]:
    """
    Parse the BCSQ INFO field into consequence dicts, then fold in the AlphaMissense INFO
    fields (per matching transcript), MANE annotations, and symbol-to-ENSG gene IDs.

    Port of annotated_vcf_into_matrixtable's csq_strings_into_hail_structs +
    annotate_all_transcript_consequences.
    """
    bcsq = variant.INFO.get('BCSQ')
    if not bcsq:
        return []

    am_transcript = variant.INFO.get('am_transcript')
    am_class = variant.INFO.get('am_class', '')
    am_score = variant.INFO.get('am_score', 0.0)

    consequences = parse_bcsq_entries(bcsq, csq_fields)
    for csq in consequences:
        transcript = csq['transcript']
        am_match = transcript == am_transcript
        csq['am_class'] = am_class if am_match else ''
        csq['am_pathogenicity'] = am_score if am_match else 0.0

        mane_entry = mane.get(transcript, {})
        csq['mane_status'] = mane_entry.get('mane_status', '')
        csq['ensp'] = mane_entry.get('ensp', '')
        csq['mane_id'] = mane_entry.get('mane_id', '')

        # unmapped symbols pass through unchanged
        csq['gene_id'] = symbol_to_ensg.get(csq['gene'], csq['gene'])
    return consequences


def read_clinvar_fields(info: Any) -> tuple[str, int, int]:
    """
    Normalise the echtvar-applied ClinvArbitration INFO fields.
    Absent annotations (or echtvar missing-value sentinels) become 'missing'/0.
    """
    significance = info.get('clinvar_significance')
    if significance in (None, '', '.'):
        significance = MISSING_STRING
    # VCF INFO values cannot contain spaces, so the echtvar zip encodes them as underscores
    significance = significance.replace('_', ' ')
    stars = info.get('clinvar_stars')
    stars = int(stars) if stars is not None and stars >= 0 else 0
    allele = info.get('clinvar_allele')
    allele = int(allele) if allele is not None and allele >= 0 else 0
    return significance, stars, allele


def is_confident_benign(significance: str, stars: int) -> bool:
    """Confidently benign variants are removed outright (annotate_clinvarbitration)."""
    return 'benign' in significance.lower() and stars > 0


def clinvar_category_flags(significance: str, stars: int, gene_ids: set[str], new_genes: set[str]) -> dict[str, int]:
    """The ClinVar-derived flags: talos (P/LP any stars), plp (1+ stars), 0star, 0star-in-new-gene."""
    pathogenic = significance == PATHOGENIC
    return {
        'clinvar_talos': int(pathogenic),
        'categorybooleanclinvarplp': int(pathogenic and stars > 0),
        'categorybooleanclinvar0star': int(pathogenic and stars == 0),
        'categorybooleanclinvar0starnewgene': int(pathogenic and stars == 0 and bool(gene_ids & new_genes)),
    }


def passes_population_rare(gnomad_af: float | None, clinvar_talos: int, af_threshold: float) -> bool:
    """gnomAD AF below threshold; missing is treated as 0; ClinVar Pathogenic overrides."""
    return bool(clinvar_talos) or (gnomad_af or 0.0) < af_threshold


def passes_ac_filter(
    info: Any,
    clinvar_talos: int,
    af_threshold: float = 0.01,
    min_ac_to_filter: int = 5,
) -> bool:
    """
    Callset-frequency filter (filter_matrix_by_ac): never removes variants with few instances,
    and is overridden by a ClinVar Pathogenic annotation.
    """
    if clinvar_talos:
        return True
    ac = first_value(info.get('AC')) or 0
    af = first_value(info.get('AF')) or 0.0
    return ac <= min_ac_to_filter or af < af_threshold


def consequences_for_gene(consequences: list[dict], gene_id: str) -> list[dict]:
    """
    The per-gene consequence reduction from split_rows_by_gene_and_filter_to_green:
    this gene only, and protein_coding / snRNA / MANE-transcript biotypes only.
    """
    return [
        csq
        for csq in consequences
        if csq['gene_id'] == gene_id
        and (csq['biotype'] in ('protein_coding', 'snRNA') or 'NM' in csq.get('mane_id', ''))
    ]


def filter_by_consequence(consequences: list[dict], relevant_csqs: set[str]) -> list[dict]:
    """Reduce consequences to critical/additional terms or snRNA - the de novo candidate gate."""
    return [
        csq for csq in consequences if relevant_csqs & set(csq['consequence'].split('&')) or csq['biotype'] == 'snRNA'
    ]


def category_alphamissense(consequences: list[dict], am_threshold: float) -> int:
    """AlphaMissense likely-pathogenic on at least one transcript."""
    return int(any((csq.get('am_pathogenicity') or 0.0) >= am_threshold for csq in consequences))


def category_high_impact(consequences: list[dict], critical_csqs: set[str]) -> int:
    """Critical protein consequence on at least one transcript."""
    return int(any(critical_csqs & set(csq['consequence'].split('&')) for csq in consequences))


def category_spliceai(delta_score: float | None, threshold: float) -> int:
    """Significant SpliceAI delta score - only ever set in deployments carrying the annotation."""
    return int(delta_score is not None and delta_score >= threshold)


def category_avi(avi_score: float | None, threshold: float | None) -> int:
    """Significant AVI score - only ever set in deployments carrying the annotation."""
    return int(threshold is not None and avi_score is not None and avi_score >= threshold)


def pm5_matches(is_snv: bool, consequences: list[dict], pm5: dict[str, str] | None) -> str:
    """
    PM5: other ClinVar-pathogenic missense variants affecting the same residue.
    Missense SNVs only; exact-allele matches are filtered out downstream.
    """
    if pm5 is None or not is_snv:
        return MISSING_STRING
    hits = set()
    for csq in consequences:
        if 'missense' not in csq['consequence'] or csq['codon'] is None:
            continue
        if alleles := pm5.get(f'{csq["transcript"]}::{csq["codon"]}'):
            hits.add(alleles)
    return '+'.join(sorted(hits)) if hits else MISSING_STRING


def variant_region(chrom: str, pos: int) -> str | None:
    """Classify a locus for the de novo configuration test, mirroring Hail's locus methods."""
    contig = normalise_chrom(chrom)
    if contig in AUTOSOMES:
        return AUTOSOME_OR_PAR
    if contig in PAR_REGIONS:
        if any(start <= pos < end for start, end in PAR_REGIONS[contig]):
            return AUTOSOME_OR_PAR
        return X_NONPAR if contig == 'X' else Y_NONPAR
    if contig == 'M':
        return MITO
    return None


def assemble_trios(pedigree_data: PedigreeParser, vcf_samples: list[str]) -> list[Trio]:
    """Collect the complete trios - child and both parents present in the (subset) VCF."""
    sample_idx = {sample: idx for idx, sample in enumerate(vcf_samples)}
    trios = []
    for participant in pedigree_data.participants.values():
        if not (
            participant.sample_id in sample_idx
            and participant.father_id in sample_idx
            and participant.mother_id in sample_idx
        ):
            continue
        father = pedigree_data.participants[participant.father_id]
        mother = pedigree_data.participants[participant.mother_id]
        is_female = True if participant.is_female else False if participant.is_male else None
        trios.append(
            Trio(
                child_id=participant.sample_id,
                child_idx=sample_idx[participant.sample_id],
                father_idx=sample_idx[participant.father_id],
                mother_idx=sample_idx[participant.mother_id],
                child_is_female=is_female,
                child_affected=participant.is_affected,
                father_affected=father.is_affected,
                mother_affected=mother.is_affected,
            ),
        )
    logger.info(f'Assembled {len(trios)} complete trios for de novo testing')
    return trios


def entry_depths(variant: Variant, n_samples: int, default_depth: int) -> list[int]:
    """Per-sample depths: DP, falling back to summed AD, falling back to a passing default."""
    if 'DP' in variant.FORMAT:
        return [int(value[0]) if value[0] >= 0 else default_depth for value in variant.format('DP')]
    if 'AD' in variant.FORMAT:
        return [int(row[row >= 0].sum()) if (row >= 0).any() else default_depth for row in variant.format('AD')]
    return [default_depth] * n_samples


def resolve_trio_entry(
    gt_type: int,
    gq: float,
    depth: int,
    affected: bool,
    conf: DeNovoConfig,
) -> tuple[int | None, float | None]:
    """
    The entry-level adjustments and filters from annotate_category_de_novo, for one sample:
    - remove entries below the all-sample GQ floor (when configured)
    - recover missing GTs as HomRef when the GQ is convincing
    - remove entries with out-of-bounds depth
    - remove affected members' entries at or below the proband GQ floor

    Returns:
        the effective (gt_type, gq) - (None, None) for a removed entry
    """
    gq_value: float | None = None if gq < 0 else float(gq)
    gt: int | None = None if gt_type == UNKNOWN else gt_type

    if conf.apply_min_all_sample_gq and gq_value is not None and gq_value < conf.min_all_sample_gq:
        return None, None
    if gt is None and gq_value is not None and gq_value > conf.min_all_sample_gq:
        gt = HOM_REF
    if depth < conf.min_depth or depth > conf.max_depth:
        return None, None
    if affected and gq_value is not None and gq_value <= conf.min_proband_gq:
        return None, None
    return gt, gq_value


def candidate_configuration(
    region: str,
    is_female: bool | None,
    kid: int | None,
    dad: int | None,
    mom: int | None,
) -> bool:
    """The Mendelian-inconsistency configurations accepted as candidate de novo calls."""
    parents_wt = dad == HOM_REF and mom == HOM_REF
    if region == AUTOSOME_OR_PAR:
        return kid == HET and parents_wt
    if region == X_NONPAR:
        # unknown child sex cannot be tested on X
        if is_female is None:
            return False
        return (kid == HET if is_female else kid in (HET, HOM_ALT)) and parents_wt
    if region == Y_NONPAR:
        return kid in (HET, HOM_ALT) and dad == HOM_REF
    return region == MITO and kid == HOM_ALT and mom == HOM_REF


def de_novo_sample_ids(variant: Variant, trios: list[Trio], conf: DeNovoConfig) -> list[str]:
    """Run the naive de novo configuration test for every complete trio at this variant."""
    region = variant_region(variant.CHROM, variant.POS)
    if region is None or not trios:
        return []

    gt_types = variant.gt_types
    gqs = variant.gt_quals
    depths = entry_depths(variant, len(gt_types), default_depth=conf.min_depth + 1)

    hits = []
    for trio in trios:
        kid_gt, kid_gq = resolve_trio_entry(
            gt_types[trio.child_idx],
            gqs[trio.child_idx],
            depths[trio.child_idx],
            trio.child_affected,
            conf,
        )
        dad_gt, _dad_gq = resolve_trio_entry(
            gt_types[trio.father_idx],
            gqs[trio.father_idx],
            depths[trio.father_idx],
            trio.father_affected,
            conf,
        )
        mom_gt, _mom_gq = resolve_trio_entry(
            gt_types[trio.mother_idx],
            gqs[trio.mother_idx],
            depths[trio.mother_idx],
            trio.mother_affected,
            conf,
        )
        if not candidate_configuration(region, trio.child_is_female, kid_gt, dad_gt, mom_gt):
            continue
        if kid_gq is None or kid_gq < conf.min_proband_gq:
            continue
        hits.append(trio.child_id)
    return sorted(hits)


def any_category_assigned(categories: dict[str, Any]) -> bool:
    """filter_to_categorised: at least one boolean flag set, or a sample/details value present."""
    for key, value in categories.items():
        if key.startswith('categoryboolean') and value == 1:
            return True
        if key.startswith(('categorysample', 'categorydetails')) and value != MISSING_STRING:
            return True
    return False


def prepare_output_header(reader: VCF):
    """Add the labelled-output INFO definitions ahead of Writer creation."""
    for header_entry in NEW_INFO_HEADERS:
        if not header_has_field(reader, header_entry['ID']):
            reader.add_info_to_header(header_entry)
    csq_contents = '|'.join(config_retrieve(['RunHailFiltering', 'csq_string']))
    if not header_has_field(reader, 'csq'):
        reader.add_info_to_header(
            {'ID': 'csq', 'Number': '.', 'Type': 'String', 'Description': f'Format: {csq_contents}'},
        )


def scrub_raw_info(variant: Variant):
    """Remove the INFO fields we re-shape: BCSQ, raw gnomAD names, and per-transcript AlphaMissense."""
    for key in [key for key, _value in variant.INFO]:
        if key == 'BCSQ' or key.startswith('gnomad_') or key in AM_INFO_FIELDS:
            del variant.INFO[key]


def write_gene_rows(
    variant: Variant,
    ctx: StreamingContext,
    consequences: list[dict],
    green_hits: list[str],
    base_flags: dict[str, int],
) -> int:
    """Apply the per-gene categories, and write one labelled row per categorised green gene."""
    written = 0
    is_snv = len(variant.REF) == 1 and len(variant.ALT[0]) == 1
    clinvar_talos = base_flags['clinvar_talos']

    # variant-level categories, shared by every gene row
    spliceai_flag = (
        category_spliceai(variant.INFO.get('splice_ai_delta'), ctx.spliceai_threshold) if ctx.has_spliceai else 0
    )
    avi_flag = category_avi(variant.INFO.get('avi_score'), ctx.avi_threshold) if ctx.has_avi else 0

    # the trio test result is identical for every gene - compute at most once
    denovo_cache: list[str] | None = None

    for gene_id in green_hits:
        gene_csqs = consequences_for_gene(consequences, gene_id)

        categories: dict[str, Any] = dict(base_flags)
        categories['categorybooleanalphamissense'] = category_alphamissense(gene_csqs, ctx.am_threshold)
        categories['categorybooleanhighimpact'] = category_high_impact(gene_csqs, ctx.critical_csqs)
        categories['categorybooleanspliceai'] = spliceai_flag
        categories['categorybooleanavi'] = avi_flag
        categories['categorydetailspm5'] = pm5_matches(is_snv, gene_csqs, ctx.pm5)

        # only rows with a relevant consequence (or a ClinVar Pathogenic rating) are trio-tested
        dn_value = MISSING_STRING
        if ctx.trios is not None and (clinvar_talos or filter_by_consequence(gene_csqs, ctx.dn_relevant_csqs)):
            if denovo_cache is None:
                denovo_cache = de_novo_sample_ids(variant, ctx.trios, ctx.dn_conf)
            if denovo_cache:
                dn_value = ','.join(denovo_cache)
        categories['categorysampledenovo'] = dn_value

        if not any_category_assigned(categories):
            continue

        for key, value in categories.items():
            variant.INFO[key] = value
        if gene_csqs:
            variant.INFO['csq'] = consequences_to_csq_string(gene_csqs)
        variant.INFO['gene_id'] = gene_id
        ctx.writer.write_record(variant)
        written += 1
    return written


def label_variant(variant: Variant, ctx: StreamingContext) -> int:
    """
    Run one variant through the full filter/label sequence.
    INFO-level predicates run first; FORMAT (genotype) data is only touched by the trio test.

    Returns:
        the number of labelled rows written for this variant
    """
    # star alleles are not currently handled
    if not variant.ALT or '*' in variant.ALT:
        return 0

    significance, stars, allele = read_clinvar_fields(variant.INFO)
    if is_confident_benign(significance, stars):
        return 0
    clinvar_talos = int(significance == PATHOGENIC)

    gnomad_af = variant.INFO.get(GNOMAD_SOURCE_FIELDS['gnomad_AF'])
    if not (
        passes_population_rare(gnomad_af, clinvar_talos, ctx.af_semi_rare)
        and (variant_is_pass(variant) or clinvar_talos)
    ):
        return 0

    if ctx.derive_af:
        ac, an = first_value(variant.INFO.get('AC')) or 0, variant.INFO.get('AN') or 0
        variant.INFO['AF'] = float(ac) / an if an else 0.0
    if not passes_ac_filter(variant.INFO, clinvar_talos):
        return 0

    consequences = structure_consequences(variant, ctx.csq_fields, ctx.mane, ctx.symbol_to_ensg)
    gene_ids = {csq['gene_id'] for csq in consequences}
    green_hits = sorted(gene_ids & ctx.green_genes)
    if not green_hits:
        return 0

    base_flags = clinvar_category_flags(significance, stars, gene_ids, ctx.new_genes)

    # take the renamed gnomAD values before scrubbing the raw fields
    gnomad_values = {new: variant.INFO.get(src) or 0 for new, src in GNOMAD_SOURCE_FIELDS.items()}
    scrub_raw_info(variant)
    variant.INFO['clinvar_significance'] = significance
    variant.INFO['clinvar_stars'] = stars
    variant.INFO['clinvar_allele'] = allele
    for key, value in gnomad_values.items():
        variant.INFO[key] = value

    return write_gene_rows(variant, ctx, consequences, green_hits, base_flags)


def de_novo_is_disabled() -> bool:
    """De novo testing is skipped for ignored-category and singleton runs."""
    ignored_categories = config_retrieve(['ValidateMOI', 'ignore_categories'], [])
    if any(to_ignore in ignored_categories for to_ignore in ('de_novo', 'denovo', '4')):
        return True
    return bool(config_retrieve('singletons', False))


def cli_main():
    parser = ArgumentParser(description='Streaming filter and label for annotated small-variant VCF shards')
    parser.add_argument('--input', help='Path to the annotated VCF shard', required=True)
    parser.add_argument('--panelapp', help='PanelApp JSON', required=True)
    parser.add_argument('--pedigree', help='Cohort pedigree', required=True)
    parser.add_argument('--output', help='Where to write the labelled VCF, written bgzipped', required=True)
    parser.add_argument('--mane', help='MANE JSON (parse_mane_into_json.py)', required=True)
    parser.add_argument('--pm5', help='JSON of ClinVar PM5 annotations, optional', default=None)
    args = parser.parse_args()
    main(
        vcf_path=args.input,
        panelapp_path=args.panelapp,
        pedigree=args.pedigree,
        vcf_out=args.output,
        mane_path=args.mane,
        pm5=args.pm5,
    )


def main(
    vcf_path: str,
    panelapp_path: str,
    pedigree: str,
    vcf_out: str,
    mane_path: str,
    *,
    pm5: str | None = None,
):
    """
    Stream the annotated VCF shard, filter, apply category labels, write a labelled VCF.

    Args:
        vcf_path (str): path to the annotated (BCSQ + echtvar, ClinVar-applied) VCF shard
        panelapp_path (str): GeneratePanelData JSON
        pedigree (str): path to the cohort pedigree
        vcf_out (str): where to write the labelled VCF
        mane_path (str): path to the MANE JSON
        pm5 (str | None): path to the PM5 JSON, optional
    """
    logger.info(
        r"""Welcome To
███████████   █████████   █████          ███████     █████████
█   ███   █  ███     ███   ███         ███     ███  ███     ███
    ███      ███     ███   ███        ███       ███ ███
    ███      ███████████   ███        ███       ███  █████████
    ███      ███     ███   ███        ███       ███         ███
    ███      ███     ███   ███      █  ███     ███  ███     ███
   █████    █████   █████ ███████████    ███████     █████████ """,
    )

    logger.info(f'Reading PanelApp data from {panelapp_path!r}')
    panelapp = read_json_from_path(panelapp_path, return_model=PanelApp)
    if not isinstance(panelapp, PanelApp):
        raise TypeError(f'PanelApp was not a PanelApp object: {panelapp}')

    green_genes = set(panelapp.genes)
    logger.info(f'Extracted {len(green_genes)} green genes')
    new_genes = {ensg for ensg, details in panelapp.genes.items() if getattr(details, 'new', None)}
    logger.info(f'Extracted {len(new_genes)} new genes')
    symbol_to_ensg = get_symbol_to_ensg_mapping(panelapp)

    pedigree_data = parse_pedigree(pedigree)

    reader = VCF(vcf_path)
    if not subset_reader_to_pedigree(reader, pedigree_data):
        raise ValueError('No samples shared between pedigree and VCF')

    csq_fields = split_csq_header(reader)
    derive_af = check_frequency_fields(reader)
    require_gnomad_fields(reader)
    prepare_output_header(reader)

    writer = Writer(vcf_out, reader, mode='wz')

    trios = None
    if de_novo_is_disabled():
        logger.info('Skipping de novo annotation, category 4 will not be used during this analysis')
    else:
        trios = assemble_trios(pedigree_data, reader.samples)

    critical_csqs = set(config_retrieve(['RunHailFiltering', 'critical_csq'], CRITICAL_CSQ_DEFAULT))
    additional_csqs = set(config_retrieve(['RunHailFiltering', 'additional_csq'], ADDITIONAL_CSQ_DEFAULT))
    ctx = StreamingContext(
        writer=writer,
        csq_fields=csq_fields,
        mane=load_json_dict(mane_path),
        symbol_to_ensg=symbol_to_ensg,
        green_genes=green_genes,
        new_genes=new_genes,
        pm5=load_pm5_json(pm5),
        trios=trios,
        dn_conf=DeNovoConfig.from_config(),
        af_semi_rare=config_retrieve(['RunHailFiltering', 'af_semi_rare'], 0.01),
        am_threshold=config_retrieve(['RunHailFiltering', 'am_pathogenicity'], 0.564),
        spliceai_threshold=config_retrieve(['RunHailFiltering', 'spliceai'], 0.5),
        avi_threshold=config_retrieve(['RunHailFiltering', 'avi'], None),
        critical_csqs=critical_csqs,
        dn_relevant_csqs=critical_csqs | additional_csqs,
        derive_af=derive_af,
        has_spliceai=header_has_field(reader, 'splice_ai_delta'),
        has_avi=header_has_field(reader, 'avi_score'),
    )

    for variant in reader:
        ctx.kept += label_variant(variant, ctx)

    writer.close()
    reader.close()
    logger.info(f'Wrote {ctx.kept} labelled rows to {vcf_out}')

    if not ctx.kept:
        logger.warning('No variants were labelled in this shard')


if __name__ == '__main__':
    cli_main()
