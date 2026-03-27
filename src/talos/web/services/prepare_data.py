"""
Prepare ResultData for JSON API responses.

Converts Pydantic models into plain JSON-serializable dicts
that the static HTML frontend can consume directly.
"""

from talos.models import ReportVariant, ResultData, SmallVariant, StructuralVariant

MAX_INDEL_LEN = 10


def get_var_change(rv: ReportVariant) -> str:
    """Compute a display string for the variant change."""
    var = rv.var_data
    ref = var.coordinates.ref
    alt = var.coordinates.alt

    if isinstance(var, SmallVariant):
        if len(ref) > MAX_INDEL_LEN or len(alt) > MAX_INDEL_LEN:
            ref_len, alt_len = len(ref), len(alt)
            if ref_len > alt_len:
                return f'del {ref_len - alt_len}bp'
            if ref_len == alt_len:
                return f'complex delins {ref_len}bp'
            return f'ins {alt_len - ref_len}bp'
        return f'{ref}->{alt}'

    if isinstance(var, StructuralVariant):
        return f'{var.info.get("svtype", "SV")} {var.info.get("svlen", "?")}bp'

    return f'{ref}->{alt}'


def parse_mane_csq(rv: ReportVariant) -> tuple[str, str]:
    """Extract MANE consequence and HGVS protein notation from a SmallVariant."""
    if not isinstance(rv.var_data, SmallVariant):
        return ('', '')

    mane_consequences: set[str] = set()
    mane_hgvsps: set[str] = set()

    for csq in rv.var_data.transcript_consequences:
        if 'consequence' not in csq:
            continue
        if csq.get('mane_id'):
            mane_consequences.update(str(csq['consequence']).split('&'))
            if aa := csq.get('amino_acid_change'):
                mane_hgvsps.add(f'{csq.get("ensp", "")}: {aa}')

    csq_str = ', '.join(c.replace('_variant', '').replace('_', ' ') for c in mane_consequences)
    hgvsp_str = ', '.join(mane_hgvsps)
    return (csq_str, hgvsp_str)


def prepare_variant(rv: ReportVariant, decision: dict | None = None) -> dict:
    """Prepare a single ReportVariant as a JSON-serializable dict."""
    coords = rv.var_data.coordinates
    var_data_dump = rv.var_data.model_dump(mode='json')

    # compute AlphaMissense max
    am_scores = []
    if isinstance(rv.var_data, SmallVariant):
        am_scores = [
            float(csq['am_pathogenicity']) for csq in rv.var_data.transcript_consequences if csq.get('am_pathogenicity')
        ]
    var_data_dump.setdefault('info', {})['alpha_missense_max'] = max(am_scores) if am_scores else 'missing'

    mane_csq, mane_hgvsps = parse_mane_csq(rv)

    # panel match display
    pheno_matches = sorted(f'{name}({pid})' for pid, name in rv.panels.matched.items())
    forced_matches = sorted(f'{name}({pid})' for pid, name in rv.panels.forced.items())

    # transcript consequences
    tx_csq = list(rv.var_data.transcript_consequences) if isinstance(rv.var_data, SmallVariant) else []

    # gene list as [gene_id, gene_id] pairs
    genes = [[g.strip(), g.strip()] for g in rv.gene.split(',')]

    return {
        'sample': rv.sample,
        'var_type': rv.var_data.__class__.__name__,
        'chrom': coords.chrom,
        'pos': coords.pos,
        'ref': coords.ref,
        'alt': coords.alt,
        'change': get_var_change(rv),
        'categories': dict(rv.categories),
        'first_tagged': rv.first_tagged,
        'evidence_last_updated': rv.evidence_last_updated,
        'found_in_current_run': rv.found_in_current_run,
        'reasons': rv.reasons,
        'gene': rv.gene,
        'genes': genes,
        'genotypes': dict(rv.genotypes),
        'flags': sorted(rv.flags),
        'labels': sorted(rv.labels),
        'support_vars': sorted(rv.support_vars),
        'mane_csq': mane_csq,
        'mane_hgvsps': mane_hgvsps,
        'pheno_matches': pheno_matches,
        'forced_matches': forced_matches,
        'phenotype_matches': sorted(rv.phenotype_labels),
        'var_data': var_data_dump,
        'transcript_consequences': tx_csq,
        'clinvar_stars': rv.clinvar_stars,
        'clinvar_increase': rv.clinvar_increase,
        'decision': decision,
    }


def prepare_sample(sample_id: str, participant, variants_with_decisions: list[dict]) -> dict:
    """Prepare a sample (participant) as a JSON-serializable dict."""
    meta = participant.metadata

    family_display = {mid: mid for mid in meta.members}
    family_members = {mid: m.model_dump(mode='json') for mid, m in meta.members.items()}
    phenotypes = [p.model_dump(mode='json') for p in meta.phenotypes]
    panel_details = {str(pid): p.model_dump(mode='json') for pid, p in meta.panel_details.items()}

    prepared_variants = []
    for vd in variants_with_decisions:
        rv = vd['report_variant']
        if not rv.found_in_current_run:
            continue
        prepared_variants.append(prepare_variant(rv, vd.get('decision')))

    return {
        'name': sample_id,
        'ext_id': sample_id,
        'family_id': meta.family_id,
        'solved': meta.solved,
        'phenotypes': phenotypes,
        'panel_details': panel_details,
        'family_members': family_members,
        'family_display': family_display,
        'variants': prepared_variants,
    }


def prepare_run_context(result_data: ResultData, joined: dict[str, list[dict]]) -> list[dict]:
    """Prepare all samples as JSON-serializable dicts. Returns a sorted list."""
    samples = []
    for sample_id, variants_with_decisions in joined.items():
        participant = result_data.results.get(sample_id)
        if not participant or not participant.variants:
            continue
        sample = prepare_sample(sample_id, participant, variants_with_decisions)
        if sample['variants']:
            samples.append(sample)
    samples.sort(key=lambda s: s['ext_id'])
    return samples


def prepare_metadata(result_data: ResultData) -> dict:
    """Prepare run metadata as a JSON-serializable dict."""
    meta = result_data.metadata
    return {
        'version': meta.version,
        'run_datetime': meta.run_datetime,
        'input_file': meta.input_file,
        'family_breakdown': dict(meta.family_breakdown),
        'variant_breakdown': {k: dict(v) for k, v in meta.variant_breakdown.items()},
        'samples_with_no_variants': list(meta.samples_with_no_variants),
        'panels': {str(pid): p.model_dump(mode='json') for pid, p in meta.panels.items()},
    }
