"""
Typed configuration structs for each pipeline step.

Each struct is read once at step entry via its from_config() classmethod,
then passed as a parameter through the step's call graph. No config reads
inside business logic.
"""

from dataclasses import dataclass, field

import toml

_CRITICAL_CSQ_DEFAULT = [
    'frameshift',
    'splice_acceptor',
    'splice_donor',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'transcript_ablation',
]
_ADDITIONAL_CSQ_DEFAULT = ['missense', 'inframe_deletion', 'inframe_insertion']


def _load_raw(config_path: str) -> dict:
    with open(config_path) as f:
        return toml.loads(f.read())


@dataclass
class HpoFlaggingConfig:
    min_similarity: int
    phenotype_match: list[str]
    semantic_match: bool

    @classmethod
    def from_config(cls, config_path: str) -> 'HpoFlaggingConfig':
        raw = _load_raw(config_path)
        hpo = raw.get('HPOFlagging', {})
        vmi = raw.get('ValidateMOI', {})
        return cls(
            min_similarity=hpo['min_similarity'],
            phenotype_match=vmi.get('phenotype_match', []),
            semantic_match=hpo.get('semantic_match', True),
        )


@dataclass
class RunHailFilteringConfig:
    af_semi_rare: float
    csq_string: list[str]
    critical_csq: list[str]
    additional_csq: list[str]
    de_novo: dict
    am_pathogenicity: float = 0.564
    spliceai: float | None = None
    avi: float | None = None
    singletons: bool = False
    ignore_categories: list[str] = field(default_factory=list)

    @classmethod
    def from_config(cls, config_path: str) -> 'RunHailFilteringConfig':
        raw = _load_raw(config_path)
        rhf = raw.get('RunHailFiltering', {})
        vmi = raw.get('ValidateMOI', {})
        return cls(
            af_semi_rare=rhf['af_semi_rare'],
            csq_string=rhf['csq_string'],
            critical_csq=rhf.get('critical_csq', _CRITICAL_CSQ_DEFAULT),
            additional_csq=rhf.get('additional_csq', _ADDITIONAL_CSQ_DEFAULT),
            de_novo=rhf['de_novo'],
            am_pathogenicity=rhf.get('am_pathogenicity', 0.564),
            spliceai=rhf.get('spliceai'),
            avi=rhf.get('avi'),
            singletons=raw.get('singletons', False),
            ignore_categories=vmi.get('ignore_categories', []),
        )


@dataclass
class HtmlConfig:
    dataset: str | None
    sequencing_type: str | None
    long_read: bool
    long_read_defined: bool
    gnomad_population: str
    hyperlinks: bool
    hyper_template: str | None
    hyper_variant_template: str | None
    hyper_external: bool
    split_reports: bool
    default_panel: int
    forbidden_genes: set[str]
    remove_solved_cases: bool

    @classmethod
    def from_config(cls, config_path: str) -> 'HtmlConfig':
        raw = _load_raw(config_path)
        rsv = raw.get('RunHailFilteringSv', {})
        panel = raw.get('GeneratePanelData', {})
        html = raw.get('CreateTalosHTML', {})
        return cls(
            dataset=raw.get('dataset'),
            sequencing_type=raw.get('sequencing_type'),
            long_read_defined=bool('long_read' in raw),
            long_read=raw.get('long_read', False),
            default_panel=panel.get('default_panel', 137),
            forbidden_genes=set(panel.get('forbidden_genes', [])),
            remove_solved_cases=html.get('remove_solved_cases', True),
            hyperlinks='hyperlinks' in html,
            hyper_template=html.get('hyperlinks', {}).get('template'),
            hyper_variant_template=html.get('hyperlinks', {}).get('variant_template'),
            hyper_external=html.get('hyperlinks', {}).get('external'),
            split_reports=html.get('split_reports', False),
            gnomad_population=rsv.get('gnomad_population', 'gnomad_v4.1'),
        )


@dataclass
class ValidateMOIConfig:
    # GlobalFilter thresholds
    min_callset_ac_to_filter: int
    callset_max_af: float
    callset_sv_max_af: float
    gnomad_max_af: float
    gnomad_max_homozygotes: int
    gnomad_max_hemizygotes: int
    gnomad_sv_max_af: float
    # DominantFilter thresholds (stricter, applied to dominant MOI only)
    dominant_callset_max_ac: int
    dominant_callset_max_af: float
    dominant_callset_sv_max_af: float
    dominant_gnomad_max_af: float
    dominant_gnomad_max_ac: int
    dominant_gnomad_max_homozygotes: int
    dominant_gnomad_sv_max_af: float
    # ClinVarFilter thresholds (lenient, applied to ClinVar P/LP variants)
    clinvar_gnomad_max_af: float
    clinvar_callset_max_af: float
    # ClinVarDominantFilter thresholds
    clinvar_dominant_gnomad_max_af: float
    clinvar_dominant_callset_max_af: float
    # BaseMoi depth thresholds (from RunHailFiltering section)
    min_alt_depth: int = 5
    minimum_depth: int = 10
    # Mitochondrial heteroplasmy threshold
    heteroplasmy_min: float = 0.2
    # validate_moi.py runtime settings
    solved_cases: list[str] = field(default_factory=list)
    ignore_categories: list[str] = field(default_factory=list)
    singletons: bool = False

    # super logging trigger
    super_logging_path: str | None = None

    @classmethod
    def from_config(cls, config_path: str) -> 'ValidateMOIConfig':
        raw = _load_raw(config_path)
        vmi = raw.get('ValidateMOI', {})
        rhf = raw.get('RunHailFiltering', {})
        return cls(
            min_callset_ac_to_filter=vmi['min_callset_ac_to_filter'],
            callset_max_af=vmi['callset_max_af'],
            callset_sv_max_af=vmi['callset_sv_max_af'],
            gnomad_max_af=vmi['gnomad_max_af'],
            gnomad_max_homozygotes=vmi['gnomad_max_homozygotes'],
            gnomad_max_hemizygotes=vmi['gnomad_max_hemizygotes'],
            gnomad_sv_max_af=vmi['gnomad_sv_max_af'],
            dominant_callset_max_ac=vmi['dominant_callset_max_ac'],
            dominant_callset_max_af=vmi['dominant_callset_max_af'],
            dominant_callset_sv_max_af=vmi['dominant_callset_sv_max_af'],
            dominant_gnomad_max_af=vmi['dominant_gnomad_max_af'],
            dominant_gnomad_max_ac=vmi['dominant_gnomad_max_ac'],
            dominant_gnomad_max_homozygotes=vmi['dominant_gnomad_max_homozygotes'],
            dominant_gnomad_sv_max_af=vmi['dominant_gnomad_sv_max_af'],
            clinvar_gnomad_max_af=vmi['clinvar_gnomad_max_af'],
            clinvar_callset_max_af=vmi['clinvar_callset_max_af'],
            clinvar_dominant_gnomad_max_af=vmi['clinvar_dominant_gnomad_max_af'],
            clinvar_dominant_callset_max_af=vmi['clinvar_dominant_callset_max_af'],
            min_alt_depth=rhf.get('min_alt_depth', 5),
            minimum_depth=rhf.get('minimum_depth', 10),
            heteroplasmy_min=vmi.get('heteroplasmy_min', 0.2),
            solved_cases=vmi.get('solved_cases', []),
            ignore_categories=vmi.get('ignore_categories', []),
            singletons=raw.get('singletons', False),
            super_logging_path=vmi.get('super_logging_path', None),
        )
