"""
Build the small-variant test cases: one trio VCF per case, each holding a single variant (or pair) of interest
inside a body of real background variation. A `single_sample` case is written instead as one VCF per trio member
(an `ss_vcf_dir` input), holding only that member's calls, to exercise the merge of single-sample VCFs.

Background comes from a real trio in a real joint-call. Example data provided here (and a good model for future
iterations) is a 1000G trio, meaning that there are no access issues with the resulting data.
The trio IDs in the input VCF are anonymised to proband/mother/father.

The sex of the real proband isn't important except for X/Y chromosomes. The created test cases use female for some test
cases, and male for some test cases - to exercise chrX/Hemizygous code paths.

For each case, every row within FETCH_WINDOW of the target(s) where at least one trio member is non-ref is pulled
from the source, then trimmed to the FLANK nearest rows either side of the target span. The synthetic variant(s)
of interest are then injected with the genotype pattern the case is testing.

INFO is dropped entirely - NormaliseVcf recomputes AC/AN/AF with `bcftools +fill-tags`. FORMAT is trimmed to the
fields run_small_filtering.py reads (GT:AD:DP:GQ:PS), so injected rows and real rows look alike. FORMAT/PL could be
estimated, but Talos ignores the field so it is dropped.

Each case directory also gets an expected.json - the machine-checkable outcome that test/test_case_expectations.py
compares against a Nextflow output directory.

Requires bcftools, bgzip, tabix on PATH, and access to SOURCE_VCF.
Run from the repository root:
    python nextflow/inputs/generate_test_cases.py [--only CASE ...]

Steps:
  - change the SOURCE_VCF path to a real VCF
  - change the sample ID translation in TRIO to match your VCF's sample IDs, always converting to proband/mother/father
  - run test data generation script
  - run main.nf on the resulting multi_case_input.tsv
  - run the evaluation script nextflow/inputs/check_test_cases.py on the 'next-run' input TSV genreated by Talos

Script output will show successful and failing test cases, along with unexpected and/or missing test cases or categories
"""

import argparse
import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent
CASES_DIR = OUTPUT_DIR / 'cases'
PEDIGREES = {
    'male': 'nextflow/inputs/pedigree_male_proband.ped',
    'female': 'nextflow/inputs/pedigree_female_proband.ped',
}

SOURCE_VCF = 'path/to/input/vcf.gzx'

# VCF name -> anonymised name
TRIO = {'PROBAND_VCF': 'proband', 'FATHER_VCF': 'father', 'MOTHER_VCF': 'mother'}
# sample column order in every output VCF, matching the anonymisation in the existing fixtures
SAMPLES = ['proband', 'mother', 'father']

# background rows kept either side of the target span
FLANK = 75
# bp either side of the target span fetched from the source before trimming to FLANK rows
FETCH_WINDOW = 250_000

# lengths from the GRCh38 reference, only the contigs used below
# if the test cases are extended to include new contigs, add them & contig length here
CONTIGS = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr6': 170805979,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr16': 90338345,
    'chrX': 156040895,
}

HEADER = """\
##fileformat=VCFv4.2
{contigs}
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the observed alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
"""

FORMAT = 'GT:AD:DP:GQ:PS'
FORMAT_KEYS = FORMAT.split(':')
# no INFO on any row - AC/AN/AF are recomputed in NormaliseVcf
INFO = '.'


@dataclass
class Call:
    """one sample's genotype call for an injected variant"""

    gt: str
    ad: list[int] | None = None
    gq: int = 60

    def format(self) -> str:
        """render as the GT:AD:DP:GQ:PL:PS FORMAT string"""
        gt = './.' if self.gt == '.' else self.gt
        # ad=[-1] means "no read-depth data" - render as missing
        missing_depth = self.ad == [-1]
        ad = '.' if missing_depth else ','.join(str(x) for x in (self.ad or [15, 15]))
        dp = '.' if missing_depth else str(sum(self.ad or [15, 15]))
        return f'{gt}:{ad}:{dp}:{self.gq}:.'


HET_CALL = Call('0/1', ad=[15, 15])
HOM_CALL = Call('1/1', ad=[0, 30])
HOMREF_CALL = Call('0/0', ad=[30, 0])
MISSING_CALL = Call('.', ad=[-1])
# AD is Number=R - one value per allele regardless of ploidy, so haploid calls still carry a REF and an ALT depth
HEMI_VAR = Call('1', ad=[0, 30])
HEMI_REF = Call('0', ad=[30, 0])

# --- shared genotype patterns for the trio ---
hom_recessive = {'proband': HOM_CALL, 'mother': HET_CALL, 'father': HET_CALL}
de_novo = {'proband': HET_CALL, 'mother': HOMREF_CALL, 'father': HOMREF_CALL}
tricky_de_novo = {'proband': HET_CALL, 'mother': MISSING_CALL, 'father': MISSING_CALL}
mat_inherited = {'proband': HET_CALL, 'mother': HET_CALL, 'father': HOMREF_CALL}
pat_inherited = {'proband': HET_CALL, 'mother': HOMREF_CALL, 'father': HET_CALL}
# hemizygous proband (male pedigree) inheriting from a het mother; checks haploid GT parsing
mat_inherited_hemi = {'proband': HEMI_VAR, 'mother': HET_CALL, 'father': HEMI_REF}
# multi-allelic site (ALT 1 and 2), each parent passing on one of the proband's ALTs - both inherited, neither de novo
split_inherited = {
    'proband': Call('1/2', ad=[0, 15, 15]),
    'mother': Call('0/2', ad=[15, 0, 15]),
    'father': Call('0/1', ad=[15, 15, 0]),
}


@dataclass
class Variant:
    """one hand-authored small variant: its locus, alleles, and each trio member's call"""

    contig: str
    pos: int
    ref: str
    alt: str
    genotypes: dict[str, Call]

    @property
    def key(self) -> str:
        return f'{self.contig}-{self.pos}-{self.ref}-{self.alt}'

    def as_row(self, samples: list[str]) -> str:
        """render as a VCF data line, restricted to the given sample columns"""
        columns = [
            self.contig,
            self.pos,
            '.',
            self.ref,
            self.alt,
            '60',
            'PASS',
            INFO,
            FORMAT,
            *(self.genotypes[sample].format() for sample in samples),
        ]
        return '\t'.join(str(column) for column in columns)


@dataclass
class Expected:
    """
    One variant the case should report for the proband.
    `categories` uses report names (models.CATEGORY_TRANSLATOR) - as it will appear in the JSON.
    """

    variant: str
    gene: str
    categories: list[str] = field(default_factory=list)
    support_vars: list[str] = field(default_factory=list)


@dataclass
class Case:
    """one test case: the variant(s) injected, the pedigree to analyse with, and what the report should contain"""

    name: str
    description: str
    pedigree: str
    contig: str
    # rows injected into the background. Empty when the variant of interest is real (already in the trio)
    targets: list[Variant] = field(default_factory=list)
    # loci the background is centred on - defaults to the injected targets' positions
    anchors: list[int] = field(default_factory=list)
    # variants the report must contain for the proband, exactly. Empty = the case must report nothing
    expected: list[Expected] = field(default_factory=list)
    # write one VCF per trio member (ss_vcf_dir input) instead of a trio VCF, exercising MergeVcfsWithBcftools
    single_sample: bool = False

    @property
    def input_path(self) -> str:
        """what the input TSV row points at: the trio VCF, or the directory of single-sample VCFs"""
        if self.single_sample:
            return f'nextflow/inputs/cases/{self.name}/single_sample_vcfs'
        return f'nextflow/inputs/cases/{self.name}/{self.name}.vcf.bgz'

    def __post_init__(self):
        if not self.anchors:
            self.anchors = [target.pos for target in self.targets]
        if not self.anchors:
            raise ValueError(f'{self.name}: a case needs targets or anchors')
        if any(target.contig != self.contig for target in self.targets):
            raise ValueError(f'{self.name}: every target must be on {self.contig}')

    @property
    def span(self) -> tuple[int, int]:
        return min(self.anchors), max(self.anchors)


CASES = [
    Case(
        name='rfx6_real_de_novo',
        description='de novo missense in RFX6 (AD)',
        pedigree='female',
        contig='chr6',
        targets=[Variant('chr6', 116927170, 'C', 'T', de_novo)],
        expected=[Expected('6-116927170-C-T', 'RFX6', ['De Novo'])],
    ),
    Case(
        name='usp48_indel_normalisation',
        description='homozygous indel needing left-alignment, no canonical consequence - checks normalisation',
        pedigree='female',
        contig='chr1',
        targets=[Variant('chr1', 21706892, 'GA', 'GAA', hom_recessive)],
        expected=[],
    ),
    Case(
        name='hfe_phenotype_match_only',
        description='homozygous ClinVar P/LP in HFE, a require_pheno_match gene with no phenotype match - filtered',
        pedigree='male',
        contig='chr6',
        targets=[Variant('chr6', 26090951, 'C', 'G', hom_recessive)],
        expected=[],
    ),
    Case(
        name='pkhd1_comp_het_pm5_clinvar',
        description='PKHD1 (AR) compound het: a de novo PM5 hit in trans with a maternal ClinVar P/LP',
        pedigree='male',
        contig='chr6',
        targets=[
            Variant('chr6', 52043699, 'T', 'A', de_novo),
            Variant('chr6', 52043102, 'C', 'G', mat_inherited),
        ],
        expected=[
            Expected('chr6-52043699-T-A', 'PKHD1', ['PM5', 'De Novo'], support_vars=['chr6-52043102-C-G']),
            Expected('chr6-52043102-C-G', 'PKHD1', ['ClinVar P/LP'], support_vars=['chr6-52043699-T-A']),
        ],
    ),
    Case(
        name='dars1_comp_het_supporting_alphamissense',
        description='DARS1 (AR) compound het: paternal high-impact reported, maternal AlphaMissense is support only',
        pedigree='female',
        contig='chr2',
        targets=[
            Variant('chr2', 135920591, 'G', 'C', pat_inherited),
            Variant('chr2', 135912503, 'G', 'A', mat_inherited),
        ],
        expected=[
            Expected(
                'chr2-135920591-G-C',
                'DARS1',
                ['AlphaMissense', 'ClinVar P/LP', 'PM5'],
                support_vars=['chr2-135912503-G-A'],
            ),
            Expected('chr2-135912503-G-A', 'DARS1', ['AlphaMissense'], support_vars=['chr2-135920591-G-C']),
        ],
    ),
    Case(
        name='daam2_alphamissense_only_filtered',
        description='de novo AlphaMissense hit in DAAM2 (AD & AR) with no phenotype match',
        pedigree='female',
        contig='chr6',
        targets=[
            Variant('chr6', 39887558, 'C', 'T', de_novo),
            Variant('chr6', 39901326, 'G', 'A', pat_inherited),
        ],
        expected=[
            Expected('chr6-39887558-C-T', 'DAAM2', ['AlphaMissense', 'De Novo'], support_vars=['chr6-39901326-G-A']),
            Expected('chr6-39901326-G-A', 'DAAM2', ['AlphaMissense'], support_vars=['chr6-39887558-C-T']),
        ],
    ),
    Case(
        name='wt1_de_novo_clinvar_high_impact',
        description='de novo ClinVar P/LP, high-impact variant in WT1 (AD)',
        pedigree='male',
        contig='chr11',
        targets=[Variant('chr11', 32392032, 'G', 'A', de_novo)],
        expected=[Expected('chr11-32392032-G-A', 'WT1', ['ClinVar P/LP', 'High Impact', 'De Novo'])],
    ),
    Case(
        name='wdr74_de_novo',
        description='de novo variant in WDR74',
        pedigree='female',
        contig='chr11',
        targets=[Variant('chr11', 62841775, 'T', 'C', de_novo)],
        expected=[Expected('chr11-62841775-T-C', 'WDR74', ['De Novo'])],
    ),
    Case(
        name='poc1b_homozygous_recessive',
        description='homozygous ClinVar P/LP, SpliceAI, high-impact variant in POC1B (AR)',
        pedigree='male',
        contig='chr12',
        targets=[Variant('chr12', 89470359, 'A', 'C', hom_recessive)],
        expected=[Expected('chr12-89470359-A-C', 'POC1B', ['ClinVar P/LP', 'High Impact'])],
    ),
    Case(
        name='il2rg_hemizygous_maternal',
        description='hemizygous ClinVar P/LP in IL2RG inherited from a het mother, haploid GT in a male proband',
        pedigree='male',
        contig='chrX',
        targets=[Variant('chrX', 71109321, 'G', 'A', mat_inherited_hemi)],
        expected=[Expected('chrX-71109321-G-A', 'IL2RG', ['ClinVar P/LP'])],
    ),
    Case(
        name='rnu4_2_de_novo_clinvar',
        description='de novo ClinVar P/LP in the non-coding gene RNU4-2 (AD)',
        pedigree='female',
        contig='chr12',
        targets=[Variant('chr12', 120291834, 'A', 'G', de_novo)],
        expected=[Expected('chr12-120291834-A-G', 'RNU4-2', ['ClinVar P/LP', 'De Novo'])],
    ),
    Case(
        name='hspa8_not_on_panel',
        description='de novo insertion in HSPA8, a gene on no panel, with parental calls missing - filtered',
        pedigree='male',
        contig='chr11',
        targets=[Variant('chr11', 123057736, 'A', 'AATC', tricky_de_novo)],
        expected=[],
    ),
    Case(
        name='ankrd11_missing_parents',
        description='high-impact deletion in ANKRD11 (AD) with parental calls missing, so not callable as de novo',
        pedigree='female',
        contig='chr16',
        targets=[Variant('chr16', 89279566, 'CCTTCGGGG', 'C', tricky_de_novo)],
        expected=[Expected('chr16-89279566-CCTTCGGGG-C', 'ANKRD11', ['High Impact'])],
    ),
    Case(
        name='wdr74_single_sample_multiallelic',
        description=(
            'single-sample VCFs: proband 1/2 at a WDR74 multi-allelic site, each parent carrying one ALT. '
            'Both alleles are inherited: T>C (ClinVar P/LP, paternal) is reported, T>G (maternal, no other category) '
            'must not become a false de novo when the merge splits it from the mother (talos#729)'
        ),
        pedigree='female',
        contig='chr11',
        targets=[Variant('chr11', 62841775, 'T', 'C,G', split_inherited)],
        expected=[Expected('chr11-62841775-T-C', 'WDR74', ['ClinVar P/LP'])],
        single_sample=True,
    ),
]


@dataclass
class BackgroundRow:
    """one real row from the source, trimmed to the FORMAT fields we keep"""

    contig: str
    pos: int
    ref: str
    alt: str
    # rendered per-sample FORMAT strings, keyed by anonymised sample name
    calls: dict[str, str]

    @property
    def key(self) -> str:
        return f'{self.contig}-{self.pos}-{self.ref}-{self.alt}'

    def as_row(self, samples: list[str]) -> str:
        columns = [self.contig, self.pos, '.', self.ref, self.alt, '.', 'PASS', INFO, FORMAT]
        columns.extend(self.calls[sample] for sample in samples)
        return '\t'.join(str(column) for column in columns)


def trim_call(format_keys: list[str], call: str) -> str:
    """keep only the FORMAT fields talos reads, in FORMAT order, missing fields rendered as '.'"""
    values = dict(zip(format_keys, call.split(':'), strict=False))
    return ':'.join(values.get(key) or '.' for key in FORMAT_KEYS)


def fetch_background(case: Case) -> list[BackgroundRow]:
    """pull every row in the fetch window where at least one trio member is non-ref"""
    start, end = case.span
    region = f'{case.contig}:{max(1, start - FETCH_WINDOW)}-{end + FETCH_WINDOW}'
    # -c 1 is applied after -s, so this keeps rows where the trio (not the whole cohort) carries an ALT
    result = subprocess.run(  # noqa: S603
        ['bcftools', 'view', '-H', '-s', ','.join(TRIO), '-c', '1', SOURCE_VCF, region],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    )

    rows = []
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        format_keys = fields[8].split(':')
        # bcftools -s emits samples in the order requested
        calls = {TRIO[src]: trim_call(format_keys, raw) for src, raw in zip(TRIO, fields[9:], strict=True)}
        rows.append(BackgroundRow(fields[0], int(fields[1]), fields[3], fields[4], calls))
    return rows


def select_background(case: Case, rows: list[BackgroundRow]) -> list[BackgroundRow]:
    """FLANK rows either side of the target span, everything inside it, minus any row an injected target replaces"""
    start, end = case.span
    injected = {target.key for target in case.targets}
    rows = sorted((row for row in rows if row.key not in injected), key=lambda row: row.pos)
    before = [row for row in rows if row.pos < start][-FLANK:]
    within = [row for row in rows if start <= row.pos <= end]
    after = [row for row in rows if row.pos > end][:FLANK]
    if len(before) < FLANK or len(after) < FLANK:
        print(f'  WARNING {case.name}: only {len(before)} rows upstream / {len(after)} downstream in the fetch window')
    return before + within + after


def render(contig: str, rows: list[BackgroundRow | Variant]) -> str:
    """build the full VCF text, one contig, rows sorted by position"""
    contigs = f'##contig=<ID={contig},length={CONTIGS[contig]},assembly=GRCh38>'
    columns = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', *SAMPLES])
    body = '\n'.join(row.as_row(SAMPLES) for row in sorted(rows, key=lambda row: (row.pos, row.ref, row.alt)))
    return f'{HEADER.format(contigs=contigs)}{columns}\n{body}\n'


def bgzip_and_tabix(text: str, stem: Path) -> None:
    """write `text` to `{stem}.vcf`, then bgzip and tabix it in place"""
    plain = stem.with_suffix('.vcf')
    plain.write_text(text)

    bgzipped = plain.with_suffix('.vcf.bgz')
    with bgzipped.open('wb') as handle:
        subprocess.run(['bgzip', '-c', str(plain)], stdout=handle, check=True)  # noqa: S603, S607
    subprocess.run(['tabix', '-p', 'vcf', str(bgzipped)], check=True)  # noqa: S603, S607
    plain.unlink()


def write_expected(case: Case, case_dir: Path, n_background: int) -> None:
    """the machine-checkable outcome for this case, consumed by test/test_case_expectations.py"""
    payload = {
        'case': case.name,
        'description': case.description,
        'pedigree': PEDIGREES[case.pedigree],
        'injected': [target.key for target in case.targets],
        'background_variants': n_background,
        'expected': {
            'proband': [
                {
                    'variant': exp.variant,
                    'gene': exp.gene,
                    'categories': exp.categories,
                    'support_vars': exp.support_vars,
                }
                for exp in case.expected
            ],
        },
    }
    (case_dir / 'expected.json').write_text(json.dumps(payload, indent=2) + '\n')


def write_input_tsv(cases: list[Case]) -> str:
    """one cohort per case - the default Nextflow test input"""
    header = ['cohort', 'path', 'type', 'pedigree', 'config']
    lines = ['\t'.join(header)]
    for case in cases:
        input_type = 'ss_vcf_dir' if case.single_sample else 'vcf'
        lines.append(
            '\t'.join([case.name, case.input_path, input_type, PEDIGREES[case.pedigree], 'nextflow/inputs/config.toml'])
        )
    (CASES_DIR / 'multi_case_input.tsv').write_text('\n'.join(lines) + '\n')
    return str(CASES_DIR / 'multi_case_input.tsv')


def split_to_single_samples(trio_vcf: Path, out_dir: Path) -> None:
    """
    replace the trio VCF with one VCF per member, holding only that member's calls - what a single-sample
    caller writes. -a trims ALTs the sample doesn't carry, so a trio multi-allelic row becomes a different
    ALT list in each sample. -c 1 drops rows the sample doesn't carry at all
    """
    out_dir.mkdir(exist_ok=True)
    for sample in SAMPLES:
        subprocess.run(  # noqa: S603
            [  # noqa: S607
                'bcftools',
                'view',
                '--no-version',
                '-s',
                sample,
                '-a',
                '-c',
                '1',
                '-Oz',
                '-W=tbi',
                '-o',
                str(out_dir / f'{sample}.vcf.bgz'),
                str(trio_vcf),
            ],
            check=True,
        )
    trio_vcf.unlink()
    trio_vcf.with_name(f'{trio_vcf.name}.tbi').unlink()


def build_case(case: Case) -> None:
    """fetch, trim, inject, write"""
    start, end = case.span
    print(f'{case.name}: fetching {case.contig}:{start}-{end} +-{FETCH_WINDOW}')
    background = select_background(case, fetch_background(case))
    if not case.targets:
        # a real variant of interest must actually be in the trio
        positions = {row.pos for row in background}
        missing = [anchor for anchor in case.anchors if anchor not in positions]
        if missing:
            raise ValueError(f'{case.name}: no trio call at {case.contig}:{missing}')

    case_dir = CASES_DIR / case.name
    case_dir.mkdir(parents=True, exist_ok=True)
    bgzip_and_tabix(render(case.contig, [*background, *case.targets]), case_dir / case.name)
    if case.single_sample:
        split_to_single_samples(case_dir / f'{case.name}.vcf.bgz', Path(case.input_path))
    write_expected(case, case_dir, len(background))
    print(f'  wrote {len(background)} background + {len(case.targets)} injected rows')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--only', nargs='*', default=None, help='case names to (re)build; default all')
    args = parser.parse_args()

    names = {case.name for case in CASES}
    if len(names) != len(CASES):
        raise ValueError('case names must be unique')
    if args.only:
        unknown = set(args.only) - names
        if unknown:
            raise SystemExit(f'unknown case(s): {sorted(unknown)}')

    for case in CASES:
        if args.only is None or case.name in args.only:
            build_case(case)

    # the TSV always covers every case, whichever were rebuilt
    input_tsv_path = write_input_tsv(CASES)
    print(f'wrote {len(CASES)} cohorts to {input_tsv_path}')


if __name__ == '__main__':
    main()
