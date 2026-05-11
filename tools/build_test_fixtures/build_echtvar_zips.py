"""Build all five echtvar zip fixtures by running `echtvar encode` inside Docker.

Outputs:
    test/fixtures/large_files/gnomad_4.1_region_merged_GRCh38_whole_genome
    test/fixtures/processed_annotations/alphamissense.zip
    test/fixtures/processed_annotations/mitimpact.zip
    test/fixtures/processed_annotations/mitotip.zip
    test/fixtures/processed_annotations/napogee.zip

Each zip is encoded from a synthetic single-record VCF whose INFO fields exactly
match the schema declared in the corresponding echtvar config. Test loci are
absent from every zip - echtvar's `anno` step then assigns each missing variant
the config-declared `missing_value` / `missing_string` defaults, which produces:
    - gnomad_AF_joint = 0  -> passes the `< 0.05` filter
    - am_score = -2147483648 (sentinel) -> AM stays unannotated
    - mito scores all unset
This matches the prep-workflow output shape but produces minimal coverage. CI
exercises every annotation step without depending on real population data.

Requires the talos Docker image to be available locally (any 10.x or 11.x tag).
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TOOLS_DIR = Path(__file__).resolve().parent
ECHTVAR_DIR = REPO_ROOT / 'echtvar'
HEADERS_DIR = REPO_ROOT / 'src' / 'talos' / 'vcf_headers'
LARGE_FIXTURES = REPO_ROOT / 'test' / 'fixtures' / 'large_files'
PROC_FIXTURES = REPO_ROOT / 'test' / 'fixtures' / 'processed_annotations'

DOCKER_IMAGE = 'talos:10.0.3'

# A throw-away placeholder variant on chr1 - well outside any test locus, used
# only to give echtvar at least one row to encode.
DUMMY_CHROM = 'chr1'
DUMMY_POS = 100
DUMMY_REF = 'A'
DUMMY_ALT = 'T'

VCF_HEADER_BASE = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrM,length=16569>
"""


def ftype_to_dummy(ftype: str, alias: str) -> str:
    """Pick a placeholder VCF INFO value for a field declared in an echtvar config."""
    if ftype == 'Integer':
        return '0'
    if ftype == 'Float':
        return '0.0'
    # String / Categorical
    if 'filters' in alias:
        return 'PASS'
    return '.'


def info_lines_for_config(cfg: list[dict]) -> list[str]:
    """Build VCF ##INFO header lines matching the fields declared in a config."""
    lines: list[str] = []
    for entry in cfg:
        field = entry['field']
        ftype = entry.get('ftype', 'String')
        # Categorical maps to String for VCF purposes.
        vcf_type = 'String' if ftype == 'Categorical' else ftype
        number = entry.get('number', '.')
        lines.append(f'##INFO=<ID={field},Number={number},Type={vcf_type},Description="">')
    return lines


def _load_lenient_json(text: str):
    """Tolerate trailing commas (the prod napogee config has them; echtvar is lenient)."""
    import re
    cleaned = re.sub(r',(\s*[\]}])', r'\1', text)
    return json.loads(cleaned)


def build_dummy_vcf(cfg_path: Path, out_vcf: Path) -> None:
    cfg = _load_lenient_json(cfg_path.read_text())
    info_pairs = [f'{e["field"]}={ftype_to_dummy(e.get("ftype", "String"), e["alias"])}' for e in cfg]
    info_field = ';'.join(info_pairs)
    lines: list[str] = [VCF_HEADER_BASE.rstrip('\n')]
    lines.extend(info_lines_for_config(cfg))
    lines.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
    lines.append(f'{DUMMY_CHROM}\t{DUMMY_POS}\t.\t{DUMMY_REF}\t{DUMMY_ALT}\t60\tPASS\t{info_field}')
    out_vcf.write_text('\n'.join(lines) + '\n')


def docker_run(work: Path, cmd: list[str]) -> None:
    """Run a command inside the talos Docker image with /work mounted to `work`."""
    full = [
        'docker', 'run', '--rm',
        '-v', f'{work}:/work',
        '-v', f'{REPO_ROOT}:/repo:ro',
        '-w', '/work',
        DOCKER_IMAGE,
        *cmd,
    ]
    subprocess.run(full, check=True)


def encode_zip(name: str, cfg_path: Path, out_zip: Path) -> None:
    """Build a tiny echtvar zip from a dummy VCF using `echtvar encode`."""
    work = TOOLS_DIR / '_tmp_echtvar'
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    vcf = work / f'{name}.vcf'
    build_dummy_vcf(cfg_path, vcf)

    # bgzip -> .vcf.gz (echtvar wants a compressed and indexed VCF)
    docker_run(work, ['bgzip', '-f', vcf.name])
    docker_run(work, ['tabix', '-p', 'vcf', f'{vcf.name}.gz'])

    out_zip.parent.mkdir(parents=True, exist_ok=True)
    docker_run(work, ['echtvar', 'encode', f'/work/{name}.zip', f'/repo/{cfg_path.relative_to(REPO_ROOT)}', f'{vcf.name}.gz'])
    shutil.move(work / f'{name}.zip', out_zip)
    shutil.rmtree(work)
    print(f'[OK] {out_zip} ({out_zip.stat().st_size:,} bytes)')


def main() -> int:
    targets = [
        ('gnomad', TOOLS_DIR / 'gnomad_echtvar_config.json',
         LARGE_FIXTURES / 'gnomad_4.1_region_merged_GRCh38_whole_genome'),
        ('alphamissense', ECHTVAR_DIR / 'am_config.json',
         PROC_FIXTURES / 'alphamissense.zip'),
        ('mitimpact', ECHTVAR_DIR / 'mitimpact_config.json',
         PROC_FIXTURES / 'mitimpact.zip'),
        ('mitotip', ECHTVAR_DIR / 'mitotip_config.json',
         PROC_FIXTURES / 'mitotip.zip'),
        ('napogee', ECHTVAR_DIR / 'napogee_config.json',
         PROC_FIXTURES / 'napogee.zip'),
    ]
    for name, cfg, out in targets:
        if not cfg.exists():
            print(f'[ERROR] missing {cfg}', file=sys.stderr)
            return 1
        encode_zip(name, cfg, out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
