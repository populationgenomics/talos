"""Subset a real panelapp_<month>.json to the genes hit by the test VCF.

Output: test/fixtures/processed_annotations/panelapp_test.json
        (the `_test` suffix matches params.fixed_month=test in the CI profile)

We retain:
  - every test gene's full record (including all panels it appears on)
  - the `versions` entries referenced by any retained gene, plus the default panel
    (137) and the forced panel (144) named in nextflow/inputs/config.toml
  - any HPO mappings for retained panels
  - the `date` field is rewritten to today so the EXPIRED_DOWNLOAD check passes

Source: any panelapp_<month>.json under nextflow/processed_annotations/.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = REPO_ROOT / 'nextflow' / 'processed_annotations'
OVERLAP_TSV = Path(__file__).resolve().parent / 'locus_gff_overlap.tsv'
OUTPUT = REPO_ROOT / 'test' / 'fixtures' / 'processed_annotations' / 'panelapp_test.json'

DEFAULT_PANEL = 137
FORCED_PANELS = {144}


def latest_source() -> Path:
    candidates = sorted(SOURCE_DIR.glob('panelapp_*.json'))
    if not candidates:
        raise FileNotFoundError(f'no panelapp_*.json found under {SOURCE_DIR}')
    return candidates[-1]


def read_target_gene_ids() -> set[str]:
    gene_ids: set[str] = set()
    lines = OVERLAP_TSV.read_text().splitlines()
    header = lines[0].split('\t')
    gid = header.index('gene_id')
    for row in lines[1:]:
        gene_ids.add(row.split('\t')[gid])
    return gene_ids


def main() -> int:
    src_path = latest_source()
    print(f'[INFO] sourcing from {src_path}')
    with open(src_path) as handle:
        data = json.load(handle)

    targets = read_target_gene_ids()
    kept_genes = {ensg: rec for ensg, rec in data['genes'].items() if ensg in targets}
    print(f'[INFO] retained {len(kept_genes)} of {len(data["genes"])} genes')

    referenced_panels: set[int] = {DEFAULT_PANEL} | FORCED_PANELS
    for rec in kept_genes.values():
        for panel_id in rec.get('panels', {}):
            referenced_panels.add(int(panel_id))

    kept_versions = [v for v in data.get('versions', []) if v['id'] in referenced_panels]
    kept_hpos = {pid: terms for pid, terms in data.get('hpos', {}).items() if int(pid) in referenced_panels}

    out = {
        'versions': kept_versions,
        'genes': kept_genes,
        'hpos': kept_hpos,
        'version': data.get('version', '2.2.0'),
        'date': datetime.now(timezone.utc).strftime('%Y-%m-%d'),
    }

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT, 'w') as handle:
        json.dump(out, handle, indent=2)
    print(f'[OK] wrote {OUTPUT}; {len(kept_versions)} panel versions, {len(kept_hpos)} hpo entries')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
