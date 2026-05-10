"""Build minimal clinvarbitration Hail Tables for the test fixtures.

Outputs:
    test/fixtures/processed_annotations/clinvarbitration_test.ht/
    test/fixtures/processed_annotations/clinvarbitration_test.pm5.ht/

The first HT mirrors the prod schema: rows of (clinical_significance, gold_stars,
allele_id, locus, alleles), keyed by [locus, alleles], plus globals creation_date,
blacklist. We mark exactly the test loci that the existing test data design intends
to fire ClinVar-driven categories on (cat 1 = Pathogenic/Likely Pathogenic with
stars > 0). Variants that are intentionally not in ClinVar are simply absent.

The second HT (PM5) is keyed by `newkey` ("ENST::codon"), value `clinvar_alleles`
("<allele_id>::<gold_stars>"). We seed exactly one entry so the file is structurally
valid; the test data's PKHD1 v2a was annotated as "pm5" but exercising that path
end-to-end requires the bcftools-csq output to use the same ENST/codon - too brittle
for a fixture. The category will simply not fire in CI; covered by unit tests.
"""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / 'test' / 'fixtures' / 'processed_annotations'
HT_PATH = OUTPUT_DIR / 'clinvarbitration_test.ht'
PM5_PATH = OUTPUT_DIR / 'clinvarbitration_test.pm5.ht'

# (chrom, pos, ref, alt) for variants we want to mark Pathogenic with stars > 0.
PATHOGENIC_LOCI = [
    ('chr6', 26090951, 'C', 'G'),    # HFE
    ('chr6', 52043102, 'C', 'G'),    # PKHD1 v2b
    ('chr11', 32392032, 'G', 'A'),   # WT1
    ('chr12', 89470359, 'A', 'C'),   # POC1B
    ('chrX', 71109321, 'G', 'A'),    # IL2RG
    ('chr12', 120291834, 'A', 'G'),  # RNU4-2
    ('chrM', 3243, 'A', 'G'),
    ('chrM', 8528, 'T', 'C'),
    ('chrM', 12278, 'T', 'C'),
]


def main() -> int:
    os.environ.setdefault('TALOS_CONFIG', str(REPO_ROOT / 'test' / 'input' / 'config.toml'))
    import hail as hl

    hl.init(quiet=True, default_reference='GRCh38')

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ----- main clinvarbitration HT -----
    rows = []
    for idx, (chrom, pos, ref, alt) in enumerate(PATHOGENIC_LOCI):
        rows.append(
            hl.Struct(
                locus=hl.Locus(chrom, pos, reference_genome='GRCh38'),
                alleles=[ref, alt],
                clinical_significance='Pathogenic/Likely Pathogenic',
                gold_stars=2,
                allele_id=1_000_000 + idx,
            )
        )
    # Pad to >= 100 entries to satisfy the startup_checks hard threshold. Filler
    # variants live on chr3 - a contig the test cohort never visits.
    for filler_idx in range(120):
        rows.append(
            hl.Struct(
                locus=hl.Locus('chr3', 1_000_000 + filler_idx * 100, reference_genome='GRCh38'),
                alleles=['A', 'G'],
                clinical_significance='VUS',
                gold_stars=0,
                allele_id=2_000_000 + filler_idx,
            )
        )
    ht = hl.Table.parallelize(
        rows,
        schema=hl.tstruct(
            locus=hl.tlocus('GRCh38'),
            alleles=hl.tarray(hl.tstr),
            clinical_significance=hl.tstr,
            gold_stars=hl.tint32,
            allele_id=hl.tint32,
        ),
        key=['locus', 'alleles'],
    )
    # Future date so the >2-months-old check never fires regardless of when CI runs.
    # The test config also sets clinvar_check_age=false; either alone is sufficient.
    creation_date = '2099-01-01'
    ht = ht.annotate_globals(creation_date=creation_date, blacklist=hl.empty_array(hl.tstr))

    if HT_PATH.exists():
        import shutil

        shutil.rmtree(HT_PATH)
    ht.write(str(HT_PATH), overwrite=True)
    print(f'[OK] wrote {HT_PATH} with {len(rows)} entries')

    # ----- PM5 HT -----
    pm5_rows = [hl.Struct(newkey='ENST00000000000::1', clinvar_alleles='999999::1')]
    pm5_ht = hl.Table.parallelize(
        pm5_rows,
        schema=hl.tstruct(clinvar_alleles=hl.tstr, newkey=hl.tstr),
        key=['newkey'],
    )
    pm5_ht = pm5_ht.annotate_globals(creation_date=creation_date)
    if PM5_PATH.exists():
        import shutil

        shutil.rmtree(PM5_PATH)
    pm5_ht.write(str(PM5_PATH), overwrite=True)
    print(f'[OK] wrote {PM5_PATH} with {len(pm5_rows)} entries')

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
