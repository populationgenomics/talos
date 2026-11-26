"""
Compare a Nextflow run over nextflow/inputs/test.tsv with each case's expected.json, printing PASS/FAIL per case.

Input is the next-run TSV the run publishes to the root of its `--outdir` (talos_input_YYYYMMDD.tsv). Each row
names a cohort, its input VCF (`path`, with expected.json alongside it) and the results JSON that run produced
for it (`history`). Cohorts with no expected.json next to their VCF are skipped.

For every case, checked against the proband:
- exactly the expected set of variants is reported (an empty expectation means nothing may be reported)
- every expected category is present on the reported variant (extra categories are fine)
- every expected support variant is present (extra support variants are fine)
and that no other family member has reported variants.

Run from the repository root, exits non-zero if any case fails:
    uv run python nextflow/inputs/check_test_cases.py /path/to/outdir/talos_input_YYYYMMDD.tsv
"""

import argparse
import csv
import json
import sys
from pathlib import Path

from talos.models import ResultData
from talos.utils import read_json_from_path


def strip_chr(variant: str) -> str:
    """results coordinates drop the chr prefix, expected.json keeps it - compare without"""
    return variant.removeprefix('chr')


def check_case(expected: dict, results_file: str) -> list[str]:
    """every way this case deviates from its expected.json - empty means PASS"""
    if not results_file:
        return ['no results JSON in the history column']

    results: ResultData = read_json_from_path(results_file, return_model=ResultData)
    failures = []

    others = {
        sample: len(res.variants) for sample, res in results.results.items() if sample != 'proband' and res.variants
    }
    if others:
        failures.append(f'variants reported for non-proband samples: {others}')

    proband = results.results.get('proband')
    proband_variants = proband.variants if proband else []
    reported = {variant.var_data.coordinates.string_format: variant for variant in proband_variants}
    wanted = {strip_chr(entry['variant']): entry for entry in expected['expected']['proband']}

    if unexpected := set(reported) - set(wanted):
        failures.append(f'unexpected variants reported: {sorted(unexpected)}')
    if missing := set(wanted) - set(reported):
        failures.append(f'expected variants not reported: {sorted(missing)}')

    for key in set(wanted) & set(reported):
        entry, variant = wanted[key], reported[key]
        if missing_cats := set(entry['categories']) - set(variant.categories):
            failures.append(f'{key}: missing categories {sorted(missing_cats)}, has {sorted(variant.categories)}')
        have_support = {strip_chr(sv) for sv in variant.support_vars}
        if missing_support := {strip_chr(sv) for sv in entry['support_vars']} - have_support:
            failures.append(f'{key}: missing support variants {sorted(missing_support)}, has {sorted(have_support)}')

    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('next_run_tsv', type=Path, help='the talos_input_YYYYMMDD.tsv published to the --outdir')
    args = parser.parse_args()

    with args.next_run_tsv.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter='\t') if row.get('cohort')]

    cases = []
    for row in rows:
        expected_file = Path(row['path']).parent / 'expected.json'
        if expected_file.exists():
            cases.append((row['cohort'], json.loads(expected_file.read_text()), row.get('history') or ''))
        else:
            print(f'SKIP  {row["cohort"]}: no expected.json beside {row["path"]}')
    if not cases:
        sys.exit(f'no cohorts in {args.next_run_tsv} have an expected.json')

    width = max(len(cohort) for cohort, _, _ in cases)
    failed = 0
    for cohort, expected, results_file in cases:
        failures = check_case(expected, results_file)
        print(f'{"FAIL" if failures else "PASS"}  {cohort.ljust(width)}')
        for failure in failures:
            print(f'      - {failure}')
        failed += bool(failures)

    print(f'\n{len(cases) - failed}/{len(cases)} cases passed')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
