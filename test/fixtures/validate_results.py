"""Validate a Talos result JSON against an expected snapshot.

Two checks:
    1. Structural: load the result via the pydantic ResultData model. Pydantic
       enforces every type / required field; if the schema drifts, this fails.
    2. Content: deep-compare the loaded model against an expected snapshot,
       ignoring date/timestamp fields that move with the calendar.

Usage:
    python test/fixtures/validate_results.py \
        --actual   test/_run_output/needs_merge_outputs/needs_merge_full_results_*.json \
        --expected test/fixtures/expected_results/needs_merge.json

Exit code 0 = match. Non-zero = either pydantic rejected the actual file (bad
schema) or the diff after scrubbing is non-empty (regression).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

# Ensure src/ is importable when invoked directly. CI runs against an installed
# talos package, so this is a no-op there; the path injection is only for
# developer convenience when running the script from a checkout.
REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / 'src'
if SRC.exists() and str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from talos.models import ResultData  # noqa: E402

# Field names whose values move every run (timestamps, ISO dates, granular
# datetimes). The scrubber drops them anywhere they appear.
DATE_FIELD_NAMES: frozenset[str] = frozenset({
    'date_of_phenotype_match',
    'evidence_last_updated',
    'first_tagged',
    'run_datetime',
    'creation_date',
    'date',
})


def _scrub(value: Any) -> Any:
    """Recursively drop date-like keys; for the `categories` dict, replace each
    value (which is the date the category fired) with a sentinel so set
    membership is preserved but the date itself is not compared."""
    if isinstance(value, dict):
        scrubbed: dict[str, Any] = {}
        for key, sub in value.items():
            if key in DATE_FIELD_NAMES:
                continue
            if key == 'categories' and isinstance(sub, dict):
                # categories: {category_name: date_assigned} -> {category_name: '<DATE>'}
                scrubbed[key] = {k: '<DATE>' for k in sub}
                continue
            scrubbed[key] = _scrub(sub)
        return scrubbed
    if isinstance(value, list):
        return [_scrub(item) for item in value]
    return value


def _diff(path: str, actual: Any, expected: Any, out: list[str]) -> None:
    if type(actual) is not type(expected):
        out.append(f'{path}: type mismatch ({type(actual).__name__} vs {type(expected).__name__})')
        return
    if isinstance(actual, dict):
        for key in sorted(set(actual) | set(expected)):
            if key not in actual:
                out.append(f'{path}.{key}: missing in actual')
            elif key not in expected:
                out.append(f'{path}.{key}: unexpected in actual ({actual[key]!r})')
            else:
                _diff(f'{path}.{key}', actual[key], expected[key], out)
        return
    if isinstance(actual, list):
        if len(actual) != len(expected):
            out.append(f'{path}: list length differs ({len(actual)} vs {len(expected)})')
        for i, (a, e) in enumerate(zip(actual, expected)):
            _diff(f'{path}[{i}]', a, e, out)
        return
    if actual != expected:
        out.append(f'{path}: {actual!r} != {expected!r}')


def _normalise_variant_lists(payload: dict) -> None:
    """Sort each participant's variant list by coordinates so list-order doesn't matter."""
    for participant in payload.get('results', {}).values():
        variants = participant.get('variants', [])
        variants.sort(key=lambda v: (
            v.get('var_data', {}).get('coordinates', {}).get('chrom', ''),
            v.get('var_data', {}).get('coordinates', {}).get('pos', 0),
            v.get('var_data', {}).get('coordinates', {}).get('ref', ''),
            v.get('var_data', {}).get('coordinates', {}).get('alt', ''),
            v.get('sample', ''),
        ))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--actual', type=Path, required=True, help='Path to the JSON produced by the test pipeline run')
    parser.add_argument('--expected', type=Path, required=True, help='Path to the committed expected JSON snapshot')
    parser.add_argument('--write-expected', action='store_true', help='Overwrite --expected with a scrubbed copy of --actual (use to bootstrap or refresh)')
    args = parser.parse_args()

    if not args.actual.exists():
        print(f'[ERROR] actual file not found: {args.actual}', file=sys.stderr)
        return 2

    raw_actual = args.actual.read_text()
    try:
        ResultData.model_validate_json(raw_actual)
    except Exception as exc:  # pydantic.ValidationError, but keep import-light
        print(f'[ERROR] {args.actual} failed pydantic ResultData validation:\n{exc}', file=sys.stderr)
        return 3

    actual_payload = json.loads(raw_actual)
    _normalise_variant_lists(actual_payload)
    actual_scrubbed = _scrub(actual_payload)

    if args.write_expected:
        args.expected.parent.mkdir(parents=True, exist_ok=True)
        args.expected.write_text(json.dumps(actual_scrubbed, indent=2, sort_keys=True) + '\n')
        print(f'[OK] wrote scrubbed expected snapshot to {args.expected}')
        return 0

    if not args.expected.exists():
        print(f'[ERROR] expected snapshot not found: {args.expected}\n'
              f'       run with --write-expected to bootstrap.', file=sys.stderr)
        return 2

    expected_payload = json.loads(args.expected.read_text())
    differences: list[str] = []
    _diff('', actual_scrubbed, expected_payload, differences)

    if differences:
        print(f'[FAIL] {len(differences)} differences between actual and expected:', file=sys.stderr)
        for line in differences[:50]:
            print(f'  {line}', file=sys.stderr)
        if len(differences) > 50:
            print(f'  ... ({len(differences) - 50} more)', file=sys.stderr)
        return 1

    print(f'[OK] {args.actual.name} matches {args.expected.name} (after date scrub)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
