"""
This is a placeholder, completely base class to prevent circular imports
"""

import zoneinfo
from datetime import date, datetime
from os import environ

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
_GRANULAR_DATE: str | None = None


def get_granular_date():
    """
    cached getter/setter
    """
    global _GRANULAR_DATE
    if _GRANULAR_DATE is None:
        _GRANULAR_DATE = datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d')
    return _GRANULAR_DATE


def get_evidence_date() -> str:
    """Return the logical evidence date, falling back to the real execution date."""
    evidence_date = environ.get('TALOS_EVIDENCE_DATE')
    if not evidence_date:
        return get_granular_date()
    try:
        return date.fromisoformat(evidence_date).isoformat()
    except ValueError as error:
        raise ValueError(
            f'Invalid TALOS_EVIDENCE_DATE {evidence_date!r}; expected a real date in YYYY-MM-DD format',
        ) from error
