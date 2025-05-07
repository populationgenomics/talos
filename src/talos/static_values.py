"""
This is a placeholder, completely base class to prevent circular imports
"""

import zoneinfo
from datetime import datetime

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
