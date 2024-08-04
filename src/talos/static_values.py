"""
This is a placeholder, completely base class to prevent circular imports
"""

import logging
import sys
import zoneinfo
from datetime import datetime
from logging import FileHandler, StreamHandler

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
_GRANULAR_DATE: str | None = None
LOGGER = None


def get_granular_date():
    """
    cached getter/setter
    """
    global _GRANULAR_DATE
    if _GRANULAR_DATE is None:
        _GRANULAR_DATE = datetime.now(tz=TIMEZONE).strftime('%Y-%m-%d')
    return _GRANULAR_DATE


def get_logger(
    logger_name: str = 'Talos-logger',
    log_level: int = logging.INFO,
    file_out: str | None = None,
) -> logging.Logger:
    """
    creates a logger instance (so as not to use the root logger)
    either logs to file or stream (exclusive)

    Args:
        logger_name (str):
        log_level (int): log level to use
        file_out (str): if required, add a filehandler for logging output

    Returns:
        a logger instance, or the global logger if already defined
    """
    global LOGGER

    if LOGGER is None:
        # create a named logger
        LOGGER = logging.getLogger(logger_name)
        LOGGER.setLevel(log_level)

        # create a stream handler to write output
        handler: FileHandler | StreamHandler = (
            logging.FileHandler(file_out) if file_out else logging.StreamHandler(sys.stdout)
        )

        handler.setLevel(log_level)

        # create format string for messages
        formatter = logging.Formatter('%(asctime)s - %(name)s %(lineno)d - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)

        # set the logger to use this handler
        LOGGER.addHandler(handler)

    return LOGGER
