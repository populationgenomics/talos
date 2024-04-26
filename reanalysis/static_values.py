"""
This is a placeholder, completely base class to prevent circular imports
"""

import logging
import sys
from datetime import datetime

from cpg_utils.config import config_retrieve

_GRANULAR_DATE: str | None = None
LOGGER = None


def get_granular_date():
    """
    cached getter/setter
    """
    global _GRANULAR_DATE
    if _GRANULAR_DATE is None:
        _GRANULAR_DATE = config_retrieve('workflow', 'fake_date', datetime.now().strftime('%Y-%m-%d'))
    return _GRANULAR_DATE


def get_logger(logger_name: str = 'AIP-logger', log_level: int = logging.INFO) -> logging.Logger:
    """
    creates a logger instance (so as not to use the root logger)

    Args:
        logger_name (str):
        log_level ():

    Returns:
        a logger instance, or the global logger if already defined
    """
    global LOGGER

    if LOGGER is None:
        # this very verbose logging is to ensure that the log level requested (INFO)
        # doesn't cause the unintentional logging of every Metamist query
        # create a named logger
        LOGGER = logging.getLogger(logger_name)
        LOGGER.setLevel(log_level)

        # create a stream handler to write output
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(log_level)

        # create format string for messages
        formatter = logging.Formatter('%(asctime)s - %(name)s %(lineno)d - %(levelname)s - %(message)s')
        stream_handler.setFormatter(formatter)

        # set the logger to use this handler
        LOGGER.addHandler(stream_handler)

    return LOGGER
