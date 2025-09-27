import logging
from enum import Enum
import os


if "MOLDRUG_VERBOSE" in os.environ:
    if os.environ["MOLDRUG_VERBOSE"].lower() in [1, 'true']:
        verbose = True
    elif os.environ["MOLDRUG_VERBOSE"].lower() in [0, 'false']:
        verbose = False
    else:
        raise ValueError(f"MOLDRUG_VERBOSE = {os.environ['MOLDRUG_VERBOSE']} is invalid. Choose from: 1, true, false (case insensitive).")
else:
    verbose = False


class LogLevel(Enum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


# Create a logger
logger = logging.getLogger("moldrug")

# Only configure if no handlers are set (prevents messing with user config)
if not logger.hasHandlers():
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # Set default level depending on verbose
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)


def log(msg: str, level: LogLevel = LogLevel.INFO):
    """Unified logging function for the library."""
    logger.log(level.value, msg)
