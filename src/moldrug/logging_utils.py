import logging
from enum import Enum
from moldrug import verbose


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
