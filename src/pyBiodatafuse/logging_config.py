# -*- coding: utf-8 -*-

"""Centralized logging configuration for pyBiodatafuse.

This module provides functions to configure logging for the entire package.
By default, logging is set to ERROR level to minimize output.

Usage:
    import pyBiodatafuse

    # Set logging level to see more information
    pyBiodatafuse.set_log_level("INFO")  # or "DEBUG", "WARNING", "ERROR", "CRITICAL"

    # Or use convenience functions
    pyBiodatafuse.enable_debug_logging()  # Sets to DEBUG level
    pyBiodatafuse.enable_info_logging()   # Sets to INFO level
    pyBiodatafuse.disable_logging()       # Sets to CRITICAL (effectively disables)
"""

import logging
import sys
from typing import IO, Optional, TextIO, Union

# Package logger name
PACKAGE_NAME = "pyBiodatafuse"

# Create package-level logger
_logger = logging.getLogger(PACKAGE_NAME)


def setup_logging(
    level: Union[int, str] = logging.ERROR,
    format_string: Optional[str] = None,
    stream: Optional[IO[str]] = None,
) -> logging.Logger:
    """Set up logging for the pyBiodatafuse package.

    :param level: Logging level (default: ERROR). Can be int or string like "INFO", "DEBUG".
    :param format_string: Custom format string for log messages.
    :param stream: Stream to output logs to (default: sys.stderr).
    :returns: The package logger.
    """
    # Convert string level to int if needed
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.ERROR)

    # Set default format if not provided
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    # Set default stream if not provided
    if stream is None:
        stream = sys.stderr

    # Configure the package logger
    _logger.setLevel(level)

    # Remove existing handlers to avoid duplicates
    _logger.handlers.clear()

    # Create and add handler
    handler = logging.StreamHandler(stream)
    handler.setLevel(level)
    formatter = logging.Formatter(format_string)
    handler.setFormatter(formatter)
    _logger.addHandler(handler)

    # Prevent propagation to root logger
    _logger.propagate = False

    return _logger


def set_log_level(level: Union[int, str]) -> None:
    """Set the logging level for pyBiodatafuse.

    :param level: Logging level. Can be int or string like "INFO", "DEBUG", "WARNING", "ERROR".

    Example:
        >>> import pyBiodatafuse
        >>> pyBiodatafuse.set_log_level("INFO")
        >>> pyBiodatafuse.set_log_level(logging.DEBUG)
    """
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.ERROR)

    _logger.setLevel(level)
    for handler in _logger.handlers:
        handler.setLevel(level)


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Get a logger for a submodule of pyBiodatafuse.

    :param name: Name of the submodule. If None, returns the package logger.
    :returns: A logger instance.

    Example:
        >>> from pyBiodatafuse.logging_config import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing data...")
    """
    if name is None:
        return _logger

    # If name already starts with package name, use as is
    if name.startswith(PACKAGE_NAME):
        return logging.getLogger(name)

    # Otherwise, prefix with package name
    return logging.getLogger(f"{PACKAGE_NAME}.{name}")


def enable_debug_logging() -> None:
    """Enable DEBUG level logging for pyBiodatafuse.

    This shows all log messages including debug information.
    """
    set_log_level(logging.DEBUG)


def enable_info_logging() -> None:
    """Enable INFO level logging for pyBiodatafuse.

    This shows informational messages, warnings, and errors.
    """
    set_log_level(logging.INFO)


def disable_logging() -> None:
    """Disable most logging for pyBiodatafuse.

    Sets level to CRITICAL, effectively silencing most messages.
    """
    set_log_level(logging.CRITICAL)


def reset_logging() -> None:
    """Reset logging to default ERROR level."""
    set_log_level(logging.ERROR)


# Initialize logging with default settings (ERROR level)
setup_logging(level=logging.ERROR)
