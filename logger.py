# kraken_tools/logger.py
import sys
import logging
from logging.handlers import RotatingFileHandler

def setup_logger(log_file=None, log_level=logging.INFO):
    """
    Set up a single logger (named 'kraken_analysis') that logs to console 
    and optionally to a file.
    """
    logger = logging.getLogger('kraken_analysis')
    logger.setLevel(log_level)
    # Remove any existing handlers to avoid duplication
    logger.handlers = []

    # Format for logs
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # Optional file handler
    if log_file:
        file_handler = RotatingFileHandler(log_file, maxBytes=10_485_760, backupCount=5)
        file_handler.setLevel(log_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def log_print(message, level='info'):
    """
    Print to console and also log with the 'kraken_analysis' logger.
    """
    logger = logging.getLogger('kraken_analysis')
    print(message)
    if level.lower() == 'debug':
        logger.debug(message)
    elif level.lower() == 'warning':
        logger.warning(message)
    elif level.lower() == 'error':
        logger.error(message)
    elif level.lower() == 'critical':
        logger.critical(message)
    else:
        logger.info(message)
