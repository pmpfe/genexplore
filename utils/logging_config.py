"""
Logging configuration for the Genetic Analysis Application.
"""

import logging
import os
from logging.handlers import RotatingFileHandler

from config import LOG_DIR, LOG_LEVEL, LOG_FILE_SIZE, LOG_BACKUP_COUNT


def setup_logging() -> logging.Logger:
    """
    Configure application logging with file and console handlers.
    
    Returns:
        logging.Logger: Configured root logger instance.
    """
    os.makedirs(LOG_DIR, exist_ok=True)
    
    log_file = os.path.join(LOG_DIR, 'app.log')
    
    root_logger = logging.getLogger()
    root_logger.setLevel(LOG_LEVEL)
    
    if root_logger.handlers:
        root_logger.handlers.clear()
    
    log_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=LOG_FILE_SIZE,
        backupCount=LOG_BACKUP_COUNT
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)
    
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)
    
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    return root_logger


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance with the specified name.
    
    Args:
        name: Logger name, typically __name__ of the calling module.
        
    Returns:
        logging.Logger: Logger instance.
    """
    return logging.getLogger(name)
