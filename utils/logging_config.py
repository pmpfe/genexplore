"""
Logging configuration for the Genetic Analysis Application.
"""

import logging
import os
import sys
import traceback
from logging.handlers import RotatingFileHandler
from datetime import datetime

from config import LOG_DIR, LOG_LEVEL, LOG_FILE_SIZE, LOG_BACKUP_COUNT, BASE_DIR

# Error log file path (in project root)
ERROR_LOG_FILE = os.path.join(BASE_DIR, 'error.log')


def setup_logging() -> logging.Logger:
    """
    Configure application logging with file and console handlers.
    Also sets up global exception handling to log all errors.
    
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
    
    # Main application log (rotating)
    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=LOG_FILE_SIZE,
        backupCount=LOG_BACKUP_COUNT
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(log_format)
    
    # Error log (all errors go here)
    error_handler = logging.FileHandler(ERROR_LOG_FILE, mode='a')
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(log_format)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(log_format)
    
    root_logger.addHandler(file_handler)
    root_logger.addHandler(error_handler)
    root_logger.addHandler(console_handler)
    
    # Set up global exception handler
    sys.excepthook = global_exception_handler
    
    return root_logger


def global_exception_handler(exc_type, exc_value, exc_traceback):
    """
    Global exception handler that logs all unhandled exceptions to error.log.
    
    Args:
        exc_type: Exception type
        exc_value: Exception value
        exc_traceback: Exception traceback
    """
    # Don't log KeyboardInterrupt
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    
    # Format the exception
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3]
    tb_lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
    tb_text = ''.join(tb_lines)
    
    # Log to error.log
    error_msg = f"{timestamp} - UNHANDLED EXCEPTION - {exc_type.__name__}: {exc_value}\n{tb_text}\n"
    
    try:
        with open(ERROR_LOG_FILE, 'a') as f:
            f.write(error_msg)
    except Exception:
        pass
    
    # Also log through logging system
    logger = logging.getLogger('unhandled')
    logger.critical(f"Unhandled exception: {exc_type.__name__}: {exc_value}", 
                   exc_info=(exc_type, exc_value, exc_traceback))
    
    # Call the default handler to print to stderr
    sys.__excepthook__(exc_type, exc_value, exc_traceback)


def log_error(message: str, exception: Exception = None) -> None:
    """
    Utility function to log an error with timestamp to error.log.
    
    Args:
        message: Error message
        exception: Optional exception object
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3]
    
    error_text = f"{timestamp} - ERROR - {message}"
    if exception:
        error_text += f"\n{traceback.format_exc()}"
    error_text += "\n"
    
    try:
        with open(ERROR_LOG_FILE, 'a') as f:
            f.write(error_text)
    except Exception:
        pass
    
    # Also log through logging system
    logger = logging.getLogger('error')
    if exception:
        logger.error(message, exc_info=True)
    else:
        logger.error(message)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance with the specified name.
    
    Args:
        name: Logger name, typically __name__ of the calling module.
        
    Returns:
        logging.Logger: Logger instance.
    """
    return logging.getLogger(name)
