"""
File utility functions for the Genetic Analysis Application.
"""

import os
from typing import Optional

from utils.logging_config import get_logger

logger = get_logger(__name__)


def validate_file_exists(filepath: str) -> bool:
    """
    Check if a file exists and is readable.
    
    Args:
        filepath: Path to the file to validate.
        
    Returns:
        bool: True if file exists and is readable, False otherwise.
    """
    if not os.path.exists(filepath):
        logger.warning(f"File does not exist: {filepath}")
        return False
    
    if not os.path.isfile(filepath):
        logger.warning(f"Path is not a file: {filepath}")
        return False
    
    if not os.access(filepath, os.R_OK):
        logger.warning(f"File is not readable: {filepath}")
        return False
    
    return True


def get_file_size(filepath: str) -> Optional[int]:
    """
    Get the size of a file in bytes.
    
    Args:
        filepath: Path to the file.
        
    Returns:
        Optional[int]: File size in bytes, or None if file doesn't exist.
    """
    try:
        return os.path.getsize(filepath)
    except OSError as e:
        logger.error(f"Error getting file size for {filepath}: {e}")
        return None


def ensure_directory_exists(dirpath: str) -> bool:
    """
    Create a directory if it doesn't exist.
    
    Args:
        dirpath: Path to the directory.
        
    Returns:
        bool: True if directory exists or was created, False on error.
    """
    try:
        os.makedirs(dirpath, exist_ok=True)
        return True
    except OSError as e:
        logger.error(f"Error creating directory {dirpath}: {e}")
        return False
