#!/usr/bin/env python3
"""
Genetic Analysis Desktop Application

Main entry point for the application. Provides upload of 23andMe data,
matching against GWAS Catalog, and visualization of genetic impact scores.

Usage:
    python main.py

Requirements:
    - Python 3.9+
    - PyQt6
    - pandas
    - numpy
    
Setup:
    1. pip install -r requirements.txt
    2. python database/setup_database.py
    3. python main.py
"""

import sys
import os

# Ensure the application directory is in the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import Qt

from frontend.main_window import MainWindow
from utils.logging_config import setup_logging, get_logger
from database.setup_database import verify_database, create_database
from config import DATABASE_PATH


def ensure_database() -> bool:
    """
    Ensure the database exists, create if necessary.
    
    Returns:
        bool: True if database is available.
    """
    if verify_database(DATABASE_PATH):
        return True
    
    print("Database not found. Creating...")
    return create_database(DATABASE_PATH)


def main() -> int:
    """
    Main entry point for the application.
    
    Returns:
        int: Exit code (0 for success).
    """
    # Setup logging
    setup_logging()
    logger = get_logger(__name__)
    
    logger.info("Starting Genetic Analysis Application")
    
    # Ensure database exists
    if not ensure_database():
        logger.error("Failed to initialize database")
        print("Error: Could not initialize database. Exiting.")
        return 1
    
    # Create Qt application
    app = QApplication(sys.argv)
    app.setApplicationName("Genetic Analysis")
    app.setApplicationVersion("1.0.0")
    
    # Enable high DPI support
    app.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )
    
    # Create and show main window
    window = MainWindow()
    window.show()
    
    logger.info("Application window displayed")
    
    # Run event loop
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
