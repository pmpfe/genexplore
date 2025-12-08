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
import traceback
from datetime import datetime

# Ensure the application directory is in the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication, QMessageBox
from PyQt6.QtCore import Qt, qInstallMessageHandler, QtMsgType

from frontend.main_window import MainWindow
from utils.logging_config import setup_logging, get_logger, log_error, ERROR_LOG_FILE
from database.setup_database import verify_database, create_database
from config import DATABASE_PATH


def qt_message_handler(mode, context, message):
    """
    Custom Qt message handler to log Qt warnings and errors.
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S,%f')[:-3]
    
    if mode == QtMsgType.QtWarningMsg:
        level = "QT_WARNING"
    elif mode == QtMsgType.QtCriticalMsg:
        level = "QT_CRITICAL"
    elif mode == QtMsgType.QtFatalMsg:
        level = "QT_FATAL"
    else:
        return  # Don't log debug/info messages
    
    log_msg = f"{timestamp} - {level} - {message}"
    if context.file:
        log_msg += f" ({context.file}:{context.line})"
    log_msg += "\n"
    
    try:
        with open(ERROR_LOG_FILE, 'a') as f:
            f.write(log_msg)
    except Exception:
        pass


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
    # Setup logging (this also sets up global exception handler)
    setup_logging()
    logger = get_logger(__name__)
    
    # Install Qt message handler
    qInstallMessageHandler(qt_message_handler)
    
    logger.info("Starting Genetic Analysis Application")
    
    try:
        # Ensure database exists
        if not ensure_database():
            logger.error("Failed to initialize database")
            log_error("Failed to initialize database")
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
        
    except Exception as e:
        # Log any unhandled exception during startup
        log_error(f"Fatal error during application startup: {e}", e)
        logger.critical(f"Fatal error: {e}", exc_info=True)
        
        # Try to show error dialog
        try:
            app = QApplication.instance() or QApplication(sys.argv)
            QMessageBox.critical(
                None,
                "Fatal Error",
                f"A fatal error occurred:\n\n{str(e)}\n\nCheck error.log for details."
            )
        except Exception:
            pass
        
        return 1


if __name__ == '__main__':
    sys.exit(main())
