import logging
import os

def configure_logger(log_file, verbose=False):
    # Get a logger instance specific to the current process
    logger = logging.getLogger()
    
    # Set the logging level
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # Check if the logger already has handlers (to avoid duplicates)
    if not logger.hasHandlers():
        # Create a file handler for logging to a specific file
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logger.level)
        
        # Set the format for log messages
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S")
        file_handler.setFormatter(formatter)
        
        # Add the handler to the logger
        logger.addHandler(file_handler)

