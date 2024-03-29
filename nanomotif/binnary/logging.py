import logging
import os
from datetime import datetime

def set_logger_config(args):
    """
    Set the configuration for the logger.
    
    params:
        args: args.output - The output directory for the log file
    """
    log_dir = os.path.join(args.out, "log")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    
    # Format the current datetime for the log filename
    current_time = datetime.now().strftime("%Y-%m-%dT%H-%M-%S")
    
    logging.basicConfig(
        level=logging.INFO,  # Set to DEBUG, INFO, WARNING, ERROR, or CRITICAL
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=os.path.join(log_dir, f"{current_time}_log.log"),  # Omit or set to None to log to the console
        filemode='w'  # Use 'a' for append; 'w' for overwrite
    )
    
    logging.info("Logger configuration set")
    

