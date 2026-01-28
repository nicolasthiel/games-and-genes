import logging
import sys
import os
from datetime import datetime

# 1. Define the Custom Formatter
class ColoredFormatter(logging.Formatter):
    """
    Formatter that adds colors to the log level names and messages.
    """
    # ANSI Escape Codes
    grey = "\x1b[38;20m"
    green = "\x1b[32;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    
    # Your desired format (Clickable)
    fmt = "%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + fmt + reset,
        logging.INFO: grey + fmt + reset,    # Kept white/grey as requested
        logging.WARNING: yellow + fmt + reset,
        logging.ERROR: red + fmt + reset,
        logging.CRITICAL: bold_red + fmt + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt="%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

def setup_logging(experiment_name: str, log_config: dict, log_dir: str = "logs"):
    # ... Create directory logic (same as before) ...
    if log_config.get('log_dir'):
        log_dir = log_config['log_dir']

    if log_config.get('log_to_file', False):
        os.makedirs(log_dir, exist_ok=True)

    logger = logging.getLogger()
    level_str = log_config.get("level", "INFO").upper()
    logger.setLevel(getattr(logging, level_str))

    # Avoid adding handlers multiple times
    if logger.handlers:
        return

    # --- HANDLER 1: CONSOLE (COLORED) ---
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(ColoredFormatter()) # <--- Use the colored class here
    logger.addHandler(console_handler)

    # --- HANDLER 2: FILE (CLEAN / NO COLORS) ---
    if log_config.get('log_to_file', False):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{experiment_name}_{timestamp}.log"
        file_path = os.path.join(log_dir, filename)
        
        file_handler = logging.FileHandler(file_path)
        
        # Use standard formatter for files so they are readable
        plain_formatter = logging.Formatter(
            "%(asctime)s - %(filename)s:%(lineno)d - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        file_handler.setFormatter(plain_formatter)
        logger.addHandler(file_handler)
        
        print(f"Logging to file: {file_path}")