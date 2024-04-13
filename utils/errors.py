# __init__.py for error_handling

# Import key error handlers and decorators from the package's modules
from .base_handler import base_error_handler
from .specific_handlers import database_error_handler, network_error_handler

# Optionally, define any package-level data or initialization logic
def init_error_logging():
    # Setup global error logging configuration here
    print("Error logging is configured.")

# Call initialization functions if necessary
init_error_logging()
[vlaufer@login004 utils]$ cat logging.py 

import logging
from functools import wraps

# Configure the logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def log_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        logging.info(f"Entering {func.__name__}")
        try:
            result = func(*args, **kwargs)
            logging.info(f"Exiting {func.__name__}")
            return result
        except Exception as e:
            logging.exception(f"Error in {func.__name__} with args {args} and kwargs {kwargs}")
            raise e
    return wrapper
