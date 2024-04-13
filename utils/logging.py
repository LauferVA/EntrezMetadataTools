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
