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

