"internal module; not for general use"

from __future__ import absolute_import
import sys

from . import ParseError

def ignore_parse_errors(msg):
    pass

def report_parse_errors(msg):
    sys.stderr.write("ERROR: %s. Skipping.\n" % (msg,))
    
def strict_parse_errors(msg):
    raise ParseError(msg)

_parse_error_handlers = {
    "ignore": ignore_parse_errors,
    "report": report_parse_errors,
    "strict": strict_parse_errors,
    }

def get_parse_error_handler(errors):
    try:
        return _parse_error_handlers[errors]
    except KeyError:
        raise ValueError("'errors' must be one of %s" % ", ".join(sorted(_parse_error_handlers)))
        
