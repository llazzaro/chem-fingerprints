import logging

logger = logging.getLogger("chemfp")

class ChemFPError(Exception):
    pass

class ChemFPParseError(ChemFPError):
    def __init__(self, msg, loc):
        self.msg = msg
        self.loc = loc
    def __str__(self):
        return self.msg + " " + self.loc.where()

def ignore_parse_errors(msg, loc):
    pass

def log_parse_errors(msg, loc):
    logger.error("%s %s" % (msg, loc.where()))

def strict_parse_errors(msg, loc):
    raise ChemFPParseError(msg, loc)

_parse_error_handlers = {
    "ignore": ignore_parse_errors,
    "log": log_parse_errors,
    "strict": strict_parse_errors,
    }

def get_parse_error_handler(name_or_callable="strict"):
    if name_or_callable is None:
        return strict_parse_errors
    if isinstance(name_or_callable, basestring):
        return _parse_error_handlers[name_or_callable]
    return name_or_callable
