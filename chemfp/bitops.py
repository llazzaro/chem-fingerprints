from __future__ import absolute_import

import os
import sys
import _chemfp
from _chemfp import (hex_isvalid, hex_popcount, hex_intersect_popcount,
                     hex_tanimoto, hex_contains)

from _chemfp import (byte_popcount, byte_intersect_popcount,
                     byte_tanimoto, byte_contains,
                     byte_intersect, byte_union, byte_difference)

__all__ = ["byte_popcount", "byte_intersect_popcount",
           "byte_tanimoto", "byte_contains",
           "hex_isvalid", "hex_popcount", "hex_intersect_popcount", 
           "hex_tanimoto", "hex_contains",
           "get_methods", "get_alignments",
           "get_alignment_methods", "set_alignment_method", "select_fastest_method"]


def get_methods():
    return [_chemfp.get_method_name(i) for i in range(_chemfp.get_num_methods())]

def get_alignments():
    return [_chemfp.get_alignment_name(i) for i in range(_chemfp.get_num_alignments())]

def get_alignment_methods():
    settings = {}
    for alignment in range(_chemfp.get_num_alignments()):
        method = _chemfp.get_alignment_method(alignment)
        settings[_chemfp.get_alignment_name(alignment)] = _chemfp.get_method_name(method)
    return settings

def get_alignment_method(alignment):
    try:
        alignment_i = get_alignments().index(alignment)
    except ValueError:
        raise ValueError("Unknown alignment %r" % (alignment,))
    return _chemfp.get_method_name(_chemfp.get_alignment_method(alignment_i))

def set_alignment_method(alignment, method):
    try:
        alignment_i = get_alignments().index(alignment)
    except ValueError:
        raise ValueError("Unknown alignment %r" % (alignment,))

    try:
        method_i = get_methods().index(method)
    except ValueError:
        raise ValueError("Unknown method %r" % (method,))

    _chemfp.set_alignment_method(alignment_i, method_i)

def select_fastest_method(repeat=10000):
    if repeat > 100000:
        raise ValueError("repeat size is meaninglessly large")
    
    for alignment_i, name in enumerate(get_alignments()):
        _chemfp.select_fastest_method(alignment_i, repeat)


def get_options():
    return [_chemfp.get_option_name(i) for i in range(_chemfp.get_num_options())]

def get_option(option):
    _chemfp.get_option(option)

def set_option(option, value):
    _chemfp.set_option(option, value)

def print_report(out=sys.stdout):
    from . import SOFTWARE
    print >>out, "== Configuration report for", SOFTWARE, "=="
    print >>out, "Available methods:", " ".join(get_methods())
    print >>out, "Alignment methods:"
    for alignment in get_alignments():
        method = get_alignment_method(alignment)
        print >>out, "  %s: %s" % (alignment, method)
    print >>out, "Option settings:"
    for option in get_options():
        print >>out, "  %s: %s" % (option, _chemfp.get_option(option))

def use_environment_variables(environ=None):
    if environ is None:
        environ = os.environ

    known = set()
    for alignment in get_alignments():
        name = "CHEMFP-" + alignment.upper()
        known.add(name)
        try:
            value = environ[name]
        except KeyError:
            pass
        else:
            try:
                set_alignment_method(alignment, value)
            except ValueError, err:
                print >>sys.stderr, "WARNING: Unable to use $%s = %r: %s" % (
                    (name, value, err))

    for option in get_options():
        name = "CHEMFP-" + option.upper()
        known.add(name)
        try:
            value = environ[name]
        except KeyError:
            pass
        else:
            try:
                set_option(option, int(value))
            except ValueError, err:
                print >>sys.stderr, "WARNING: Unable to use $%s = %r: %s" % (
                    (name, value, err))

    known.add("CHEMFP-REPORT")
    report = environ.get("CHEMFP-REPORT", "0") == "1"

    if report:
        set_option("report-popcount", 1)
        set_option("report-intersect", 1)

    
    known.add("CHEMFP-PRINT-CONFIG")
    if (environ.get("CHEMFP-PRINT-CONFIG", "0") == "1" or report):
        print_report(sys.stderr)


    for k in environ:
        if not k.startswith("CHEMFP-"):
            continue
        if k not in known:
            print >>sys.stderr, "WARNING: Unknown chemfp environment variable %r" % (k,)
