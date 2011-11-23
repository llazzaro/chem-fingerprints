import _chemfp
from _chemfp import (hex_isvalid, hex_popcount, hex_intersect_popcount,
                     hex_tanimoto, hex_contains)

from _chemfp import (byte_popcount, byte_intersect_popcount,
                     byte_tanimoto, byte_contains)

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
