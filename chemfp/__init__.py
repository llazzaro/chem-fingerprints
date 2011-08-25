# Library for working with cheminformatics fingerprints

# All chem-fingerprint software is distributed with the following license:

# Copyright (c) 2010-2011 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

__version__ = "1.0a1"
__version_info = (1, 0, 0)
SOFTWARE = "chemfp/" + __version__

import os
import __builtin__

from .error_handlers import ChemFPError

def read_structure_fingerprints(type, source=None, format=None):
    from . import types
    return types.read_structure_fingerprints(type, source, format)
    
# Low-memory, forward-iteration, or better
def open(source, format=None, type=None):
    from . import io
    format_name, compression = io.normalize_format(source, format)

    if format_name == "fps":
        from . import readers
        return readers.open_fps(source, format_name+compression)

    if format_name == "fpb":
        raise NotImplementedError

    # Otherwise it's a structure input.
    if type is None:
        raise TypeError("'type' is required to read structure fingerprints")
    return read_structure_fingerprints(type, source, format_name+compression)

def open_fps(source):
    from . import readers
    return readers.open_fps(source)

def load_library(reader, sort=True):
    if isinstance(reader, basestring):
        reader = open(reader)
    # See if it has its own way to generate an in-memory search library
    f = getattr(reader, "_chemfp_load_library_", None)
    if f is not None:
        return f(sort=sort)

    # Nope. Use the basic forward-iteration algorithm
    from chemfp import library
    return readers.fps_to_library(reader, sort=sort)

##

# Emulate multiple dispatch

# Q is an iterator    ; fine
# Q,T is an iterator  ; can't do NxM
#   T is an iterator  ; must 
#     is an iterator 

# return a list of (id, score) pairs
# These are called "hits"

def tanimoto_count_fp(query_fp, targets, threshold=0.9):
    return targets._tanimoto_count_fp_(query_fp, targets, threshold)


# Only read the queries once
# Return a list of (query_id, hits)

def tanimoto_count(queries, targets, threshold=0.9):
    return targets._tanimoto_count_(queries, targets, threshold)

# Only read the targets once
# Return a list of (query_id, hits)

def tanimoto_count_once(queries, targets, threshold=0.9):
    return tanimoto_count_once._tanimoto_count_once_(queries, targets, threshold)

# Must read

def tanimoto_count_self(fingerprints, threshold=0.9):
    raise NotImplementedError("Someday")

##########

def threshold_tanimoto_search_fp(query_fp, targets, threshold=0.9):
    return targets._threshold_tanimoto_search_fp_(query_fp, targets, threshold)

def threshold_tanimoto_search(queries, targets, threshold=0.9):
    """Find all targets at least 'threshold' similar to the query fingerprint

    The 'query' fingerprint string must be in binary representation, not hex.
    If 'targets' implements '_chemfp_tanimoto_search' then the search is delegated to it
      as targets._chemfp_tanimoto_search(query, threshold)
    Otherwise, 'targets' must be iterable, returning binary fingerprint strings.
    'threshold' is the minimum allowed Tanimoto score and must be between 0.0 and 1.0 .

    """
    return targets._threshold_tanimoto_search_(queries, targets, threshold)

def threshold_tanimoto_search_once(queries, targets, threshold=0.9):
    return targets._threshold_tanimoto_search_once_(queries, targets, threshold)

def threshold_tanimoto_search_self(fingerprints, threshold=0.9):
    raise NotImplementedError("Someday")

##

def knearest_tanimoto_search_fp(query_fp, targets, k=3, threshold=0.9):
    return targets._knearest_tanimoto_search_fp_(query_fp, targets, k, threshold)

def knearest_tanimoto_search(queries, targets, k=3, threshold=0.9):
    return targets._knearest_tanimoto_search_(queries, targets, k, threshold)

def knearest_tanimoto_search_once(queries, targets, k=3, threshold=0.9):
    return targets._knearest_tanimoto_once_(queries, targets, k, threshold)

def knearest_tanimoto_search_self(fingerprints, k=3, threshold=0.9):
    raise NotImplementedError("Someday")
