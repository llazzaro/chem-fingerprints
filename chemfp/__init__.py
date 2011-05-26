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

__version__ = "0.9.1"
__version_info = (0, 9, 1)
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

def read_into_memory(reader):
    # See if it has its own way to generate an in-memory search
    chemfp_in_memory = getattr(reader, "_chemfp_in_memory_", None)
    if chemfp_in_memory is not None:
        return chemfp_in_memory()

    # Nope. Use the basic forward-iteration algorithm
    from chemfp import readers
    return readers.fps_to_in_memory(reader)
    
def tanimoto_count_batch(queries, targets, threshold=0.0):
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    f = getattr(targets, "_chemfp_tanimoto_count_batch", None)
    if f is not None:
        return f(queries, threshold)
    
    from chemfp import search
    return search.generic_tanimoto_count_batch(query, targets, threshold)

def tanimoto_count(query, targets, threshold=0.0):
    return tanimoto_count_batch([queries], targets, threshold=0.0)[0]

# There is no non-batch mode because this requires indefinite memory
# and I don't see a batch search will lead to any savings.
# If I implement a batch mode, do I return a list of iterators?

def tanimoto_search(query, targets, threshold=0.0):
    """Find all targets at least 'threshold' similar to the query fingerprint

    The 'query' fingerprint string must be in binary representation, not hex.
    If 'targets' implements '_chemfp_tanimoto_search' then the search is delegated to it
      as targets._chemfp_tanimoto_search(query, threshold)
    Otherwise, 'targets' must be iterable, returning binary fingerprint strings.
    'threshold' is the minimum allowed Tanimoto score and must be between 0.0 and 1.0 .

    """
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    # Allow the targets to override the default search code
    f = getattr(targets, "_chemfp_tanimoto_search", None)
    if f is not None:
        return f(query, threshold)
    
    from chemfp import search
    return search.generic_tanimoto(query, targets, k, threshold)

def tanimoto_knearest_search_batch(queries, targets, k, threshold=0.0):
    if not (k > 0):
        raise TypeError("k must be a postive integer")
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    # Allow the targets to override the default search code
    f = getattr(targets, "_chemfp_tanimoto_knearest_search_batch", None)
    if f is not None:
        return f(queries, k, threshold)

    results = []
    from chemfp import search

    return search.generic_tanimoto_knearest_search_batch(query, targets, k, threshold)
                       
    
def tanimoto_knearest_search(query, targets, k, threshold=0.0):
    return tanimoto_knearest_search([query], targets, k, threshold)[0]
