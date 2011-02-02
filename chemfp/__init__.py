# Library for working with cheminformatics fingerprints

# All chem-fingerprint software is distributed with the following license:

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
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

import os
import __builtin__

from . import decompressors
from .error_handlers import ChemFPError


def open_reader(source):
    from . import readers
    return readers.open_fps(source)

def open_in_memory(source):
    ext = source[-4:].lower()
    if ext == ".fps":
        from chemfp import readers
        return readers.fps_to_in_memory(readers.open_fps(source))
    elif ext == ".fpb":
        raise NotImplementedError("No support yet for .fpb files")
    # Should I open and sniff?
    raise NotImplementedError("Unknown fingerprint format extension %r" % (ext,))
    
    raise NotImplementedError

def open_mmap(source):
    raise NotImplementedError

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
