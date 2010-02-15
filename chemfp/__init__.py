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


def open(filename):
    ext = filename[-4:].lower()
    if ext == ".fps":
        from chemfp import fps_reader
        return fps_reader.open(filename)
    elif ext == ".fpb":
        raise NotImplementedError("No support yet for .fpb files")
    # Should I open and sniff?
    raise NotImplementedError("Unknown fingerprint format")

def tanimoto_count(query, targets, threshold=0.0):
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    f = getattr(targets, "_chemfp_tanimoto_count", None)
    if f is not None:
        return f(query, threshold)
    
    from chemfp import search
    return search.generic_tanimoto_count(query, targets, k, threshold)


def tanimoto(query, targets, threshold=0.0):
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    f = getattr(targets, "_chemfp_tanimoto", None)
    if f is not None:
        return f(query, threshold)
    
    from chemfp import search
    return search.generic_tanimoto(query, targets, k, threshold)

def tanimoto_knearest(query, targets, k, threshold=0.0):
    if not (k > 0):
        raise TypeError("k must be a postive integer")
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.0, inclusive")

    # Allow the targets to override the default search code
    f = getattr(targets, "_chemfp_tanimoto_knearest", None)
    if f is not None:
        return f(query, k, threshold)

    from chemfp import search
    return search.generic_tanimoto_knearest(query, targets, k, threshold)

