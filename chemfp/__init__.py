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

__version__ = "1.0b1"
__version_info = (1, 0, 0)
SOFTWARE = "chemfp/" + __version__

import os
import __builtin__
import itertools

_K = 3
_THRESHOLD = 0.7

class ChemFPError(Exception):
    pass

def read_structure_fingerprints(type, source=None, format=None, id_tag=None, aromaticity=None):
    """Read structures from `source` and return the corresponding ids and fingerprints

    This returns a FingerprintReader which can be iterated over to get
    the id and fingerprint for each read structure record. The
    fingerprint generated depends on the `type` string. Structures are
    read from `source`, which can either be the structure filename, or
    None to read from stdin.

    `type` is a format string. Examples are "OpenBabel-FP2/1", "OpenEye-Path", and
    "OpenEye-Path/1 min_bonds=0 max_bonds=5 atype=DefaultAtom btype=DefaultBond".
    Default values are used for unspecified parameters.

    If `format` is None then the structure file format and compression
    are determined by the filename's extension(s), defaulting to
    uncompressed SMILES if that is not possible. Otherwise `format` may
    be "smi" or "sdf" optionally followed by ".gz" or "bz2" to indicate
    compression.

    Note that some toolkits support additional structure formats,
    which are available to fingerprint types which use those formats.
    
    If `id_tag` is None, then the record id is based on the title
    field for the given format. If the input type is "sdf" then `id_tag`
    specifies the tag field containing the identifier. (Only the first
    line is used for multi-line values.) For example, ChEBI omits the
    title from the SD files and stores the id after the ">  <ChEBI ID>"
    line. In that case, use id_tag = "ChEBI ID".

    `aromaticity` specifies the aromaticity model, and is only appropriate for
    OEChem. It must be a string like "openeye" or "daylight".

    Here is an example of using fingerprints generated from structure file:
    
        fp_reader = read_structure_fingerprints("OpenBabel-FP4/1", "example.sdf.gz")
        print "Each fingerprint has", fps.metadata.num_bits, "bits"
        for (id, fp) in fp_reader:
           print id, repr(fp)


    :param type: The fingerprint type string. Example: 'OpenBabel-FP2/1'
    :type type: string
    :param source: The structure data source.
    :type source: A filename (as a string), a file object, or None to read from stdin.
    :param format: The file format and optional compression.
            Examples: 'smi' and 'sdf.gz'
    :type format: string, or None to autodetect based on the source
    :param id_tag: The tag containing the record id. Example: 'ChEBI ID'.
            Only valid for SD files.
    :type id_tag: string, or None to use the default title for the given format
    :param aromaticity: The aromaticity perception name (only valid for OEChem). Example: 'openeye'
    :type aromaticity: string, or None for the toolkit's default aromaticity

    :returns: an FingerprintReader

    """ # ' # emacs cruft
    from . import types
    return types.read_structure_fingerprints(type, source, format, id_tag, aromaticity)
    
# Low-memory, forward-iteration, or better
def open(source, format=None):
    """Read fingerprints from a fingerprint file

    Read fingerprints from `source`, using the given format. If
    `source` is a string then it is treated as a filename. If `source`
    is None then fingerprints are read from stdin. Otherwise, `source`
    must be a Python file object supporting 'read' and 'readline'.
    (It must also support seek() and tell() for 'FPSReader.reset()' to
    work.)

    If `format` is None then the format and compression type are
    derived from the source filename, or the name attribute of the
    source file object if it exists. If the source is None then the
    stdin is assumed to be uncompressed data in "fps" format.

    Use a `format` string to specify the format type and compression
    The supported format strings are:

       fps, fps.gz, fps.bz2  - fingerprints are in FPS format
    
    The result is an FPSReader. Here's an example of printing the
    contents of the file.
    
        reader = open("example.fps.gz")
        for id, fp in reader:
            print id, fp.encode("hex")
        
    If the `source` supports tell() and seek() then this resets
    the reader to the start of the file and finds the number of
    fingerprints in the source which are within 0.7 of the query.

        reader.reset()
        print reader.tanimoto_count_fp(query_fp)

    :param source: The fingerprint source.
    :type source: A filename string, a file object, or None.
    :param format: The file format and optional compression.
    :type format: string, or None

    :returns: an FPSReader
    """
    from . import io
    format_name, compression = io.normalize_format(source, format)

    if format_name == "fps":
        from . import readers
        return readers.open_fps(source, format_name+compression)

    if format_name == "fpb":
        raise NotImplementedError("fpb format support not implemented")

    raise TypeError("Unable to determine fingerprint format type from %r" % (source,))


def load_fingerprints(reader, metadata=None, sort=True):
    """Load all of the fingerprints into an in-memory FingerprintArena data structure

    The FingerprintArena data structure reads all of the fingerprints and
    identifers from `reader` and stores them into an in-memory data
    structure which supports fast similarity searches.

    `reader` is either an iterator over (id, fingerprint) pairs or the name
    of a file containing fingerprints. `metadata` contains the metadata metadata
    for the arena. If not specified then `reader` must contain a "metadata"
    attribute.

    The loader may reorder the fingerprints by for better search performance.
    To prevent ordering, use `sort`=False.

    :param reader: An iterator over (id, fingerprint) pairs
    :type reader: FingerprintReader or any iterator
    :param metadata: The metadata for the arena, if other than reader.metadata
    :type metadata: Metadata
    :param sort: Specify if fingerprints should be reordered for better performance
    :type sort: True or False
    :returns: FingerprintArena
    """
    if isinstance(reader, basestring):
        reader = open(reader)
    if metadata is None:
        metadata = reader.metadata

    from . import arena
    return arena.fps_to_arena(reader, metadata=metadata, sort=sort)

##### High-level search interfaces

def count_tanimoto_hits(queries, targets, threshold=_THRESHOLD, arena_size=100):
    """Count the number of targets within `threshold` of each query term

    For each query in `queries`, count the number of targets in `targets`
    which are at least `threshold` similar to the query. This function
    returns an iterator containing the (query_id, count) pairs.

    Example:

      queries = chemfp.open("queries.fps")
      targets = chemfp.load_fingerprints("targets.fps.gz")
      for (query_id, count) in chemfp.count_tanimoto_hits(queries, targets, threshold=0.9):
          print query_id, "has", count, "neighbors with at least 0.9 similarity"

    Internally, queries are processed in batches of size `arena_size`.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: `targets` must implement reset() to support multiple batches,
    and file-based targets, like an FPSReader, are slower than an Arena.

    :param queries: The query fingerprints.
    :type queries: FingerprintReader
    :param targets: The target fingerprints.
    :type targets: FingerprintArena or the slower FPSReader
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param arena_size: The number of queries to process in a batch
    :type arena_size: a positive integer, or None
    :returns:
       An iterator containing (query_id, score) pairs, one for each query
    """
    if arena_size == 1:
        for (query_id, query_fp) in queries:
            targets.reset()
            yield query_id, targets.tanimoto_count_arena(query_fp, threshold)
        return

    for query_arena in queries.iter_arenas(arena_size):
        targets.reset()
        results = targets.count_tanimoto_hits_arena(query_arena, threshold)
        for item in zip(query_arena.ids, results):
            yield item
    

def threshold_tanimoto_search(queries, targets, threshold=_THRESHOLD, arena_size=100):
    """Find all targets within `threshold` of each query term

    For each query in `queries`, find all the targets in `targets` which
    are at least `threshold` similar to the query. This function returns
    an iterator containing the (query_id, hits) pairs, where hits is a
    list of (target_id, score) pairs.

    Example:

      queries = chemfp.open("queries.fps")
      targets = chemfp.load_fingerprints("targets.fps.gz")
      for (query_id, hits) in chemfp.threshold_tanimoto_search(queries, targets, threshold=0.8):
          print query_id, "has", len(hits), "neighbors with at least 0.8 similarity"
          non_identical = [target_id for (target_id, score) in hits if score != 1.0]
          print "  The non-identical hits are:", non_identical

    Internally, queries are processed in batches of size `arena_size`.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: `targets` must implement reset() to support multiple batches,
    and file-based targets, like an FPSReader, are slower than an Arena.

    :param queries: The query fingerprints.
    :type queries: FingerprintReader
    :param targets: The target fingerprints.
    :type targets: FingerprintArena or the slower FPSReader
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param arena_size: The number of queries to process in a batch
    :type arena_size: positive integer, or None
    :returns:
      An iterator containing (query_id, hits) pairs, one for each query.
      `hits` contains a list of (target_id, score) pairs.
    """
    if arena_size == 1:
        for (query_id, query_fp) in queries:
            targets.reset()
            yield query_id, targets.threshold_tanimoto_search_fp(query_fp, threshold)
        return
    
    for query_arena in queries.iter_arenas(arena_size):
        targets.reset()
        results = targets.threshold_tanimoto_search_arena(query_arena, threshold)
        for item in zip(query_arena.ids, results):
            yield item

def knearest_tanimoto_search(queries, targets, k=_K, threshold=_THRESHOLD, arena_size=100):
    """Find the k-nearest targets within `threshold` of each query term

    For each query in `queries`, find the k-nearest of all the targets
    in `targets` which are at least `threshold` similar to the query. Ties
    are broken arbitrarily and hits with scores equal to the smallest value
    may have been omitted.
    
    This function returns an iterator containing the (query_id, hits) pairs,
    where hits is a list of (target_id, score) pairs, sorted so that the
    highest scores are first.

    Example:

      # Use the first 5 fingerprints as the queries 
      queries = next(chemfp.open("pubchem_subset.fps").iter_arenas(5))
      targets = chemfp.load_fingerprints("pubchem_subset.fps")
      
      # Find the 3 nearest hits with a similarity of at least 0.8
      for (query_id, hits) in chemfp.knearest_tanimoto_search(queries, targets, k=3, threshold=0.8):
          print query_id, "has", len(hits), "neighbors with at least 0.8 similarity"
          if hits:
              target_id, score = hits[-1]
              print "    The least similar is", target_id, "with score", score

    Internally, queries are processed in batches of size `arena_size`.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: `targets` must implement reset() to support multiple batches,
    and file-based targets, like an FPSReader, are slower than an Arena.

    :param queries: The query fingerprints.
    :type queries: FingerprintReader
    :param targets: The target fingerprints.
    :type targets: FingerprintArena or the slower FPSReader
    :param k: The maximum number of nearest neighbors to find.
    :type k: positive integer
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param arena_size: The number of queries to process in a batch
    :type arena_size: positive integer, or None
    :returns:
      An iterator containing (query_id, hits) pairs, one for each query.
      `hits` contains a list of (target_id, score) pairs.
    """
    if arena_size == 1:
        for (query_id, query_fp) in queries:
            targets.reset()
            yield query_id, targets.knearest_tanimoto_search_fp(query_fp, k, threshold)
        return

    for query_arena in queries.iter_arenas(arena_size):
        targets.reset()
        results = targets.knearest_tanimoto_search_arena(query_arena, k, threshold)
        for item in zip(query_arena.ids, results):
            yield item

def check_fp_problems(fp, metadata):
    if len(fp) != metadata.num_bytes:
        msg = ("{fp} fingerprint contains %d bytes but {metadata} has %d byte fingerprints" %
               (len(fp), metadata.num_bytes))
        return [("error", "num_bytes mismatch", msg)]
    return []

def check_metadata_problems(metadata1, metadata2):
    messages = []
    compared_num_bits = False
    if (metadata1.num_bits is not None and metadata2.num_bits is not None):
        compared_num_bits = True
        if metadata1.num_bits != metadata2.num_bits:
            msg = ("{metadata1} has %d bit fingerprints but {metadata2} has %d bit fingerprints" %
                   (metadata1.num_bits, metadata2.num_bits))
            messages.append( ("error", "num_bits mismatch", msg) )

    if (not compared_num_bits and
        metadata1.num_bytes is not None and
        metadata2.num_bytes is not None and
        metadata1.num_bytes != metadata2.num_bytes):
        
        msg = ("{metadata1} has %d byte fingerprints but {metadata2} has %d byte fingerprints" %
               (metadata1.num_bytes, metadata2.num_bytes))
        messages.append( ("error", "num_bytes mismatch", msg) )


    if (metadata1.type is not None and
        metadata2.type is not None and
        metadata1.type != metadata2.type):
        
        msg = ("{metadata1} has fingerprints of type %r but {metadata2} has fingerprints of type %r" %
               (metadata1.type, metadata2.type))
        messages.append( ("warning", "type mismatch", msg) )

    if (metadata1.aromaticity is not None and
        metadata2.aromaticity is not None and
        metadata1.aromaticity != metadata2.aromaticity):

        msg = ("{metadata1} uses aromaticity %r but {metadata2} uses aromaticity %r" %
               (metadata1.aromaticity, metadata2.aromaticity))
        messages.append( ("warning", "aromaticity mismatch", msg) )

    if (metadata1.software is not None and
        metadata2.software is not None and
        metadata1.software != metadata2.software):

        msg = ("{metadata1} comes from software %r but {metadata2} comes from software %r" %
               (metadata1.software, metadata2.software))
        messages.append( ("info", "software mismatch", msg) )

    return messages

class Metadata(object):
    def __init__(self, num_bits=None, num_bytes=None, software=None, type=None,
                 sources=None, date=None, aromaticity=None):
        if num_bytes is None:
            if num_bits is None:
                pass
            else:
                num_bytes = (num_bits + 7)//8
        elif num_bits is not None:
            if (num_bits + 7)//8 != num_bytes:
                raise TypeError("num_bits of %d is incompatible with num_bytes of %d" %
                                (num_bits, num_bytes))
            
        self.num_bits = num_bits
        self.num_bytes = num_bytes
        self.software = software
        self.type = type
        if sources is None:
            self.sources = []
        elif isinstance(sources, basestring):
            raise TypeError("sources must be a list, not a string")
        else:
            self.sources = sources
        self.date = date
        self.aromaticity = aromaticity


class FingerprintReader(object):
    """Base class for all chemfp objects holding fingerprint data

    All FingerprintReader instances have a `metadata` attribute
    containing a Metadata.
    """
    def __init__(self, metadata):
        """initialize with a Metadata instance"""
        self.metadata = metadata

    def __iter__(self):
        """iterate over the (id, fingerprint) pairs"""
        raise NotImplementedError
    
    def reset(self):
        """restart any internal iterators

        This is only relevant for fingerprint containers which have
        only one iterator. An example is the FPSReader, which uses
        stream-based file I/O to read fingerprint data.

        Calling reset() resets the iterator to its initial state.
        Iterators must allow reset() if data has not yet been read.
        Otherwise, if a reset is not possible then reset() will
        raise a TypeError.

        :returns: None
        :raises: TypeError
        """

    def iter_arenas(self, arena_size):
        """iterate through `arena_size` fingerprints at a time

        This iterates through the fingerprints `arena_size` at a time,
        yielding a FingerprintArena for each group. Working with
        arenas is often faster than processing one fingerprint at a
        time, and more memory efficient than processing all
        fingerprints at once.

        If arena_size=None then this makes an iterator containing
        a single arena containing all of the input.
        
        :param arena_size: The number of fingerprints to put into an arena.
        :type arena_size: positive integer, or None
        """
        if arena_size is None:
            yield load_fingerprints(self, self.metadata, sort=False)
            return

        if arena_size < 1:
            raise ValueError("arena_size cannot be zero")
            return
        
        it = iter(self)
        while 1:
            slice = itertools.islice(it, 0, arena_size)
            arena = load_fingerprints(slice, self.metadata, sort=False)
            if not arena:
                break
            yield arena

class FingerprintIterator(FingerprintReader):
    """A FingerprintReader for an iterator of (id, fingerprint) pairs

    This is most used used when reading files or using containers
    which only support forward iteration.
    """
    def __init__(self, metadata, id_fp_iterator):
        """initialize with a Metadata instance and the (id, fingerprint) iterator"""
        super(FingerprintIterator, self).__init__(metadata)
        self._id_fp_iterator = id_fp_iterator
        self._at_start = True

    def __iter__(self):
        """iterate over the (id, fingerprint) pairs"""
        for x in self._id_fp_iterator:
            self._at_start = False
            yield x

    def reset(self):
        """raise TypeError except if the iterator has not been used"""
        if not self._at_start:
            raise TypeError("It is not possible to reset a FingerprintIterator once it is in use")

    
class Fingerprints(FingerprintReader):
    """A FingerprintReader contining a list of (id, fingerprint) pairs

    This class add implementations for len() and []
    """
    def __init__(self, metadata, id_fp_pairs):
        """initialize with a Metadata instance and the (id, fingerprint) pair list"""
        super(Fingerprints, self).__init__(metadata)
        self._id_fp_pairs = id_fp_pairs
    def __len__(self):
        """return the number of available (id, fingerprint) pairs"""
        return len(self._id_fp_pairs)
    def __iter__(self):
        """iterate over the (id, fingerprint) pairs"""
        return iter(self._id_fp_pairs)
    
    def __repr__(self):
        return "FingerprintList(%r, %r)" % (self.metadata, self._id_fp_pairs)
    
    def __getitem__(self, i):
        """return the given (id, fingerprint) pair"""
        return self._id_fp_pairs[i]

    # Question: should I support other parts of the list API?
