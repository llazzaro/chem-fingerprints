# Library for working with cheminformatics fingerprints

# All chem-fingerprint software is distributed with the following license:

# Copyright (c) 2010-2013 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
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

__version__ = "1.1"
__version_info = (1, 1, 0)
SOFTWARE = "chemfp/" + __version__

import os
import itertools

__all__ = ["open", "load_fingerprints", "read_structure_fingerprints",
           "Metadata", "FingerprintIterator", "Fingerprints"]
           

class ChemFPError(Exception):
    pass

class ParseError(ChemFPError):
    pass

def read_structure_fingerprints(type, source=None, format=None, id_tag=None, errors="strict"):
    """Read structures from 'source' and return the corresponding ids and fingerprints

    This returns a FingerprintReader which can be iterated over to get
    the id and fingerprint for each read structure record. The
    fingerprint generated depends on the value of 'type'. Structures
    are read from 'source', which can either be the structure
    filename, or None to read from stdin.

    'type' contains the information about how to turn a structure
    into a fingerprint. It can be a string or a metadata instance.
    String values look like "OpenBabel-FP2/1", "OpenEye-Path", and
    "OpenEye-Path/1 min_bonds=0 max_bonds=5 atype=DefaultAtom btype=DefaultBond".
    Default values are used for unspecified parameters. Use a
    Metadata instance with 'type' and 'aromaticity' values set
    in order to pass aromaticity information to OpenEye.

    If 'format' is None then the structure file format and compression
    are determined by the filename's extension(s), defaulting to
    uncompressed SMILES if that is not possible. Otherwise 'format' may
    be "smi" or "sdf" optionally followed by ".gz" or "bz2" to indicate
    compression. The OpenBabel and OpenEye toolkits also support
    additional formats.
    
    If 'id_tag' is None, then the record id is based on the title
    field for the given format. If the input format is "sdf" then 'id_tag'
    specifies the tag field containing the identifier. (Only the first
    line is used for multi-line values.) For example, ChEBI omits the
    title from the SD files and stores the id after the ">  <ChEBI ID>"
    line. In that case, use id_tag = "ChEBI ID".

    'aromaticity' specifies the aromaticity model, and is only appropriate for
    OEChem. It must be a string like "openeye" or "daylight".

    Here is an example of using fingerprints generated from structure file::
    
        fp_reader = read_structure_fingerprints("OpenBabel-FP4/1", "example.sdf.gz")
        print "Each fingerprint has", fps.metadata.num_bits, "bits"
        for (id, fp) in fp_reader:
           print id, fp.encode("hex")


    :param type: information about how to convert the input structure into a fingerprint
    :type type: string or Metadata
    :param source: The structure data source.
    :type source: A filename (as a string), a file object, or None to read from stdin
    :param format: The file format and optional compression.
            Examples: 'smi' and 'sdf.gz'
    :type format: string, or None to autodetect based on the source
    :param id_tag: The tag containing the record id. Example: 'ChEBI ID'.
            Only valid for SD files.
    :type id_tag: string, or None to use the default title for the given format
    :returns: a FingerprintReader

    """ # ' # emacs cruft
    from . import types
    if isinstance(type, basestring):
        metadata = None
    else:
        metadata = type
        if metadata.type is None:
            raise ValueError("Missing fingerprint type information in metadata")
        type = metadata.type

    structure_fingerprinter = types.parse_type(type)
    return structure_fingerprinter.read_structure_fingerprints(source, format, id_tag, errors, metadata=metadata)
    
# Low-memory, forward-iteration, or better
def open(source, format=None):
    """Read fingerprints from a fingerprint file

    Read fingerprints from 'source', using the given format. If
    'source' is a string then it is treated as a filename. If 'source'
    is None then fingerprints are read from stdin. Otherwise, 'source'
    must be a Python file object supporting 'read' and 'readline'.

    If 'format' is None then the fingerprint file format and
    compression type are derived from the source filename, or from the
    name attribute of the source file object. If the source is None
    then the stdin is assumed to be uncompressed data in "fps" format.

    The supported format strings are:

       fps, fps.gz  - fingerprints are in FPS format
    
    The result is an FPSReader. Here's an example of printing the
    contents of the file::
    
        reader = open("example.fps.gz")
        for id, fp in reader:
            print id, fp.encode("hex")
        
    :param source: The fingerprint source.
    :type source: A filename string, a file object, or None
    :param format: The file format and optional compression.
    :type format: string, or None

    :returns: an FPSReader
    """
    from . import io
    format_name, compression = io.normalize_format(source, format)

    if format_name == "fps":
        from . import fps_io
        return fps_io.open_fps(source, format_name+compression)

    if format_name == "fpb":
        raise NotImplementedError("fpb format support not implemented")
    if format is None:
        raise ValueError("Unable to determine fingerprint format type from %r" % (source,))
    else:
        raise ValueError("Unknown fingerprint format %r" % (format,))


def load_fingerprints(reader, metadata=None, reorder=True, alignment=None):
    """Load all of the fingerprints into an in-memory FingerprintArena data structure
    
    The FingerprintArena data structure reads all of the fingerprints and
    identifers from 'reader' and stores them into an in-memory data
    structure which supports fast similarity searches.
    
    If 'reader' is a string or implements "read" then the contents will be
    parsed with the 'chemfp.open' function. Otherwise it must support
    iteration returning (id, fingerprint) pairs. 'metadata' contains the
    metadata the arena. If not specified then 'reader.metadata' is used.
    
    The loader may reorder the fingerprints for better search performance.
    To prevent ordering, use reorder=False.

    The 'alignment' option specifies the alignment data alignment and
    padding size for each fingerprint. A value of 8 means that each
    fingerprint will start on a 8 byte alignment, and use storage
    space which a multiple of 8 bytes long. The default value of None
    determines the best alignment based on the fingerprint size and
    available popcount methods.

    :param reader: An iterator over (id, fingerprint) pairs
    :type reader: a string, file object, or (id, fingerprint) iterator
    :param metadata: The metadata for the arena, if other than reader.metadata
    :type metadata: Metadata
    :param reorder: Specify if fingerprints should be reordered for better performance
    :type reorder: True or False
    :returns: FingerprintArena
    :param alignment: Alignment size (both data alignment and padding) 
    """
    if isinstance(reader, basestring):
        reader = open(reader)
    elif hasattr(reader, "read"):
        reader = open(reader)
    if metadata is None:
        metadata = reader.metadata

    from . import arena
    return arena.fps_to_arena(reader, metadata=metadata, reorder=reorder,
                              alignment=alignment)

##### High-level search interfaces

def count_tanimoto_hits(queries, targets, threshold=0.7, arena_size=100):
    """Count the number of targets within 'threshold' of each query term

    For each query in 'queries', count the number of targets in 'targets'
    which are at least 'threshold' similar to the query. This function
    returns an iterator containing the (query_id, count) pairs.

    Example::

        queries = chemfp.open("queries.fps")
        targets = chemfp.load_fingerprints("targets.fps.gz")
        for (query_id, count) in chemfp.count_tanimoto_hits(queries, targets, threshold=0.9):
            print query_id, "has", count, "neighbors with at least 0.9 similarity"

    Internally, queries are processed in batches of size 'arena_size'.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: the FPSReader may be used as a target but it can only process
    one batch, and searching a FingerprintArena is faster if you have more
    than a few queries.

    :param queries: The query fingerprints.
    :type queries: any fingerprint container
    :param targets: The target fingerprints.
    :type targets: FingerprintArena or the slower FPSReader
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param arena_size: The number of queries to process in a batch
    :type arena_size: a positive integer, or None
    :returns:
       An iterator containing (query_id, score) pairs, one for each query
    """
    from . import fps_io
    if isinstance(targets, fps_io.FPSReader):
        from . import fps_search
        count_hits = fps_search.count_tanimoto_hits_arena
    else:
        from . import search
        count_hits = search.count_tanimoto_hits_arena

    ### Start the search now so compatibility errors are raised eagerly

    # Start iterating through the subarenas, and get the first of those
    subarenas = queries.iter_arenas(arena_size)
    try:
        first_query_arena = subarenas.next()
    except StopIteration:
        # There are no subarenas; return an empty iterator
        return iter([])

    # Get the first result, and hold on to it for the generator
    first_counts = count_hits(first_query_arena, targets, threshold=threshold)
    
    def count_tanimoto_hits():
        # Return results for the first arena
        for query_id, count in zip(first_query_arena.ids, first_counts):
            yield query_id, count
        # Return results for the rest of the arenas
        for query_arena in subarenas:
            counts = count_hits(query_arena, targets, threshold=threshold)
            for query_id, count in zip(query_arena.ids, counts):
                yield query_id, count
    return count_tanimoto_hits()


def threshold_tanimoto_search(queries, targets, threshold=0.7, arena_size=100):
    """Find all targets within 'threshold' of each query term

    For each query in 'queries', find all the targets in 'targets' which
    are at least 'threshold' similar to the query. This function returns
    an iterator containing the (query_id, hits) pairs. The hits are stored
    as a list of (target_id, score) pairs.

    Example::

      queries = chemfp.open("queries.fps")
      targets = chemfp.load_fingerprints("targets.fps.gz")
      for (query_id, hits) in chemfp.id_threshold_tanimoto_search(queries, targets, threshold=0.8):
          print query_id, "has", len(hits), "neighbors with at least 0.8 similarity"
          non_identical = [target_id for (target_id, score) in hits if score != 1.0]
          print "  The non-identical hits are:", non_identical

    Internally, queries are processed in batches of size 'arena_size'.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: the FPSReader may be used as a target but it can only process
    one batch, and searching a FingerprintArena is faster if you have more
    than a few queries.

    :param queries: The query fingerprints.
    :type queries: any fingerprint container
    :param targets: The target fingerprints.
    :type targets: FingerprintArena or the slower FPSReader
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param arena_size: The number of queries to process in a batch
    :type arena_size: positive integer, or None
    :returns:
      An iterator containing (query_id, hits) pairs, one for each query.
      'hits' contains a list of (target_id, score) pairs.
    """
    from . import fps_io
    if isinstance(targets, fps_io.FPSReader):
        from . import fps_search
        threshold_search = fps_search.threshold_tanimoto_search_arena
    else:
        from . import search
        threshold_search = search.threshold_tanimoto_search_arena

    ### Start the search now so compatibility errors are raised eagerly

    # Start iterating through the subarenas, and get the first of those
    subarenas = queries.iter_arenas(arena_size)
    try:
        first_query_arena = subarenas.next()
    except StopIteration:
        # There are no subarenas; return an empty iterator
        return iter([])

    # Get the first result, and hold on to it for the generator
    first_results = threshold_search(first_query_arena, targets, threshold=threshold)
    ## Here's a thought; allow a 'result_order' parameter so I can do:
    # if result_order is not None:
    #    first_results.reorder(reorder)

    def threshold_tanimoto_search():
        # Return results for the first arena
        for query_id, row in zip(first_query_arena.ids, first_results):
            yield query_id, row
        
        for query_arena in subarenas:
            results = threshold_search(query_arena, targets, threshold=threshold)
            ## I would also need to do
            #if result_order is not None:
            #    first_results.reorder(reorder)
                
            for query_id, row in zip(query_arena.ids, results):
                yield (query_id, row)
    return threshold_tanimoto_search()

def knearest_tanimoto_search(queries, targets, k=3, threshold=0.7, arena_size=100):
    """Find the 'k'-nearest targets within 'threshold' of each query term

    For each query in 'queries', find the 'k'-nearest of all the targets
    in 'targets' which are at least 'threshold' similar to the query. Ties
    are broken arbitrarily and hits with scores equal to the smallest value
    may have been omitted.
    
    This function returns an iterator containing the (query_id, hits) pairs,
    where hits is a list of (target_id, score) pairs, sorted so that the
    highest scores are first. The order of ties is arbitrary.

    Example::

      # Use the first 5 fingerprints as the queries 
      queries = next(chemfp.open("pubchem_subset.fps").iter_arenas(5))
      targets = chemfp.load_fingerprints("pubchem_subset.fps")
      
      # Find the 3 nearest hits with a similarity of at least 0.8
      for (query_id, hits) in chemfp.id_knearest_tanimoto_search(queries, targets, k=3, threshold=0.8):
          print query_id, "has", len(hits), "neighbors with at least 0.8 similarity"
          if hits:
              target_id, score = hits[-1]
              print "    The least similar is", target_id, "with score", score

    Internally, queries are processed in batches of size 'arena_size'.
    A small batch size uses less overall memory and has lower
    processing latency, while a large batch size has better overall
    performance. Use arena_size=None to process the input as a single batch.

    Note: the FPSReader may be used as a target but it can only process
    one batch, and searching a FingerprintArena is faster if you have more
    than a few queries.

    :param queries: The query fingerprints.
    :type queries: any fingerprint container
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
      'hits' contains a list of (target_id, score) pairs, sorted by score.
    """
    from . import fps_io
    if isinstance(targets, fps_io.FPSReader):
        from . import fps_search
        knearest_search = fps_search.knearest_tanimoto_search_arena
    else:
        from . import search
        knearest_search = search.knearest_tanimoto_search_arena
        
    ### Start the search now so compatibility errors are raised eagerly

    # Start iterating through the subarenas, and get the first of those
    subarenas = queries.iter_arenas(arena_size)
    try:
        first_query_arena = subarenas.next()
    except StopIteration:
        # There are no subarenas; return an empty iterator
        return iter([])

    # Get the first result, and hold on to it for the generator
    first_results = knearest_search(first_query_arena, targets, k=k, threshold=threshold)

    def knearest_tanimoto_search():
        # Return results for the first arena
        for query_id, row in zip(first_query_arena.ids, first_results):
            yield query_id, row

        # and for the subarenas
        for query_arena in subarenas:
            results = knearest_search(query_arena, targets, k=k, threshold=threshold)
            for query_id, row in zip(query_arena.ids, results):
                yield (query_id, row)
        
    return knearest_tanimoto_search()

def count_tanimoto_hits_symmetric(fingerprints, threshold=0.7):
    """Find the number of other fingerprints within `threshold` of each fingerprint
    
    For each fingerprint in the `fingerprints` arena, find the number
    of other fingerprints in the same arena which are at least
    `threshold` similar to it. The arena must have pre-computed
    popcounts. A fingerprint never matches itself.

    This function returns an iterator of (fingerprint_id, count) pairs.

    Example::

      arena = chemfp.load_fingerprints("targets.fps.gz")
      for (fp_id, count) in chemfp.count_tanimoto_hits_symmetric(arena, threshold=0.6):
          print fp_id, "has", count, "neighbors with at least 0.6 similarity"
    
    :param fingerprints: The arena containing the fingerprints.
    :type fingerprints: a FingerprintArena with precomputed popcount_indices
    :param threshold: The minimum score threshold.
    :type threshod: float between 0.0 and 1.0, inclusive
    :returns:
      An iterator of (fp_id, count) pairs, one for each fingerprint
    """
    from . import fps_io, search
    if (isinstance(fingerprints, fps_io.FPSReader) or
        not getattr(fingerprints, "popcount_indices", None)):
        raise ValueError("`fingerprints` must be a FingerprintArena with pre-computed popcount indices")

    # Start the search now so the errors are caught early
    results = search.count_tanimoto_hits_symmetric(fingerprints, threshold)
    def count_tanimoto_hits_symmetric_internal():
        for id, count in zip(fingerprints.ids, results):
            yield id, count
    return count_tanimoto_hits_symmetric_internal()

def threshold_tanimoto_search_symmetric(fingerprints, threshold=0.7):
    """Find the other fingerprints within `threshold` of each fingerprint

    For each fingerprint in the `fingerprints` arena, find the other
    fingerprints in the same arena which hare at least `threshold`
    similar to it. The arena must have pre-computed popcounts. A
    fingerprint never matches itself.

    This function returns an iterator of (fingerprint, SearchResult) pairs.
    The SearchResult hit order is arbitrary.

    Example::
    
      arena = chemfp.load_fingerprints("targets.fps.gz")
      for (fp_id, hits) in chemfp.threshold_tanimoto_search_symmetric(arena, threshold=0.75):
          print fp_id, "has", len(hits), "neighbors:"
          for (other_id, score) in hits.get_ids_and_scores():
              print "   %s  %.2f" % (other_id, score)

    :param fingerprints: The arena containing the fingerprints.
    :type fingerprints: a FingerprintArena with precomputed popcount_indices
    :param threshold: The minimum score threshold.
    :type threshod: float between 0.0 and 1.0, inclusive
    :returns: An iterator of (fp_id, SearchResult) pairs, one for each fingerprint
    """
    from . import fps_io, search
    if (isinstance(fingerprints, fps_io.FPSReader) or
        not getattr(fingerprints, "popcount_indices", None)):
        raise ValueError("`fingerprints` must be a FingerprintArena with pre-computed popcount indices")

    # Start the search now so the errors are caught early
    results = search.threshold_tanimoto_search_symmetric(fingerprints, threshold)
    def threshold_tanimoto_search_symmetric_internal():
        for id, hits in zip(fingerprints.ids, results):
            yield id, hits
    return threshold_tanimoto_search_symmetric_internal()

def knearest_tanimoto_search_symmetric(fingerprints, k=3, threshold=0.7):
    """Find the nearest `k` fingerprints within `threshold` of each fingerprint

    For each fingerprint in the `fingerprints` arena, find the nearest
    `k` fingerprints in the same arena which hare at least `threshold`
    similar to it. The arena must have pre-computed popcounts. A
    fingerprint never matches itself.

    This function returns an iterator of (fingerprint, SearchResult) pairs.
    The SearchResult hits are ordered from highest score to lowest, with
    ties broken arbitrarily.

    Example::
    
      arena = chemfp.load_fingerprints("targets.fps.gz")
      for (fp_id, hits) in chemfp.knearest_tanimoto_search_symmetric(arena, k=5, threshold=0.5):
          print fp_id, "has", len(hits), "neighbors, with scores", 
          print ", ".join("%.2f" % x for x in hits.get_scores())

    :param fingerprints: The arena containing the fingerprints.
    :type fingerprints: a FingerprintArena with precomputed popcount_indices
    :param k: The maximum number of nearest neighbors to find.
    :type k: positive integer
    :param threshold: The minimum score threshold.
    :type threshod: float between 0.0 and 1.0, inclusive
    :returns: An iterator of (fp_id, SearchResult) pairs, one for each fingerprint
    """
    from . import fps_io, search
    if (isinstance(fingerprints, fps_io.FPSReader) or
        not getattr(fingerprints, "popcount_indices", None)):
        raise ValueError("`fingerprints` must be a FingerprintArena with pre-computed popcount indices")

    # Start the search now so the errors are caught early
    results = search.knearest_tanimoto_search_symmetric(fingerprints, k, threshold)
    def knearest_tanimoto_search_symmetric_internal():
        for id, hits in zip(fingerprints.ids, results):
            yield id, hits
    return knearest_tanimoto_search_symmetric_internal()
        

def check_fp_problems(fp, metadata):
    "This interface is not documented and may change in the future"
    if len(fp) != metadata.num_bytes:
        msg = ("%%(fp)s fingerprint contains %d bytes but %%(metadata)s has %d byte fingerprints" %
               (len(fp), metadata.num_bytes))
        return [("error", "num_bytes mismatch", msg)]
    return []

def check_metadata_problems(metadata1, metadata2):
    "This interface is not documented and may change in the future"
    messages = []
    compared_num_bits = False
    if (metadata1.num_bits is not None and metadata2.num_bits is not None):
        compared_num_bits = True
        if metadata1.num_bits != metadata2.num_bits:
            msg = ("%%(metadata1)s has %d bit fingerprints but %%(metadata2)s has %d bit fingerprints" %
                   (metadata1.num_bits, metadata2.num_bits))
            messages.append( ("error", "num_bits mismatch", msg) )

    if (not compared_num_bits and
        metadata1.num_bytes is not None and
        metadata2.num_bytes is not None and
        metadata1.num_bytes != metadata2.num_bytes):
        
        msg = ("%%(metadata1)s has %d byte fingerprints but %%(metadata2) has %d byte fingerprints" %
               (metadata1.num_bytes, metadata2.num_bytes))
        messages.append( ("error", "num_bytes mismatch", msg) )


    if (metadata1.type is not None and
        metadata2.type is not None and
        metadata1.type != metadata2.type):
        
        msg = ("%%(metadata1)s has fingerprints of type %r but %%(metadata2)s has fingerprints of type %r" %
               (metadata1.type, metadata2.type))
        messages.append( ("warning", "type mismatch", msg) )

    if (metadata1.aromaticity is not None and
        metadata2.aromaticity is not None and
        metadata1.aromaticity != metadata2.aromaticity):

        msg = ("%%(metadata1)s uses aromaticity %r but %%(metadata2)s uses aromaticity %r" %
               (metadata1.aromaticity, metadata2.aromaticity))
        messages.append( ("warning", "aromaticity mismatch", msg) )

    if (metadata1.software is not None and
        metadata2.software is not None and
        metadata1.software != metadata2.software):

        msg = ("%%(metadata1)s comes from software %r but %%(metadata2)s comes from software %r" %
               (metadata1.software, metadata2.software))
        messages.append( ("info", "software mismatch", msg) )

    return messages

class Metadata(object):
    """Store information about a set of fingerprints

    The metadata attributes are:
      num_bits:
        number of bits in the fingerprint
      num_bytes:
        number of bytes in the fingerprint
      type:
        fingerprint type
      aromaticity:
        aromaticity model (only used with OEChem)
      software:
        software used to make the fingerprints
      sources:
        list of sources used to make the fingerprint
      date:
        timestamp of when the fingerprints were made

    """
    def __init__(self, num_bits=None, num_bytes=None, type=None, aromaticity=None,
                 software=None, sources=None, date=None):
        if num_bytes is None:
            if num_bits is None:
                pass
            else:
                num_bytes = (num_bits + 7)//8
        elif num_bits is None:
            num_bits = num_bytes * 8
        else:
            if (num_bits + 7)//8 != num_bytes:
                raise ValueError("num_bits of %d is incompatible with num_bytes of %d" %
                                (num_bits, num_bytes))
            
        self.num_bits = num_bits
        self.num_bytes = num_bytes
        self.type = type
        self.aromaticity = aromaticity
        self.software = software
        if sources is None:
            self.sources = []
        elif isinstance(sources, basestring):
            self.sources = [sources]
            #raise TypeError("sources must be a list, not a string")
        else:
            self.sources = sources
        self.date = date

    def __repr__(self):
        return "Metadata(num_bits=%(num_bits)r, num_bytes=%(num_bytes)r, type=%(type)r, aromaticity=%(aromaticity)r, sources=%(sources)r, software=%(software)r, date=%(date)r)" % self.__dict__

    def __str__(self):
        from cStringIO import StringIO
        from . import io
        f = StringIO()
        io.write_fps1_header(f, self)
        return f.getvalue()

class FingerprintReader(object):
    """Base class for all chemfp objects holding fingerprint records

    All FingerprintReader instances have a 'metadata' attribute
    containing a Metadata and can be iteratated over to get the (id,
    fingerprint) for each record.
    
    """
    def __init__(self, metadata):
        """Initialize with a Metadata instance"""
        self.metadata = metadata

    def __iter__(self):
        """iterate over the (id, fingerprint) pairs"""
        raise NotImplementedError
    
    def reset(self):
        """restart any internal iterators

        NOTE: method is likely to be removed in the future

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

    def iter_arenas(self, arena_size=1000):
        """iterate through 'arena_size' fingerprints at a time

        This iterates through the fingerprints 'arena_size' at a time,
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
            yield load_fingerprints(self, self.metadata, reorder=False)
            return

        if arena_size < 1:
            raise ValueError("arena_size cannot be zero")
            return
        
        it = iter(self)
        while 1:
            slice = itertools.islice(it, 0, arena_size)
            arena = load_fingerprints(slice, self.metadata, reorder=False)
            if not arena:
                break
            yield arena

    def save(self, destination):
        from . import io
        io.write_fps1_output(self, destination, self.metadata)


class FingerprintIterator(FingerprintReader):
    """A FingerprintReader for an iterator of (id, fingerprint) pairs

    This is often used as an adapter so that something which reads
    the id and fingerprint data can be used as a query source.
    
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

    This is often used as an adapter so that something which contains
    the id and fingerprint data can be used as a query source.

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
    # I almost certainly want to support slice syntax like x[:5]

def get_num_threads():
    """Return the number of OpenMP threads to use in searches

    Initially this is the value returned by omp_get_max_threads(),
    which is generally 4 unless you set the environment variable
    OMP_NUM_THREADS to some other value. 
    
    It may be any value in the range 1 to get_max_threads(), inclusive.
    """
    # I don't want the top-level chemfp module import to import a submodule.
    import _chemfp

    return _chemfp.get_num_threads()

def set_num_threads(num_threads):
    """Set the number of OpenMP threads to use in searches

    If `num_threads` is less than one then it is treated as one, and a
    value greater than get_max_threads() is treated as get_max_threads().
    """
    # I don't want the top-level chemfp module import to import a submodule.
    import _chemfp

    return _chemfp.set_num_threads(num_threads)

def get_max_threads():
    """Return the maximum number of threads available.

    If OpenMP is not available then this will return 1. Otherwise it
    returns the maximum number of threads available, as reported by
    omp_get_num_threads().
    
    """
    # I don't want the top-level chemfp module import to import a submodule.
    import _chemfp

    return _chemfp.get_max_threads()


