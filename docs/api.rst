.. _chemfp-api:

==========
chemfp API
==========

This chapter contains the docstrings for the public portion of the
chemfp API.

=============
chemfp module
=============

The following functions and classes are in the chemfp module.

.. py:module:: chemfp

open
====

.. py:function:: open(source, format=None)

Read fingerprints from a fingerprint file

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
:type source: A filename string, a file object, or None.
:param format: The file format and optional compression.
:type format: string, or None

:returns: an FPSReader

.. _chemfp_load_fingerprints:

load_fingerprints
=================

.. py:function:: load_fingerprints(reader, metadata=None, reorder=True)

Load all of the fingerprints into an in-memory FingerprintArena data structure

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


.. _chemfp_read_structure_fingerprints:

read_structure_fingerprints
===========================

.. py:function:: read_structure_fingerprints(type, source=None, format=None, id_tag=None, errors="strict"):

Read structures from 'source' and return the corresponding ids and fingerprints

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
:type source: A filename (as a string), a file object, or None to read from stdin.
:param format: The file format and optional compression.
        Examples: 'smi' and 'sdf.gz'
:type format: string, or None to autodetect based on the source
:param id_tag: The tag containing the record id. Example: 'ChEBI ID'.
        Only valid for SD files.
:type id_tag: string, or None to use the default title for the given format
:returns: a FingerprintReader


.. _chemfp_count_tanimoto_hits:

count_tanimoto_hits
===================

.. py:function:: count_tanimoto_hits(queries, targets, threshold=0.7, arena_size=100)

Count the number of targets within 'threshold' of each query term

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


.. _chemfp_threshold_tanimoto_search:

threshold_tanimoto_search
=========================

.. py:function:: threshold_tanimoto_search (queries, targets, threshold=0.7, arena_size=100)

Find all targets within 'threshold' of each query term

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

.. _chemfp_knearest_tanimoto_search:

knearest_tanimoto_search
========================

.. py:function:: knearest_tanimoto_search (queries, targets, k=3, threshold=0.7, arena_size=100)

Find the 'k'-nearest targets within 'threshold' of each query term

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


.. _chemfp_metadata:

Metadata
========

.. py:class:: Metadata(num_bits=None, num_bytes=None, type=None, aromaticity=None, software=None, sources=None, date=None)

Store information about a set of fingerprints

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

.. _chemfp_fingerprintreader:

FingerprintReader (base class)
==============================

.. py:class:: chemfp.FingerprintReader(metadata)

Initialize with a Metadata instance

Base class for all chemfp objects holding fingerprint records

All FingerprintReader instances have a 'metadata' attribute
containing a Metadata and can be iteratated over to get the (id,
fingerprint) for each record.

iter(arena)
-----------

.. py:method:: __iter__()

iterate over the (id, fingerprint) pairs

iter_arenas
-----------

.. py:method:: iter_arenas(arena_size=1000)

iterate through 'arena_size' fingerprints at a time

This iterates through the fingerprints 'arena_size' at a time,
yielding a FingerprintArena for each group. Working with
arenas is often faster than processing one fingerprint at a
time, and more memory efficient than processing all
fingerprints at once.

If arena_size=None then this makes an iterator containing
a single arena containing all of the input.

:param arena_size: The number of fingerprints to put into an arena.
:type arena_size: positive integer, or None


===================
chemfp.arena module
===================

FingerprintArena instances are returned as part of the public API but
should not be constructed directly.

.. _chemfp_arena_fingerprintarena:

.. py:module:: chemfp.arena

FingerprintArena
================

Implements the FingerprintReader interface.

.. py:class:: FingerprintArena(... do not call directly ...)

Stores fingerprints in a contiguous block of memory

The public attributes are:
   metadata
       `Metadata` about the fingerprints
   ids
       list of identifiers, ordered by position

arena.ids
---------

A list of the fingerprint identifiers, in the same order as the
fingerprints.


len(arena)
----------

.. py:method:: __len__()

Number of fingerprint records in the FingerprintArena

arena[i]
--------

.. py:method:: __getitem__(i)

Return the (id, fingerprint) at position i


copy
----

.. py:method:: copy(indices=None, reorder=None)

None

get_by_id
---------

.. py:method:: get_by_id(id)

Given the record identifier, return the (id, fingerprint) tuple or None if not present

get_fingerprint_by_id
---------------------

.. py:method:: get_fingerprint_by_id(id)

Given the record identifier, return its fingerprint or None if not present

get_index_by_id
---------------

.. py:method:: get_index_by_id(id)

Given the record identifier, return the record index or None if not present


iter(arena)
-----------

.. py:method:: __iter__()

Iterate over the (id, fingerprint) contents of the arena


iter_arenas
-----------

.. py:method:: iter_arenas(arena_size=1000)

iterate through `arena_size` fingerprints at a time

This iterates through the fingerprints `arena_size` at a time,
yielding a FingerprintArena for each group. Working with
arenas is often faster than processing one fingerprint at a
time, and more memory efficient than processing all
fingerprints at once.

If arena_size=None then this makes an iterator containing
a single arena containing all of the input.

:param arena_size: The number of fingerprints to put into an arena.
:type arena_size: positive integer, or None


save
----

.. py:method:: save(destination)

Save the arena contents to the given filename or file object


count_tanimoto_hits_fp
----------------------

(Deprecated. Please use chemfp.search.count_tanimoto_hits_fp.)

.. py:method:: count_tanimoto_hits_fp(query_fp, threshold=0.7)

Count the fingerprints which are similar enough to the query fingerprint

XXX

Return the number of fingerprints in this arena which are
at least `threshold` similar to the query fingerprint `query_fp`.

:param query_fp: query fingerprint
:type query_fp: byte string
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: integer count


count_tanimoto_hits_arena
-------------------------

(Deprecated. Please use chemfp.search.count_tanimoto_hits_arena
or chemfp.search.count_tanimoto_hits_symmetric.)

.. py:method:: count_tanimoto_hits_arena(query_arena, threshold=0.7)

Count the fingerprints which are similar enough to each query fingerprint

XXX

Returns an iterator containing the (query_id, count) for each
fingerprint in `queries`, where `query_id` is the query
fingerprint id and `count` is the number of fingerprints found
which are at least `threshold` similar to the query.

The order of results is the same as the order of the
queries. For efficiency reasons, `arena_size` queries are
processed at a time.

:param queries: query fingerprints
:type query_fp: FingerprintArena or FPSReader (must implement iter_arenas())
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:param arena_size: number of queries to process at a time (default: 100)
:type arena_size: positive integer
:returns: list of (query_id, integer count) pairs, one for each query

threshold_tanimoto_search_fp
----------------------------

(Deprecated. Please use chemfp.search.threshold_tanimoto_search_fp.)

.. py:method:: threshold_tanimoto_search_fp(query_fp, threshold=0.7)

Find the fingerprints which are similar enough to the query fingerprint

XXX

Find all of the fingerprints in this arena which are at least
`threshold` similar to the query fingerprint `query_fp`.
The hits are returned as a list containing (id, score) tuples
in arbitrary order.

:param query_fp: query fingerprint
:type query_fp: byte string
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: list of (int, score) tuples


threshold_tanimoto_search_arena
-------------------------------

(Deprecated. Please use chemfp.search.threshold_tanimoto_search_arena
or chemfp.search.threshold_tanimoto_search_symmetric.)

.. py:method:: threshold_tanimoto_search_arena(query_arena, threshold=0.7)

Find the fingerprints which are similar to each of the query fingerprints

XXX

For each fingerprint in the `query_arena`, find all of the
fingerprints in this arena which are at least `threshold`
similar. The hits are returned as a `SearchResults` instance.

:param query_arena: query arena
:type query_arena: FingerprintArena
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: SearchResults

knearest_tanimoto_search_fp
----------------------------

(Deprecated. Please use chemfp.search.knearest_tanimoto_search_fp.)

.. py:method:: knearest_tanimoto_search_fp(query_fp, k=3, threshold=0.7)

Find the k-nearest fingerprints which are similar to the query fingerprint

XXX

Find the `k` fingerprints in this arena which are most similar
to the query fingerprint `query_fp` and which are at least `threshold`
similar to the query. The hits are returned as a list of
(id, score) tuples sorted with the highest similarity first.
Ties are broken arbitrarily.

:param query_fp: query fingerpring
:type query_fp: byte string
:param k: number of nearest neighbors to find (default: 3)
:type k: positive integer
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: SearchResults


knearest_tanimoto_search_arena
-------------------------------

(Deprecated. Please use chemfp.search.knearest_tanimoto_search_arena
or chemfp.search.knearest_tanimoto_search_symmetric.)

.. py:method:: knearest_tanimoto_search_arena(query_arena, k=3, threshold=0.7)

Find the k-nearest fingerprint which are similar to each of the query fingerprints

XXX

For each fingerprint in the `query_arena`, find the `k`
fingerprints in this arena which are most similar and which
are at least `threshold` similar to the query fingerprint.
The hits are returned as a SearchResult where the hits are
sorted with the highest similarity first. Ties are broken
arbitrarily.

:param query_arena: query arena
:type query_arena: FingerprintArena
:param k: number of nearest neighbors to find (default: 3)
:type k: positive integer
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: SearchResult


====================
chemfp.search module
====================

The following functions and classes are in the chemfp.search module.

.. py:module:: chemfp.search

Module functions
================

The `*_fp` functions search a query fingerprint against a target
arena. The `*_arena` functions search a query arena against a target
arena. The `*_symmetric` functions use the same arena as query and
target, and exclude matching a fingerprint against itself.

count_tanimoto_hits_fp
----------------------

.. py:method:: count_tanimoto_hits_fp(query_fp, target_arena, threshold=0.7)

Count the number of hits in `target_arena` at least `threshold` similar to the `query_fp`

Example::

    query_id, query_fp = chemfp.load_fingerprints("queries.fps")[0]
    targets = chemfp.load_fingerprints("targets.fps")
    print chemfp.search.count_tanimoto_hits_fp(query_fp, targets, threshold=0.1)
    

:param query_fp: the query fingerprint
:type query_fp: a byte string
:param target_arena: the target arena
:type target_fp: a FingerprintArena
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: an integer count

count_tanimoto_hits_arena
-------------------------

.. py:method:: count_tanimoto_hits_arena(query_arena, target_arena, threshold=0.7)

For each fingerprint in `query_arena`, count the number of hits in `target_arena` at least `threshold` similar to it

Example::

    queries = chemfp.load_fingerprints("queries.fps")
    targets = chemfp.load_fingerprints("targets.fps")
    counts = chemfp.search.count_tanimoto_hits_arena(queries, targets, threshold=0.1)
    print counts[:10]

The result is implementation specific. You'll always be able to
get its length and do an index lookup to get an integer
count. Currently it's a ctype array of longs, but it could be an
array.array or Python list in the future.

:param query_arena: The query fingerprints.
:type query_arena: a FingerprintArena
:param target_arena: The target fingerprints.
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: an array of counts

count_tanimoto_hits_symmetric
-----------------------------

.. py:method:: count_tanimoto_hits_symmetric(arena, threshold=0.7, batch_size=100)

For each fingerprint in the `arena`, count the number of other fingerprints at least `threshold` similar to it

A fingerprint never matches itself.

The computation can take a long time. Python won't check check for
a ^C until the function finishes. This can be irritating. Instead,
process only `batch_size` rows at a time before checking for a ^C.

Example::

    arena = chemfp.load_fingerprints("targets.fps")
    counts = chemfp.search.count_tanimoto_hits_symmetric(arena, threshold=0.2)
    print counts[:10]

The result object is implementation specific. You'll always be able to
get its length and do an index lookup to get an integer
count. Currently it's a ctype array of longs, but it could be an
array.array or Python list in the future.

:param arena: the set of fingerprints
:type arena: a FingerprintArena
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:param batch_size: the number of rows to process before checking for a ^C
:type batch_size: integer
:returns: an array of counts


threshold_tanimoto_search_fp
----------------------------

.. py:method:: threshold_tanimoto_search_fp(query_fp, target_arena, threshold=0.7)

Search for fingerprint hits in `target_arena` which are at least `threshold` similar to `query_fp`

The hits in the returned `SearchResult` are in arbitrary order.

Example::

    query_id, query_fp = chemfp.load_fingerprints("queries.fps")[0]
    targets = chemfp.load_fingerprints("targets.fps")
    print list(chemfp.search.threshold_tanimoto_search_fp(query_fp, targets, threshold=0.15))

:param query_fp: the query fingerprint
:type query_fp: a byte string
:param target_arena: the target arena
:type target_fp: a FingerprintArena
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: a SearchResult

threshold_tanimoto_search_arena
-------------------------------

.. py:method:: threshold_tanimoto_search_arena(query_arena, target_arena, threshold=0.7)

Search for the hits in the `target_arena` at least `threshold` similar to the fingerprints in `query_arena`

The hits in the returned `SearchResults` are in arbitrary order.

Example::

    queries = chemfp.load_fingerprints("queries.fps")
    targets = chemfp.load_fingerprints("targets.fps")
    results = chemfp.search.threshold_tanimoto_search_arena(queries, targets, threshold=0.5)
    for query_id, query_hits in zip(queries.ids, results):
        if len(query_hits) > 0:
            print query_id, "->", ", ".join(query_hits.get_ids())

:param query_arena: The query fingerprints.
:type query_arena: a FingerprintArena
:param target_arena: The target fingerprints.
:type target_arena: a FingerprintArena
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: a SearchResults instance

threshold_tanimoto_search_symmetric
-----------------------------------

.. py:method:: threshold_tanimoto_search_symmetric(arena, threshold=0.7, include_lower_triangle=True, batch_size=100)

Search for the hits in the `arena` at least `threshold` similar to the fingerprints in the arena

When `include_lower_triangle` is True, compute the upper-triangle
similarities, then copy the results to get the full set of
results. When `include_lower_triangle` is False, only compute the
upper triangle.

The computation can take a long time. Python won't check check for
a ^C until the function finishes. This can be irritating. Instead,
process only `batch_size` rows at a time before checking for a ^C.

The hits in the returned `SearchResults` are in arbitrary order.

Example::

    arena = chemfp.load_fingerprints("queries.fps")
    full_result = chemfp.search.threshold_tanimoto_search_symmetric(arena, threshold=0.2)
    upper_triangle = chemfp.search.threshold_tanimoto_search_symmetric(
              arena, threshold=0.2, include_lower_triangle=False)
    assert sum(map(len, full_result)) == sum(map(len, upper_triangle))*2
              
:param arena: the set of fingerprints
:type arena: a FingerprintArena
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:param include_lower_triangle:
    if False, compute only the upper triangle, otherwise use symmetry to compute the full matrix
:type include_lower_triangle: boolean
:param batch_size: the number of rows to process before checking for a ^C
:type batch_size: integer
:returns: a SearchResults instance


knearest_tanimoto_search_fp
---------------------------

.. py:method:: knearest_tanimoto_search_fp(query_fp, target_arena, k=3, threshold=0.7)

Search for `k`-nearest hits in `target_arena` which are at least `threshold` similar to `query_fp`

The hits in the `SearchResults` are ordered by decreasing similarity score.

Example::

    query_id, query_fp = chemfp.load_fingerprints("queries.fps")[0]
    targets = chemfp.load_fingerprints("targets.fps")
    print list(chemfp.search.knearest_tanimoto_search_fp(query_fp, targets, k=3, threshold=0.0))

:param query_fp: the query fingerprint
:type query_fp: a byte string
:param target_arena: the target arena
:type target_fp: a FingerprintArena
:param k: the number of nearest neighbors to find.
:type k: positive integer
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: a SearchResult

knearest_tanimoto_search_arena
------------------------------

.. py:method:: knearest_tanimoto_search_arena(query_arena, target_arena, k=3, threshold=0.7)

Search for the `k` nearest hits in the `target_arena` at least `threshold` similar to the fingerprints in `query_arena`

The hits in the `SearchResults` are ordered by decreasing similarity score.

Example::

    queries = chemfp.load_fingerprints("queries.fps")
    targets = chemfp.load_fingerprints("targets.fps")
    results = chemfp.search.knearest_tanimoto_search_arena(queries, targets, k=3, threshold=0.5)
    for query_id, query_hits in zip(queries.ids, results):
        if len(query_hits) >= 2:
            print query_id, "->", ", ".join(query_hits.get_ids())

:param query_arena: The query fingerprints.
:type query_arena: a FingerprintArena
:param target_arena: The target fingerprints.
:type target_arena: a FingerprintArena
:param k: the number of nearest neighbors to find.
:type k: positive integer
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:returns: a SearchResults instance

knearest_tanimoto_search_symmetric
----------------------------------

.. py:method:: knearest_tanimoto_search_symmetric(arena, k=3, threshold=0.7, batch_size=100)

Search for the `k`-nearest hits in the `arena` at least `threshold` similar to the fingerprints in the arena

The computation can take a long time. Python won't check check for
a ^C until the function finishes. This can be irritating. Instead,
process only `batch_size` rows at a time before checking for a ^C.

The hits in the `SearchResults` are ordered by decreasing similarity score.

Example::

    arena = chemfp.load_fingerprints("queries.fps")
    results = chemfp.search.knearest_tanimoto_search_symmetric(arena, k=3, threshold=0.8)
    for (query_id, hits) in zip(arena.ids, results):
        print query_id, "->", ", ".join(("%s %.2f" % hit) for hit in  hits.get_ids_and_scores())

:param arena: the set of fingerprints
:type arena: a FingerprintArena
:param k: the number of nearest neighbors to find.
:type k: positive integer
:param threshold: The minimum score threshold.
:type threshold: float between 0.0 and 1.0, inclusive
:param include_lower_triangle:
    if False, compute only the upper triangle, otherwise use symmetry to compute the full matrix
:type include_lower_triangle: boolean
:param batch_size: the number of rows to process before checking for a ^C
:type batch_size: integer
:returns: a SearchResults instance



SearchResults
=============

.. py:class:: SearchResults(... do not call directly ...)

Search results for a list of query fingerprints against a target arena

This acts like a list of SearchResult elements, with the ability
to iterate over each search results, look them up by index, and
get the number of scores.

In addition, there are helper methods to iterate over each hit and
to get the hit indicies, scores, and identifiers directly as Python
lists, sort the list contents, and more.

len(results)
------------

.. py:method:: __len__()

The number of rows in the SearchResults

results[i]
----------

.. py:method:: __getitem__(i)

Get the `i`th SearchResult

clear_all
---------

.. py:method:: clear_all()

Remove all hits from all of the search results

count_all
---------

.. py:method:: count_all(min_score=None, max_score=None, interval="[]")

Remove all hits from all of the search results

cumulative_score_all
--------------------

.. py:method:: cumulative_score_all(min_score=None, max_score=None, interval="[]")

The sum of all scores in all rows which are between `min_score` and `max_score`

Using the default parameters this returns the sum of all of
the scores in all of the results. With a specified range this
returns the sum of all of the scores in that range. The
cumulative score is also known as the raw score.

The default `min_score` of None is equivalent to -infinity.
The default `max_score` of None is equivalent to +infinity.

The `interval` parameter describes the interval end
conditions. The default of "[]" uses a closed interval,
where min_score <= score <= max_score. The interval "()"
uses the open interval where min_score < score < max_score.
The half-open/half-closed intervals "(]" and "[)" are
also supported.

:param min_score: the minimum score in the range.
:type min_score: a float, or None for -infinity
:param max_score: the maximum score in the range.
:type max_score: a float, or None for +infinity
:param interval: specify if the end points are open or closed.
:type interval: one of "[]", "()", "(]", "[)"
:returns: an floating point count
 
iter(results)
-------------

.. py:method:: __iter__()

Iterate over each SearchResult hit

iter_ids
--------

.. py:method:: iter_ids()

For each hit, yield the list of target identifiers

iter_ids_and_scores
-------------------

.. py:method:: iter_ids_and_scores()

For each hit, yield the list of (target id, score) tuples

iter_indices
------------

.. py:method:: iter_indices()

For each hit, yield the list of target indices

iter_indices_and_scores
-----------------------

.. py:method:: iter_indices_and_scores()

For each hit, yield the list of (target index, score) tuples

iter_scores
-----------

.. py:method:: iter_scores()

For each hit, yield the list of target scores

iter_hits
---------

REMOVED: Renamed to iter_ids_and_scores for 1.1.

reorder_all
-----------

.. py:method:: reorder_all()

Reorder the hits for all of the rows based on the requested `order`.

The available orderings are:
  increasing-score: sort by increasing score
  decreasing-score: sort by decreasing score
  increasing-index: sort by increasing target index
  decreasing-index: sort by decreasing target index
  move-closest-first: move the hit with the highest score to the first position
  reverse: reverse the current ordering

:param ordering: the name of the ordering to use



SearchResult
============

.. py:class:: SearchResult(... do not call directly ...)

Search results for a query fingerprint against a target arena.

The results contains a list of hits. Hits contain a target index,
score, and optional target ids. The hits can be reordered based on
score or index.

len(result)
------------

.. py:method:: __len__()

The number of hits

iter(result)
------------

.. py:method:: __iter__()

Iterate through the pairs of (target index, score) using the current ordering

clear
-----

.. py:method:: clear()

Remove all hits from this result

count
-----

.. py:method:: count(min_score=None, max_score=None, interval="[]")

Count the number of hits with a score between `min_score` and `max_score`

Using the default parameters this returns the number of
hits in the result.

The default `min_score` of None is equivalent to -infinity.
The default `max_score` of None is equivalent to +infinity.

The `interval` parameter describes the interval end
conditions. The default of "[]" uses a closed interval,
where min_score <= score <= max_score. The interval "()"
uses the open interval where min_score < score < max_score.
The half-open/half-closed intervals "(]" and "[)" are
also supported.

:param min_score: the minimum score in the range.
:type min_score: a float, or None for -infinity
:param max_score: the maximum score in the range.
:type max_score: a float, or None for +infinity
:param interval: specify if the end points are open or closed.
:type interval: one of "[]", "()", "(]", "[)"
:returns: an integer count

cumulative_score
----------------

.. py:method:: cumulative_score(min_score=None, max_score=None, interval="[]")

The sum of the scores which are between `min_score` and `max_score`

Using the default parameters this returns the sum of all of
the scores in the result. With a specified range this returns
the sum of all of the scores in that range. The cumulative
score is also known as the raw score.

The default `min_score` of None is equivalent to -infinity.
The default `max_score` of None is equivalent to +infinity.

The `interval` parameter describes the interval end
conditions. The default of "[]" uses a closed interval,
where min_score <= score <= max_score. The interval "()"
uses the open interval where min_score < score < max_score.
The half-open/half-closed intervals "(]" and "[)" are
also supported.

:param min_score: the minimum score in the range.
:type min_score: a float, or None for -infinity
:param max_score: the maximum score in the range.
:type max_score: a float, or None for +infinity
:param interval: specify if the end points are open or closed.
:type interval: one of "[]", "()", "(]", "[)"
:returns: an floating point count

get_ids
-------

.. py:method:: get_ids()

The list of target identifiers (if available), in the current ordering

get_ids_and_scores
------------------

.. py:method:: get_ids_and_scores()

The list of (target identifier, target score) pairs, in the current ordering

Raises a TypeError if the target IDs are not available.

get_indices
-----------

.. py:method:: get_indices()

The list of target indices, in the current ordering.

get_indices_and_scores
----------------------

.. py:method:: get_indices_and_scores()

The list of (target index, score) pairs, in the current ordering

get_scores
----------

.. py:method:: get_scores()

The list of target scores, in the current ordering

reorder
-------

.. py:method:: reorder(ordering="decreasing-score")

Reorder the hits based on the requested ordering.

The available orderings are:
  increasing-score: sort by increasing score
  decreasing-score: sort by decreasing score
  increasing-index: sort by increasing target index
  decreasing-index: sort by decreasing target index
  move-closest-first: move the hit with the highest score to the first position
  reverse: reverse the current ordering

:param ordering: the name of the ordering to use





.. _chemfp.bitops:

=====================
chemfp.bitopts module
=====================

.. py:module:: chemfp.bitops

The following functions are in the chemfp.bitops module. They
provide low-level bit operations on byte and hex fingerprints.


byte_popcount
=============

.. py:function:: byte_popcount()

byte_popcount(fp)

Return the number of bits set in a byte fingerprint

byte_intersect_popcount
=======================

.. py:function:: byte_intersect_popcount()

byte_intersect_popcount(fp1, fp2)

Return the number of bits set in the instersection of the two byte fingerprints

byte_tanimoto
=============

.. py:function:: byte_tanimoto()

byte_tanimoto(fp1, fp2)

Compute the Tanimoto similarity between two byte fingerprints

byte_contains
=============

.. py:function:: byte_contains()

byte_contains(super_fp, sub_fp)

Return 1 if the on bits of sub_fp are also 1 bits in super_fp

hex_isvalid
===========

.. py:function:: hex_isvalid()

hex_isvalid(s)

Return 1 if the string is a valid hex fingerprint, otherwise 0

hex_popcount
============

.. py:function:: hex_popcount()

hex_popcount(fp)

Return the number of bits set in a hex fingerprint, or -1 for non-hex strings

hex_intersect_popcount
======================

.. py:function:: hex_intersect_popcount()

hex_intersect_popcount(fp1, fp2)

Return the number of bits set in the intersection of the two hex fingerprint,
or -1 if either string is a non-hex string


hex_tanimoto
============

.. py:function:: hex_tanimoto()

hex_tanimoto(fp1, fp2)

Compute the Tanimoto similarity between two hex fingerprints.
Return a float between 0.0 and 1.0, or -1.0 if either string is not a hex fingerprint


hex_contains
============

.. py:function	:: hex_contains()

hex_contains(super_fp, sub_fp)

Return 1 if the on bits of sub_fp are also 1 bits in super_fp, otherwise 0.
Return -1 if either string is not a hex fingerprint