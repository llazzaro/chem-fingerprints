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

.. py:function:: open(source[, format=None])

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

.. _chemfp_read_structure_fingerprints:

read_structure_fingerprints
===========================

.. py:function:: load_fingerprints(reader[, metadata=None][, reorder=True])

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

:param reader: An iterator over (id, fingerprint) pairs
:type reader: a string, file object, or (id, fingerprint) iterator
:param metadata: The metadata for the arena, if other than reader.metadata
:type metadata: Metadata
:param reorder: Specify if fingerprints should be reordered for better performance
:type reorder: True or False
:returns: FingerprintArena

.. _chemfp_count_tanimoto_hits:

count_tanimoto_hits
===================

.. py:function:: count_tanimoto_hits(queries, targets[, threshold=0.7][, arena_size=100])

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

.. py:function:: threshold_tanimoto_search (queries, targets[, threshold=0.7][, arena_size=100])

Find all targets within 'threshold' of each query term

For each query in 'queries', find all the targets in 'targets' which
are at least 'threshold' similar to the query. This function returns
an iterator containing the (query_id, hits) pairs. The hits are stored
as a list of (target_id, score) pairs.

Example::

  queries = chemfp.open("queries.fps")
  targets = chemfp.load_fingerprints("targets.fps.gz")
  for (query_id, hits) in chemfp.threshold_tanimoto_search(queries, targets, threshold=0.8):
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

.. py:function:: knearest_tanimoto_search (queries, targets[, k=3][, threshold=0.7][, arena_size=100])

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
  for (query_id, hits) in chemfp.knearest_tanimoto_search(queries, targets, k=3, threshold=0.8):
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

.. py:class:: Metadata([num_bits=None][, num_bytes=None][, type=None][, aromaticity=None][, software=None][, sources=None][, date=None])

Store information about a set of fingerprints

num_bits = number of bits in the fingerprint
num_bytes = number of bytes in the fingerprint
type = fingerprint type
aromaticity = aromaticity model (only used with OEChem)
software = software used to make the fingerprints
sources = list of sources used to make the fingerprint
date = timestamp of when the fingerprints were made

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

.. py:method:: iter_arenas([arena_size=1000])

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

The following classes are returned as part of the public API but
should not be constructed directly.

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

len(arena)
----------

.. py:method:: __len__()

Number of fingerprint records in the FingerprintArena

arena[i]
--------

.. py:method:: __getitem__(i)

Return the (id, fingerprint) at position i


iter(arena)
-----------

.. py:method:: __iter__()

Iterate over the (id, fingerprint) contents of the arena


iter_arenas
-----------

.. py:method:: iter_arenas([arena_size=1000])

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


count_tanimoto_hits_fp
----------------------

.. py:method:: count_tanimoto_hits_fp(query_fp[, threshold=0.7])

Count the fingerprints which are similar enough to the query fingerprint

Return the number of fingerprints in this arena which are
at least `threshold` similar to the query fingerprint `query_fp`.

:param query_fp: query fingerprint
:type query_fp: byte string
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: integer count


count_tanimoto_hits_arena
-------------------------

.. py:method:: count_tanimoto_hits_arena(query_arena[, threshold=0.7])

Count the fingerprints which are similar enough to each query fingerprint

For each fingerprint in the `query_arena`, count the number of
fingerprints in this arena with Tanimoto similarity of at
least `threshold`. The resulting list order is the same as the
query fingerprint order.

:param query_fp: query arena
:type query_fp: FingerprintArena
:param threshold: minimum similarity threshold (default: 0.7)
:type threshold: float between 0.0 and 1.0, inclusive
:returns: list of integer counts

threshold_tanimoto_search_fp
----------------------------

.. py:method:: threshold_tanimoto_search_fp(query_fp[, threshold=0.7])

Find the fingerprints which are similar enough to the query fingerprint

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

.. py:method:: threshold_tanimoto_search_arena(query_arena[, threshold=0.7])

Find the fingerprints which are similar to each of the query fingerprints

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

.. py:method:: knearest_tanimoto_search_fp(query_fp[, k=3][, threshold=0.7])

Find the k-nearest fingerprints which are similar to the query fingerprint

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

.. py:method:: knearest_tanimoto_search_arena(query_arena[, k=3][, threshold=0.7])

Find the k-nearest fingerprint which are similar to each of the query fingerprints

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


save
----

.. py:method:: save(destination)

Save the arena contents to the given filename or file object



SearchResults
=============

.. py:class:: SearchResults(... do not call directly ...)

Contains the result of a Tanimoto threshold or k-nearest search

Each result contains a list of hits, where the hit is a
two-element tuple. If you iterate over the SearchResult then
you'll get the hits as (target_id, target_score) pairs.
tuples. If you iterate using the method `iter_hits()` then you'll
get the hits as (target_index, target_score) pairs.

iter(results)
-------------

.. py:method:: __iter__()

Iterate over the named hits for each result

Each term is a list of hits. A hit contains (id, score) tuples.
The order of the hits depends on the search algorithm.

iter_hits
---------

.. py:method:: iter_hits()

Iterate over the indexed hits for each result

Each term is a list of hits. A hit contains (index, score) tuples.
The order of the hits depends on the search algorithm.


len(results)
------------

.. py:method:: __len__()

Number of search results

results[i]
----------

.. py:method:: __getitem__(i)

The list of hits for result at position i

Each hit contains a (id, score) tuple.

size(i)
-------

.. py:method:: size()

The number of hits for result at position i

:param i: index into the search results
:type i: int
:returns: int


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