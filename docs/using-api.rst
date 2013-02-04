.. highlight:: python

=========================
The chemfp Python library
=========================

The chemfp command-line programs use a Python library called
chemfp. Portions of the API are in flux and subject to change. The
stable portions of the API which are open for general use are
documented in :ref:`chemfp-api`.

The API includes:

 - low-level Tanimoto and popcount operations
 - Tanimoto search algorithms based on threshold and/or k-nearest neighbors
 - a cross-toolkit interface for reading fingerprints from a structure file

The following chapters give examples of how to use the API.

Byte and hex fingerprints
=========================

In this section you'll learn how chemfp stores fingerprints and some
of the low-level bit operations on those fingerprints.

chemfp stores fingerprints as byte strings. Here are two 8 bit
fingerprints::

    >>> fp1 = "A"
    >>> fp2 = "B"

The :ref:`chemfp.bitops <chemfp.bitops>` module contains functions which work on byte
fingerprints. Here's the Tanimoto of those two fingerprints::

    >>> from chemfp import bitops
    >>> bitops.byte_tanimoto(fp1, fp2)
    0.33333333333333331

To understand why, you have to know that ASCII character "A" has the
value 65, and "B" has the value 66. The bit representation is::

     "A" = 01000001   and   "B" = 01000010

so their intersection has 1 bit and the union has 3, giving a Tanimoto
of 1/3 or 0.33333333333333331 when represented as a 64 bit floating
point number on the computer.

You can compute the Tanimoto between any two byte strings with the
same length, as in::

    >>> bitops.byte_tanimoto("apples&", "oranges")
    0.58333333333333337

You'll get a chemfp exception if they have different lengths.

.. highlight:: none

Most fingerprints are not as easy to read as the English ones I showed
above. They tend to look more like::


    P1@\x84K\x1aN\x00\n\x01\xa6\x10\x98\\\x10\x11

which is hard to read. I usually show hex-encoded fingerprints. The above
fingerprint in hex is::

    503140844b1a4e000a01a610985c1011

which is simpler to read, though you still need to know your hex
digits.

.. highlight:: python

The bitops module includes other low-level functions which work on
byte fingerprints, as well as corresponding functions which work on
hex fingerprints. (Hex-encoded fingerprints are decidedly second-class
citizens in chemfp, but they are citizens.)


Fingerprint collections and metadata
====================================

In this section you'll learn the basic operations on a fingerprint
collection and the fingerprint metadata.

A fingerprint record is the fingerprint plus an identifier. In chemfp,
a fingerprint collection is a object which contains fingerprint
records and which follows the common API providing access to those
records.

That's rather abstract, so let's work with a few real examples. You'll
need to create a copy of the "pubchem_targets.fps" file generated in
:ref:`pubchem_fingerprints` in order to follow along.

Here's how to open an FPS file::

    >>> import chemfp
    >>> reader = chemfp.open("pubchem_targets.fps")

Every fingerprint collection has a metadata attribute with details
about the fingerprints. It comes from the header of the FPS file. You
can view the metadata in Python repr format:

    >>> reader.metadata
    Metadata(num_bits=881, num_bytes=111, type='CACTVS-E_SCREEN/1.0 extend
    ed=2', aromaticity=None, sources=['Compound_014550001_014575000.sdf.gz
    '], software=u'CACTVS/unknown', date='2011-09-14T12:10:34')

but I think it's easier to view it in string format, which matches the
format of the FPS header:

    >>> print reader.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    

All fingerprint collections support iteration. Each step of the
iteration returns the fingerprint identifier and its score. Since I
know the 6th record has the id 14550045, I can write a simple loop
which stops with that record::

    >>> for (id, fp) in reader:
    ...   print id, "starts with", fp.encode("hex")[:20]
    ...   if id == "14550045":
    ...     break
    ... 
    14550001 starts with 034e1c00020000000000
    14550002 starts with 034e0c00020000000000
    14550003 starts with 034e0400020000000000
    14550005 starts with 010e1c00000600000000
    14550010 starts with 034e1c40000000000000
    14550045 starts with 071e8c03000000000000


Fingerprint collections also support iterating via arenas, and several
support Tanimoto search functions.


FingerprintArena
================

In this section you'll learn about the FingerprintArena fingerprint
collection and how to iterate through arenas in a collection.

The FPSReader reads through or searches a fingerprint file once. If
you want to read the file again you have to reopen it.

Reading from disk is slow, and the FPS format is designed for
ease-of-use and not performance. If you want to do many queries then
it's best to store everything in memory. The
:ref:`FingerprintArena <chemfp_arena_fingerprintarena>` is a
fingerprint collection which does that.

Here's how to load fingerprints into an arena::

    >>> import chemfp
    >>> arena = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> print arena.metadata
    #num_bits=881
    #type=CACTVS-E_SCREEN/1.0 extended=2
    #software=CACTVS/unknown
    #source=Compound_014550001_014575000.sdf.gz
    #date=2011-09-14T12:10:34
    

This implements the fingerprint collection API, so you can do things
like iterate over an arena and get the id/fingerprint pairs.::

    >>> from chemfp import bitops
    >>> for id, fp in arena:
    ...     print id, "with popcount", bitops.byte_popcount(fp)
    ...     if id == "14574718":
    ...         break
    ... 
    14550474 with popcount 2
    14574635 with popcount 2
    14550409 with popcount 4
    14550416 with popcount 6
    14574551 with popcount 7
    14550509 with popcount 8
    14550423 with popcount 10
    14550427 with popcount 10
    14574637 with popcount 10
    14574890 with popcount 11
    14574718 with popcount 12

If you look closely you'll notice that the fingerprint record order
has changed from the previous section, and that the population counts
are suspiciously non-decreasing. By default ref:`load_fingerprints`
reorders the fingerprints into a data structure which is faster to
search, although you can disable that if you want the fingerprints to
be the same as the input order.

The :ref:`FingerprintArena <chemfp_arena_fingerprintarena>` has new capabilities. You can ask it how many
fingerprints it contains, get the list of identifiers, and look up a
fingerprint record given an index, as in::

    >>> len(arena)
    3119
    >>> arena.ids[:5]
    ['14550474', '14574635', '14550409', '14550416', '14574551']
    >>> id, fp = arena[6]
    >>> id
    '14550423'
    >>> arena[-1][0]
    '14566760'
    >>> bitops.byte_popcount(arena[-1][1])
    231

An arena supports iterating through subarenas. This is like having a
long list and being able to iterate over sublists. Here's an example
of iterating over the arena to get subarenas of size 1000 (excepting
the last), and print information about each subarena.::

    >>> for subarena in arena.iter_arenas(1000):
    ...   print subarena.ids[0], len(subarena)
    ... 
    14550474 1000
    14573373 1000
    14555885 1000
    14560068 119
    >>> arena[0][0]
    '14550474'
    >>> arena[1000][0]
    '14573373'

To help demonstrate what's going on, I showed the first id of each
record along with the main arena ids for records 0 and 1000, so you
can verify that they are the same.

Arenas are a core part of chemfp. Processing one fingerprint at a time
is slow, so the main search routines expect to iterate over query
arenas, rather than query fingerprints.

Thus, the FPSReaders - and all chemfp fingerprint collections - also
support the `iter_arenas` interface. Here's an example of reading the
targets file 25 records at a time::

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for arena in queries.iter_arenas(25):
    ...   print len(arena)
    ... 
    25
    25
    25
    25
    25
    25
    25
    25
    24

Those add up to 224, which you can verify is the number of structures
in the original source file.

If you have a `FingerprintArena` then you can also use Python's slice
notation to make a subarena::

    >>> queries = chemfp.load_fingerprints("pubchem_queries.fps")
    >>> queries[10:15]
    <chemfp.arena.FingerprintArena object at 0x552c10>
    >>> queries[10:15].ids
    ['27599116', '27599118', '27599120', '27583411', '27599082']
    >>> queries.ids[10:15]
    ['27599116', '27599118', '27599120', '27583411', '27599082']

The big restriction is that slices can only have a step size
of 1. Slices like `[10:20:2]` and `[::-1]` aren't supported. If you
want something like that then you'll need to make a new arena instead
of using a subarena slice.

In case you were wondering, yes, you can use `iter_arenas` or the other
FingerprintArena methods on a subarena::

    >>> queries[10:15][1:3].ids
    ['27599118', '27599120']
    >>> queries.ids[11:13]
    ['27599118', '27599120']




How to use query fingerprints to search for similar target fingerprints
=======================================================================

In this section you'll learn how to do a Tanimoto search using the
previously created PubChem fingerprint files for the queries and the
targets.

It's faster to search an arena, so I'll load the target fingerprints:

    >>> import chemfp
    >>> targets = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> len(targets)
    3119

and open the queries as an FPSReader.

    >>> queries = chemfp.open("pubchem_queries.fps")

I'll use :ref:`threshold_tanimoto_search <chemfp_threshold_tanimoto_search>` to find, for each query,
all hits which are at least 0.7 similar to the query.

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, hits) in chemfp.threshold_tanimoto_search(queries, targets, threshold=0.7):
    ...   print query_id, len(hits), hits[:2]
    ... 
    27575433 0 []
    27575577 18 [('14570945', 0.74874371859296485), ('14570946', 0.73762376237623761)]
    27575602 3 [('14572463', 0.72560975609756095), ('14553070', 0.75935828877005351)]
    27575603 3 [('14572463', 0.72560975609756095), ('14553070', 0.75935828877005351)]
    27575880 9 [('14569876', 0.72307692307692306), ('14567856', 0.73076923076923073)]
    27575897 0 []
    27577227 1 [('14570135', 0.7142857142857143)]
    27577234 0 []
          # ... many lines omitted ...

I'm only showing the first two hits for the sake of space. It seems
rather pointless, after all, to show all 18 hits of query id 27575577.

What you don't see is that the implementation uses the iter_arenas()
interface on the queries so that it processes only a subarena at a
time. There's a tradeoff between a large arena, which is faster
because it doesn't often go back to Python code, or a small arena,
which uses less memory and is more responsive. You can change the
tradeoff using the `arena_size` parameter.


If all you care about is the count of the hits within a given
threshold then use :ref:`chemfp.count_tanimoto_hits <chemfp_count_tanimoto_hits>`

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, count) in chemfp.count_tanimoto_hits(queries, targets, threshold=0.7):
    ...     print query_id, count
    ...   break
    ... 
    27575433 0
    27575577 18
    27575602 3
    27575603 3
    27575880 9
    27575897 0
    27577227 1
    27577234 0
    27577237 1
    27577250 4
         # ... many lines omitted ...

Or, if you only want the k=2 nearest neighbors to each target within
that same threshold of 0.7 then use ref:`chemfp.knearest_tanimoto_search`::

    >>> queries = chemfp.open("pubchem_queries.fps")
    >>> for (query_id, hits) in chemfp.knearest_tanimoto_search(query_arena, targets, k=2, threshold=0.7):
    ...     print query_id, hits
    ...   break
    ... 
    27575433 []
    27575577 [('14570945', 0.74874371859296485), ('14570951', 0.73853211009174313)]
    27575602 [('14553070', 0.75935828877005351), ('14572463', 0.72560975609756095)]
    27575603 [('14553070', 0.75935828877005351), ('14572463', 0.72560975609756095)]
    27575880 [('14569866', 0.77272727272727271), ('14567856', 0.73076923076923073)]
    27575897 []
    27577227 [('14570135', 0.7142857142857143)]
    27577234 []
    27577237 [('14569555', 0.73711340206185572)]
    27577250 [('14569555', 0.74742268041237114), ('14550456', 0.72131147540983609)]
         # ... many lines omitted ...



How to search an FPS file
=========================

In this section you'll learn how to search an FPS file directly,
without loading it into a FingerprintArena.

The previous example loaded the fingerprints into a
FingerprintArena. That's the fastest way to do multiple
searches. Sometimes though you only want to do one or a couple of
queries. It seems rather excessive to read the entire targets file
into an in-memory data structure before doing the search when you
could search will processing the file.

For that case, use an FPSReader as the target file. Here I'll get the
first record from the queries file and use it to search the targets
file::

    >>> query_arena = next(chemfp.open("pubchem_queries.fps").iter_arenas(1))

This line opens the file, iterates over its fingerprint records, and
return the first one.

(Note: the next() function was added after Python 2.5 so the above
won't work for that version. Instead, use::

    >>> query_arena = chemfp.open("pubchem_queries.fps").iter_arenas(1).next()

which is the older form. Or you can use the equally bewildering

    >>> for query_arena in chemfp.open("pubchem_queries.fps").iter_arenas(1):
    ...   break
.)

Here are the k=5 closest hits against the targets file::

    >>> targets = chemfp.open("pubchem_targets.fps")
    >>> for query_id, hits in chemfp.knearest_tanimoto_search(query_arena, targets, k=5, threshold=0.0):
    ...   print "Hits for", query_id
    ...   for hit in hits:
    ...     print "", hit
    ... 
    Hits for 27575433
     ('14568234', 0.69035532994923854)
     ('14550456', 0.64921465968586389)
     ('14572463', 0.64444444444444449)
     ('14566364', 0.63953488372093026)
     ('14573723', 0.63247863247863245)

Remember that the FPSReader is based on reading an FPS file. Once
you've done a search, the file is read, and you can't do another
search. You'll need to reopen the file.

Each search processes arena_size query fingerprints at a time. You
will need to increase that value if you want to search more than that
number of fingerprints with this method. The search performance
tradeoff between a FPSReader search and loading the fingerprints into
a FingerprintArena occurs with under 10 queries, so there should be
little reason to worry about this.


FingerprintArena searches returning indices instead of ids
===========================================================

In this section you'll learn how to search a FingerprintArena and use
hits based on integer indices rather than string ids.

The previous sections used a high-level interface to the Tanimoto
search code. Those are designed for the common case where you just
want the query id and the hits, where each hit includes the target id.

Working with strings is actually rather inefficient in both speed and
memory. It's usually better to work with indices if you can, and in
the next section I'll show how to make a distance matrix using this
interface.

NOTE: up until the final 1.1 release, this document said to use the
FingerprintArena methods. This is no longer recommended.  Use the
`chemfp.search` functions instead. Most of the search methods, except
perhaps the single fingerprint methods, will issue a
DeprecationWarning with 1.2 and be removed with 1.3.

The index-based search functions are in the `chemfp.search` module.
They can be categorized into three groups:

  1. Count the number of hits:
    * count_tanimoto_hits_fp - search an arena using a single fingerprint
    * count_tanimoto_hits_arena - search an arena using an arena
    * count_tanimoto_hits_symmetric - search an arena using itself
    * partial_count_tanimoto_hits_symmetric - (advanced use; see the doc string)

  2. Find all hits at or above a given threshold, sorted arbitrarily:
    * threshold_tanimoto_search_fp - search an arena using a single fingerprint
    * threshold_tanimoto_search_arena - search an arena using an arena
    * threshold_tanimoto_search_symmetric - search an arena using itself
    * partial_threshold_tanimoto_search_symmetric - (advanced use; see the doc string)
    * fill_lower_triangle - copy the upper triangle terms to the lower triangle

  3. Find the k-nearest hits at or above a given threshold, sorted by decreasing similarity:
    * knearest_tanimoto_search_fp - search an arena using a single fingerprint
    * knearest_tanimoto_search_arena - search an arena using an arena
    * knearest_tanimoto_search_symmetric - search an arena using itself

The functions ending '_fp' take a query fingerprint and a target
arena. The functions ending '_arena' take a query arena and a target
arena. The functions ending '_symmetric' use the same arena as both
the query and target.

In the following example, I'll use the first 5 fingerprints of a data
set to search the entire data set. To do this, I load the data set as
an arena, extract the first 5 records as a sub-arena, and do the
search.

    >>> import chemfp
    >>> from chemfp import search
    >>> targets = chemfp.load_fingerprints("pubchem_queries.fps")
    >>> queries = targets[:5]
    >>> results = search.threshold_tanimoto_search_arena (queries, targets, threshold=0.7)

The threshold_tanimoto_search_arena search finds the target
fingerprints which have a similarity score of at least 0.7 compared to
the query.

You can iterate over the results to get the list of hits for each of
the queries. The order of the results is the same as the order of the
records in the query.::

    >>> for hits in results:
    ...   print len(hits), hits.get_ids_and_scores()[:3]
    ... 
    2 [('27581954', 1.0), ('27581957', 1.0)]
    2 [('27581954', 1.0), ('27581957', 1.0)]
    3 [('27580389', 1.0), ('27580394', 0.88235294117647056), ('27581637', 0.75)]
    2 [('27584917', 1.0), ('27585106', 0.89915966386554624)]
    2 [('27584917', 0.89915966386554624), ('27585106', 1.0)]


This result is like what you saw earlier, except that it doesn't have
the query id. You can get that from the arena's `id` attribute, which
contains the list of fingerprint identifiers.

    >>> for query_id, hits in zip(queries.ids, results):
    ...   print "Hits for", query_id
    ...   for hit in hits.get_ids_and_scores()[:3]:
    ...     print "", hit
    Hits for 27581954
     ('27581954', 1.0)
     ('27581957', 1.0)
    Hits for 27581957
     ('27581954', 1.0)
     ('27581957', 1.0)
         ...


What I really want to show is that you can get the same data only
using the offset index for the target record instead of its id. The
result from a Tanimoto search is a `SearchResults` object, with the
methods `get_indices_and_scores()`::

    >>> for hits in results:
    ...   print len(hits), hits.get_indices_and_scores()[:3]
    ... 
    2 [(0, 1.0), (1, 1.0)]
    2 [(0, 1.0), (1, 1.0)]
    3 [(2, 1.0), (5, 0.88235294117647056), (20, 0.75)]
    2 [(3, 1.0), (4, 0.89915966386554624)]
    2 [(3, 0.89915966386554624), (4, 1.0)]
    >>> 
    >>> targets.ids[0]
    '27581954'
    >>> targets.ids[1]
    '27581957'
    >>> targets.ids[5]
    '27580394'

I did a few id lookups given the target dataset to show you that the
index corresponds to the identifiers from the previous code.

These examples iterated over each `SearchResult` to fetch the ids and
scores, or indices and scores. Another possibility is to ask the
`SearchResults` to iterate directly over the list of fields you want.

    >>> for row in results.iter_indices_and_scores():
    ...   print len(row), row[:3]
    ... 
    2 [(0, 1.0), (1, 1.0)]
    2 [(0, 1.0), (1, 1.0)]
    3 [(2, 1.0), (5, 0.88235294117647056), (20, 0.75)]
    2 [(3, 1.0), (4, 0.89915966386554624)]
    2 [(3, 0.89915966386554624), (4, 1.0)]

This was added to get a bit more performance out of chemfp and because
the API is sometimes cleaner one way and sometimes cleaner than the
other. Yes, I know that the Zen of Python recommends that "there
should be one-- and preferably only one --obvious way to do it." Oh
well.

NOTE: The API has changed slightly from 1.0 to 1.1. Previously the
`SearchResults` had the methods `iter_hits` and iteration over the
`SearchResult` returned a "hit." However, I couldn't remember if a hit
used the identifier or the index. You must now be explicit and use
`iter_ids*` or `iter_indices*` on the `SearchResults`, and use
`get_ids*` or `get_indices*` on the `SearchResult`.



Computing a distance matrix for clustering
==========================================

In this section you'll learn how to compute a distance matrix using
the chemfp API.

chemfp does not do clustering. There's a huge number of tools which
already do that. A goal of chemfp in the future is to provide some
core components which clustering algorithms can use.

That's in the future. Right now you can use the following to build a
distance matrix and pass that to one of those tools.

Since we're using the same fingerprint arena for both queries and
targets, we know the distance matrix will be symmetric along the
diagonal, and the diagonal terms will be 1.0. The
`threshold_tanimoto_search_symmetric` functions can take advantage of
the symmetry for a factor of two performance gain. There's also a way
to limit it to just the upper triangle, which gives a factor of two
memory gain as well.


Most of those tools use `NumPy <http://numpy.scipy.org/>`_, which is a
popular third-party package for numerical computing. You will need to
have it installed for the following to work.

::

    import numpy  # NumPy must be installed
    from chemfp import search
    
    # Compute distance[i][j] = 1-Tanimoto(fp[i], fp[j])
    
    def distance_matrix(arena):
        n = len(arena)
        
        # Start off a similarity matrix with 1.0s along the diagonal
        similarities = numpy.identity(n, "d")
        
        ## Compute the full similarity matrix.
        # The implementation computes the upper-triangle then copies
        # the upper-triangle into lower-triangle. It does not include
        # terms for the diagonal.
        results = search.threshold_tanimoto_search_symmetric(arena, threshold=0.0)
        
        # Copy the results into the NumPy array.
        for row_index, row in enumerate(results.iter_indices_and_scores()):
            for target_index, target_score in row:
                similarities[row_index, target_index] = target_score

        # Return the distance matrix using the similarity matrix
        return 1.0 - similarities


Once you've computed the distance matrix, clustering is easy. I
installed the `hcluster <http://code.google.com/p/scipy-cluster/>`_
package, as well as `matplotlib <http://matplotlib.sourceforge.net/>`_,
then ran the following to see the hierarchical clustering::

    import chemfp
    import hcluster # Clustering package from http://code.google.com/p/scipy-cluster/
    
    # ... insert the 'distance_matrix' function definition here ...

    dataset = chemfp.load_fingerprints("pubchem_queries.fps")
    distances  = distance_matrix(dataset)
    
    linkage = hcluster.linkage(distances, method="single", metric="euclidean")
    
    # Plot using matplotlib, which you must have installed
    hcluster.dendrogram(linkage, labels=dataset.ids)
    
    import pylab
    pylab.show()



Taylor-Butina clustering
========================

For the last clustering example, here's my (non-validated) variation
of the `Butina algorithm from JCICS 1999, 39, 747-750 <http://www.chemomine.co.uk/dbclus-paper.pdf>`_.
See also http://www.redbrick.dcu.ie/~noel/R_clustering.html . You
might know it as Leader clustering.


First, for each fingerprint find all other fingerprints with a
threshold of 0.8::

    import chemfp
    from chemfp import search
    
    arena = chemfp.load_fingerprints("pubchem_targets.fps")
    results = search. threshold_tanimoto_search_symmetric (arena, threshold = 0.8)


Sort the results so that fingerprints with more hits come first. This
is more likely to be a cluster centroid. Break ties arbitrarily by the
fingerprint id; since fingerprints are ordered by the number of bits
this likely makes larger structures appear first.::

    # Reorder so the centroid with the most hits comes first.
    # (That's why I do a reverse search.)
    # Ignore the arbitrariness of breaking ties by fingerprint index
    results = sorted( (  (len(indices), i, indices)
                              for (i,indices) in enumerate(results.iter_indices())  ),
                      reverse=True)


Apply the leader algorithm to determine the cluster centroids and the singletons::


    # Determine the true/false singletons and the clusters
    true_singletons = []
    false_singletons = []
    clusters = []
    
    seen = set()
    for (size, fp_idx, members) in results:
        if fp_idx in seen:
            # Can't use a centroid which is already assigned
            continue
        seen.add(fp_idx)
    
        # Figure out which ones haven't yet been assigned
        unassigned = set(members) - seen
    
        if not unassigned:
            false_singletons.append(fp_idx)
            continue
            
        # this is a new cluster
        clusters.append( (fp_idx, unassigned) )
        seen.update(unassigned)

Once done, report the results::

    print len(true_singletons), "true singletons"
    print "=>", " ".join(sorted(arena.ids[idx] for idx in true_singletons))
    print
    
    print len(false_singletons), "false singletons"
    print "=>", " ".join(sorted(arena.ids[idx] for idx in false_singletons))
    print
    
    # Sort so the cluster with the most compounds comes first,
    # then by alphabetically smallest id
    def cluster_sort_key(cluster):
        centroid_idx, members = cluster
        return -len(members), arena.ids[centroid_idx]
        
    clusters.sort(key=cluster_sort_key)
    
    print len(clusters), "clusters"
    for centroid_idx, members in clusters:
        print arena.ids[centroid_idx], "has", len(members), "other members"
        print "=>", " ".join(arena.ids[idx] for idx in members)


The algorithm is quick for this small data set.

Out of curiosity, I tried this on 100,000 compounds selected
arbitrarily from PubChem. It took 35 seconds on my desktop (a 3.2 GHZ
Intel Core i3) with a threshold of 0.8. In the Butina paper, it took
24 hours to do the same, although that was with a 1024 bit fingerprint
instead of 881. It's hard to judge the absolute speed differences of a
MIPS R4000 from 1998 to a desktop from 2011, but it's less than the
factor of about 2000 you see here.

More relevent is the comparison between these numbers for the 1.1
release compared to the original numbers for the 1.0 release. On my
old laptop, may it rest it peace, it took 7 minutes to compute the
same benchmark. Where did the roughly 16-fold peformance boost come
from? Money. After 1.0 was released, Roche funded me to add various
optimizations, including taking advantage of the symmetery (2x) and
using hardware POPCNT if available (4x). Roche and another company
helped fund the OpenMP support, and my desktop ran this with 4 cores
instead of 1.

The wary among you might notice that 2*4*4 = 32x faster, while I
said the overall code was only 16x faster. Where's the factor of 2x
slowdown? It's in the Python code! The
`threshold_tanimoto_search_symmetric` step took only 13 seconds. The
remaining 22 seconds was in the leader code written in Python. To
make the analysis more complicated, improvements to the chemfp API
sped up the clustering step by about 40%.

With chemfp 1.0 version, the clustering performance overhead was minor
compared to the full similarity search, so I didn't keep track of
it. With chemfp 1.1, those roles have reversed! 


Reading structure fingerprints using a toolkit
==============================================

In this section you'll learn how to use a chemistry toolkit in order
to compute fingerprints from a given structure file.

What happens if you're given a structure file and you want to find the
two nearest matches in an FPS file? You'll have to generate the
fingerprints for the structures in the structure file, then do the
comparison.

For this section you'll need to have a chemistry toolkit. I'll use the
"chebi_maccs.fps" file you generated earlier as the targets, and the
PubChem file "Compound_027575001_027600000.sdf.gz as the source of
query structures.::

    >>> import chemfp
    >>> from chemfp import search
    >>> targets = chemfp.load_fingerprints("chebi_maccs.fps")
    >>> queries = chemfp.read_structure_fingerprints(targets.metadata, "Compound_027575001_027600000.sdf.gz")
    >>> for (query_id, hits) in chemfp.knearest_tanimoto_search(queries, targets, k=2, threshold=0.4):
    ...   print query_id, "=>",
    ...   for (target_id, score) in hits.get_ids_and_scores():
    ...     print "%s %.3f" % (target_id, score),
    ...   print
    ... 
    27575433 => CHEBI:280152 0.667 CHEBI:3176 0.662
    27575577 => CHEBI:6375 0.600 CHEBI:46068 0.600
    27575602 => CHEBI:3090 0.683 CHEBI:6790 0.635
    27575603 => CHEBI:3090 0.683 CHEBI:6790 0.635
    27575880 => CHEBI:59736 0.725 CHEBI:8887 0.617
    27575897 => CHEBI:8887 0.632 CHEBI:51491 0.622
    27577227 => CHEBI:59007 0.831 CHEBI:59120 0.721
    27577234 => CHEBI:59007 0.809 CHEBI:9398 0.722
    27577237 => CHEBI:59007 0.789 CHEBI:52890 0.741
    27577250 => CHEBI:59007 0.753 CHEBI:4681 0.722
         # ... many lines omitted ...

That's it! Pretty simple, wasn't it? You didn't even need to explictly
specify which toolkit you wanted to use.

The only new thing here is :ref:`read_structure_fingerprints <chemfp_read_structure_fingerprints>`. The
first parameter of this is the metadata used to configure the
reader. In my case it's::

    >>> print targets.metadata
    #num_bits=166
    #type=OpenEye-MACCS166/1
    #software=OEGraphSim/1.0.0 (20100809)
    #aromaticity=openeye
    #source=ChEBI_lite.sdf.gz
    #date=2011-09-14T17:50:28

The "type" told chemfp which toolkit to use to read molecules, and how
to generate fingerprints from those molecules, while "aromaticity"
told it which aromaticity model to use when reading the molecule file.

You can of course pass in your own metadata as the first parameter to
read_structure_fingerprints, and as a shortcut, if you pass in a
string then it will be used as the fingerprint type.

For examples, if you have OpenBabel installed then you can do::

   >>> reader = chemfp.read_structure_fingerprints("OpenBabel-MACCS", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 800404000840549e848189cca1f132aedfab6eff1b
    27575577 800400000000449e850581c22190022f8a8baadf1b
    27575602 000000000000449e840191d820a0122eda9abaff1b
    27575603 000000000000449e840191d820a0122eda9abaff1b

If you have OEChem and OEGraphSim installed then you can do::

    >>> reader = chemfp.read_structure_fingerprints("OpenEye-MACCS166", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 000000080840448e8481cdccb1f1b216daaa6a7e3b
    27575577 000000080000448e850185c2219082178a8a6a5e3b
    27575602 000000080000448e8401d14820a01216da983b7e3b
    27575603 000000080000448e8401d14820a01216da983b7e3b

And if you have RDKit installed then you can do::

    >>> reader = chemfp.read_structure_fingerprints("RDKit-MACCS166", "Compound_027575001_027600000.sdf.gz")
    >>> for i, (id, fp) in enumerate(reader):
    ...   print id, fp.encode("hex")
    ...   if i == 3:
    ...     break
    ... 
    27575433 000000000840549e84818dccb1f1323cdfab6eff1f
    27575577 000000000000449e850185c22190023d8a8beadf1f
    27575602 000000000000449e8401915820a0123eda98bbff1f
    27575603 000000000000449e8401915820a0123eda98bbff1f


Select a random fingerprint sample
==================================

In this section you'll learn how to make a new arena where the
fingerprints are randomly selected from the old arena.

A FingerprintArena slice creates a subarena. Technically speacking,
this is a "view" of the original data. The subarena doesn't actually
copy its fingerprint data from the original arena. Instead, it uses
the same fingerprint data, but keeps track of the start and end
position of the range it needs. This is why it's not possible to slice
with a step size other than +1.

This also means that memory for a large arena won't be freed until
all of its subarenas are also removed.

You can see some evidence for this because a `FingerprintArena` stores
the entire fingerprint data as a set of bytes named `arena`::

    >>> import chemfp
    >>> targets = chemfp.load_fingerprints("pubchem_targets.fps") 
    >>> subset = targets[10:20]
    >>> targets.arena is subset.arena
    True

This shows that the `targets` and `subset` share the same raw data
set. At least it does to me, the person who wrote the code.

You can ask an arena or subarena to make a `copy`_. This allocates
new memory for the new arena and copies all of its fingerprints there.

    >>> new_subset = subset.copy()
    >>> len(new_subset) == len(subset)
    >>> new_subset.arena is subset.arena
    False
    >>> subset[7][0]
    '14554484'
    >>> new_subset[7][0]
    '14554484'


The `copy` method can do more than just copy the arena. You can give
it a list of indices and it will only copy those fingerprints::

    >>> three_targets = targets.copy([3112, 0, 1234])
    >>> three_targets.ids
    ['14550474', '14564466', '14564904']
    >>> [targets.ids[3112], targets.ids[0], targets.ids[1234]]
    ['14564904', '14550474', '14564466']

Are you confused about why the identifiers aren't in the same order?
That's because when you specify indicies, the copy automatically
reorders them by popcount and stores the popcount information. This
extra work help makes future searches faster. Use `reorder=False`
to leave the order unchanged

   >>> my_ordering = targets.copy([3112, 0, 1234], reorder=False)
   >>> my_ordering.ids
   ['14564904', '14550474', '14564466']

This interesting, in a boring sort of way. Let's get back to the main
goal of getting a random subset of the data. I want to select `m`
records at random, without replacement, to make a new data set. You
can see this just means making a list with `m` different index
values. Python's built-in `random.sample` function makes this easy::

    >>> import random
    >>> random.sample("abcdefgh", 3)
    ['b', 'h', 'f']
    >>> random.sample("abcdefgh", 2)
    ['d', 'a']
    >>> random.sample([5, 6, 7, 8, 9], 2)
    [7, 9]
    >>> help(random.sample)
    sample(self, population, k) method of random.Random instance
       Chooses k unique random elements from a population sequence.
       ...
       To choose a sample in a range of integers, use xrange as an argument.
       This is especially fast and space efficient for sampling from a
       large population:   sample(xrange(10000000), 60)

The last line of the help points out what do next!::

    >>> random.sample(xrange(len(targets)), 5)
    [610, 2850, 705, 1402, 2635]
    >>> random.sample(xrange(len(targets)), 5)
    [1683, 2320, 1385, 2705, 1850]

Putting it all together, and here's how to get a new arena containing
100 randomly selected fingerprints, without replacement, from the
`targets` arena::

    >>> sample_indices = random.sample(xrange(len(targets)), 100)
    >>> sample = targets.copy(indices=sample_indices)
    >>> len(sample)
    100


Look up a fingerprint with a given id
=====================================

In this section you'll learn how to get a fingerprint record with a
given id.

All fingerprint records have an identifier and a
fingerprint. Identifiers should be unique. (Duplicates are allowed, and
if they exist then the lookup code described in this section will
arbitrarily decide which record to return. Once made, the choice will
not change.)

Let's find the fingerprint for the record in "pubchem_targets.fps"
which has the identifier `14564126`. One solution is to iterate
over all of the records in a file, using the FPS reader::

    >>> import chemfp
    >>> for id, fp in chemfp.open("pubchem_targets.fps"):
    ...   if id == "14564126":
    ...     break
    ... else:
    ...   raise KeyError("%r not found" % (id,))
    ... 
    >>> fp[:5]
    '\x07\x1e\x1c\x00\x00'

I used the somewhat obscure `else` clause to the `for` loop. If the
`for` finishes without breaking, which would happen if the identifier
weren't present, then it will raise an exception saying that it
couldn't find the given identifier.

If the fingerprint records are already in a `FingerprintArena` then
there's a better solution. Use the `get_fingerprint_with_id` method to
get the fingerprint byte string, or `None` if the identifier doesn't exist::

    >>> arena = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> fp = arena.get_fingerprint_by_id("14564126")
    >>> fp[:5]
    '\x07\x1e\x1c\x00\x00'
    >>> missing_fp = arena.get_fingerprint_by_id("does-not-exist")
    >>> missing_fp
    >>> missing_fp is None
     True

Internally this does about what you think it would. It uses the
arena's `id` list to make a lookup table mapping identifier to
index, and caches the table for later use. Given the index, it's very
easy to get the fingerprint.

In fact, you can get the index and do the record lookup yourself::

    >>> fp_index = arena.get_index_by_id("14564126")
    >>> arena.get_index_by_id("14564126")
     1559
    >>> arena[1559]
     ('14564126', '\x07\x1e\x1c\x00\x00 ...')


Sorting search results
======================

In this section you'll learn how to sort the search results.

The k-nearest searches return the hits sorted from highest score to
lowest, and break ties arbitrarily. This is usually what you want, and
the extra cost to sort is small (k*log(k)) compared to the time needed
to maintain the internal heap (N*log(k)).

By comparison, the threshold searches return the hits in arbitrary
order. Sorting takes up to N*log(N) time, which is extra work for
those cases where you don't want sorted data. Use the `reorder` method
of a `SearchResult` if you want the hits sorted in-place::

    >>> import chemfp
    >>> arena = chemfp.load_fingerprints("pubchem_queries.fps")
    >>> query_fp = arena.get_fingerprint_by_id("27599116")
    >>> from chemfp import search
    >>> result = search.threshold_tanimoto_search_fp(query_fp, arena, threshold=0.90)
    >>> len(result)
    6
    >>> result.get_ids_and_scores()
    [('27599092', 0.96153846153846156), ('27599115', 1.0), ('27599116', 1.0),
    ('27599118', 1.0), ('27599120', 1.0), ('27599082', 0.92537313432835822)]

    >>> result.reorder("decreasing-score")
    >>> result.get_ids_and_scores()
    [('27599115', 1.0), ('27599116', 1.0), ('27599118', 1.0), ('27599120', 1.0),
    ('27599092', 0.96153846153846156), ('27599082', 0.92537313432835822)]

    >>> result.reorder("increasing-score")
    >>> result.get_ids_and_scores()
    [('27599082', 0.92537313432835822), ('27599092', 0.96153846153846156),
    ('27599115', 1.0), ('27599116', 1.0), ('27599118', 1.0), ('27599120', 1.0)]

There are currently six different sort methods, all specified by
name. These are

      * increasing-score: sort by increasing score
      * decreasing-score: sort by decreasing score
      * increasing-index: sort by increasing target index
      * decreasing-index: sort by decreasing target index
      * reverse: reverse the current ordering
      * move-closest-first: move the hit with the highest score to the first position

The first two should be obvious from the examples. If you find
something useful for the next two then let me know. The `reverse`
option reverses the current ordering, and is most useful if you want
to reverse the sorted results from a k-nearest search.

The `move-closest-first` option exists to improve the leader algorithm
stage used by the Taylor-Butina algorithm. The newly seen compound is
either in the same cluster as its nearest neighbor or it is the new
centroid. I felt it best to implement this as a special reorder term,
rather than one of the other possible options.

If you are interested in other ways to help improve your clustering
performance, let me know.

Each `SearchResult` has a `reorder` method. If you want to reorder all
of the hits of a `SearchResults` then use its `reorder_all`
method::

    >>> similarity_matrix = search.threshold_tanimoto_search_symmetric(
    ...                         arena, threshold=0.8)
    >>> for query_id, row in zip(arena.ids, similarity_matrix):
    ...   print query_id, "->", row.get_ids_and_scores()[:3]
    ... 
    27581954 -> [('27581957', 1.0)]
    27581957 -> [('27581954', 1.0)]
    27580389 -> [('27580394', 0.88235294117647056)]
    27584917 -> [('27585106', 0.89915966386554624)]
    27585106 -> [('27584917', 0.89915966386554624)]
    27580394 -> [('27580389', 0.88235294117647056)]
    27593061 -> []
           ...

It takes the same set of ordering names as `reorder`_.



Working with raw scores and counts in a range
=============================================

In this section you'll learn how to get the hit counts and raw scores
for a interval.

The length of the `SearchResult` is the number of hits it contains::

    >>> import chemfp
    >>> from chemfp import search
    >>> arena = chemfp.load_fingerprints("pubchem_targets.fps")
    >>> fp = arena.get_fingerprint_by_id("14564126")
    >>> result = search.threshold_tanimoto_search_fp(fp, arena, threshold=0.2)
    >>> len(result)
    2836

This gives you the number of hits at or above a threshold of 0.2,
which you can also get by doing
`count_threhsold_tanimoto_search_fp`_. The result also stores the
hits, and you can get the number of hits which are within a specified
interval. Here are the hits counts at or above 0.5, 0.80, and 0.95::

    >>> result.count(0.5)
    735
    >>> result.count(0.8)
    5
    >>> result.count(0.95)
    2

The first parameter, `min_score`, specifies the minimum
threshold. The second, `max_score`, specifies the maximum. Here's
how to get the number of hits with a score of at most 0.95 and 0.5::

    >>> result.count(max_score=0.95)
    2834
    >>> result.count(max_score=0.5)
    2118

If you work out the math, you add 2118+735 and realize that
2853!=2836. There's a difference of 17. This is because the default
interval uses a closed range, and there are 17 hits with a score of
exactly 0.5::

    >>> result.count(0.5, 0.5)
    17

The third parameter, `interval`, specifies the end conditions. The
default is "[]" which means that both ends are closed. The interval
"()" means that both ends are open, and "[)" and "(]" are the two
half-open/half-closed ranges. To get the number of hits below 0.5 and
the number of hits at or above 0.5 then you might use:

    >>> result.count(None, 0.5, "[)")
    2101
    >>> result.count(0.5, None, "[]")
    735

at get the expected results. (A min or max of `None` means that there
is respectively no lower or no upper bound.)


Now for something a bit fancier. Suppose you have two sets of
structures. How well do they compare to each other? I can think of
various ways to do it. One is to look at a comparison profile. Find
all NxM comparisons between the two sets. How many of the hits have a
threshold of 0.2? How many at 0.5? 0.95?

If there are "many", then the two sets are likely more similar than
not. If the answer is "few", then they are likely rather distinct.

I'll be more specific. Are the coenzyme A-like structures in ChEBI
more similar to the penicillin-like structures than you would expect
by comparing two randomly chosen subsets? By similar, I'll use
Tanimoto similarity of the "chebi_maccs.fps" file created in the
`Generating fingerprints with ...` command-line tool example.

The CHEBI id for coenzyme A is CHEBI:15346 and for penicillin is
CHEBI:17334. I'll define the "coenzyme A-like" structures as the 117
structures where the fingerprint is at least 0.95 similar to coenzyme
A, and "penicillin-like" as the 15 structures at least 0.90 similar to
penicillin. This gives 1755 total comparisons.

You know enough to do this, but there's a nice optimization I haven't
told you about. You can get the total count of all of the threshold
hits using the `SearchResults.count_all` method, instead of looping
over each `SearchResult` and calling `count`_::

    import chemfp
    from chemfp import search
    
    def get_neighbors_as_arena(arena, id, threshold):
        fp = arena.get_fingerprint_by_id(id)
        neighbor_results =  search.threshold_tanimoto_search_fp(fp, chebi, threshold=threshold)
        neighbor_arena = arena.copy(neighbor_results.get_indices())
        return neighbor_arena
    
    chebi = chemfp.load_fingerprints("chebi_maccs.fps")
    
    # coenzyme A
    coA_arena = get_neighbors_as_arena(chebi, "CHEBI:15346", threshold=0.95)
    print len(coA_arena), "coenzyme A-like structures"
    
    # penicillin
    penicillin_arena = get_neighbors_as_arena(chebi, "CHEBI:17334", threshold=0.9)
    print len(penicillin_arena), "penicillin-like structures"
    
    # I'll compute a profile at different thresholds
    thresholds = [0.25, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    
    # Compare the two sets. (For this case the speed difference between a threshold
    # of 0.25 and 0.0 is not noticible, but having it makes me feel better.)
    coA_against_penicillin_result= search.threshold_tanimoto_search_arena(
        coA_arena, penicillin_arena, threshold=min(thresholds))
    
    # Show a similarity profile
    print "Counts  coA/penicillin"
    for threshold in thresholds:
        print " %.2f      %5d" % (threshold,
                                  coA_against_penicillin_result.count_all(min_score=threshold))

This gives a not very useful output:

    117 coenzyme A-like structures
    15 penicillin-like structures
    Counts  coA/penicillin
     0.25       1755
     0.50        445
     0.60          0
     0.70          0
     0.80          0
     0.90          0
     0.95          0

It's not useful because it's not possible to make any decisions from
this. Are the numbers high or low? It should be low, because these are
two quite different structure classes, but there's nothing to compare
it against.

I need some sort of background reference. What I'll two is construct
two randomly chosen sets, one with 117 fingerprints and the other with
15, and generate the same similarity profile with them. That isn't
quite fair, since randomly chosen sets will most likely be
diverse. Instead, I'll pick one fingerprint at random, then get its
117 or 15, respectively, nearest neighbors as the set members::

    # Get background statistics for random similarity groups of the same size
    import random
    
    # Find a fingerprint at random, get its k neighbors, return them as a new arena
    def get_random_fp_and_its_k_neighbors(arena, k):
        fp = arena[random.randrange(len(arena))][1]
        similar_search = search.knearest_tanimoto_search_fp(fp, arena, k)
        return arena.copy(similar_search.get_indices())

I'll construct 1000 pairs of sets this way, accumulate the threshold
profile, and compare the CoA/penicillin profile to it::

    # Initialize the threshold counts to 0
    total_background_counts = dict.fromkeys(thresholds, 0)
    
    REPEAT = 1000
    for i in range(REPEAT):
        # Select background sets of the same size and accumulate the threshold count totals
        set1 = get_random_fp_and_its_k_neighbors(chebi, len(coA_arena))
        set2 = get_random_fp_and_its_k_neighbors(chebi, len(penicillin_arena))
        background_search = search.threshold_tanimoto_search_arena(set1, set2, threshold=min(thresholds))
        for threshold in thresholds:
            total_background_counts[threshold] += background_search.count_all(min_score=threshold)
    
    print "Counts  coA/penicillin  background"
    for threshold in thresholds:
        print " %.2f      %5d          %5d" % (threshold,
                                               coA_against_penicillin_result.count_all(min_score=threshold),
                                               total_background_counts[threshold] / (REPEAT+0.0))

Your output should look something like::

    Counts  coA/penicillin  background
     0.25       1755            423
     0.50        445             82
     0.60          0             38
     0.70          0             17
     0.80          0              6
     0.90          0              4
     0.95          0              1


This is a bit hard to interpret. Clearly the coenzyme A and penicillin
sets are not closely similar, but for low Tanimoto scores the
similarity is higher than expected.


That difficulty is okay for now because I mostly wanted to show an
example of how to use the chemfp API. If you want to dive deeper into
this sort of analysis, then look into the `Similarity Ensemble Approach`
(SEA) work of Keiser, Roth, Armbruster, Ernsberger, and Irwin. The
paper online and is available from http://sea.bkslab.org/ .

The paper actually wants you to use the `raw score`_. This is the sum
of the hit scores in a given range, and not just the number of
hits. No problem! Use `SearchResult.cumulative_score` for an
individual result or `SearchResults.cumulative_score_all` for the
entire set of results::

    >>> sum(row.cumulative_score(min_score=0.5, max_score=0.9)
    ...             for row in coA_against_penicillin_result)
    224.83239025119906
    >>> coA_against_penicillin_result.cumulative_score_all(min_score=0.5, max_score=0.9)
    224.83239025119866

These also take the `interval` parameter if you don't want the default
of `[]`_.

You may wonder why these two values aren't exactly the same. Addition
of floating point numbers isn't associative. You can see that I get
still different results if I sum up the values in reverse order::

    >>> sum(list(row.cumulative_score(min_score=0.5, max_score=0.9)
    ...                for row in coA_against_penicillin_result)[::-1])
    224.83239025119875

