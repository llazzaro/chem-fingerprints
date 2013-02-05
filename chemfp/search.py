"""Search a FingerprintArena and work with the search results

This module implements the different ways to search a
`FingerprintArena`. The search functions are:

  Count the number of hits:
    count_tanimoto_hits_fp - search an arena using a single fingerprint
    count_tanimoto_hits_arena - search an arena using an arena
    count_tanimoto_hits_symmetric - search an arena using itself
    partial_count_tanimoto_hits_symmetric - (advanced use; see the doc string)

  Find all hits at or above a given threshold, sorted arbitrarily:
    threshold_tanimoto_search_fp - search an arena using a single fingerprint
    threshold_tanimoto_search_arena - search an arena using an arena
    threshold_tanimoto_search_symmetric - search an arena using itself
    partial_threshold_tanimoto_search_symmetric - (advanced use; see the doc string)
    fill_lower_triangle - copy the upper triangle terms to the lower triangle

  Find the k-nearest hits at or above a given threshold, sorted by
  decreasing similarity:
    knearest_tanimoto_search_fp - search an arena using a single fingerprint
    knearest_tanimoto_search_arena - search an arena using an arena
    knearest_tanimoto_search_symmetric - search an arena using itself

The threshold and k-nearest search results use a `SearchResult` when
a fingerprint is used as a query, or a `SearchResults` when an arena
is used as a query. These internally use a compressed sparse row format.
"""

import _chemfp
import ctypes
import array

# 
__all__ = ["SearchResult", "SearchResults",
           "count_tanimoto_hits_fp", "count_tanimoto_hits_arena",
           "count_tanimoto_hits_symmetric", "partial_count_tanimoto_hits_symmetric",

           "threshold_tanimoto_search_fp", "threshold_tanimoto_search_arena",
           "threshold_tanimoto_search_symmetric", "partial_threshold_tanimoto_search_symmetric",
           "fill_lower_triangle",

           "knearest_tanimoto_search_fp", "knearest_tanimoto_search_arena",
           "knearest_tanimoto_search_symmetric"
           ]
           

class SearchResult(object):
    """Search results for a query fingerprint against a target arena.

    The results contains a list of hits. Hits contain a target index,
    score, and optional target ids. The hits can be reordered based on
    score or index.
    
    """
    def __init__(self, search_results, row):
        self._search_results = search_results
        self._row = row

    def __len__(self):
        """The number of hits"""
        return self._search_results._size(self._row)

    def __iter__(self):
        """Iterate through the pairs of (target index, score) using the current ordering"""
        return iter(self._search_results._get_indices_and_scores(self._row))
        
    def clear(self):
        """Remove all hits from this result"""
        self._search_results._clear_row(self._row)

    def get_indices(self):
        """The list of target indices, in the current ordering."""
        return self._search_results._get_indices(self._row)

    def get_ids(self):
        """The list of target identifiers (if available), in the current ordering"""
        ids = self._search_results.target_ids
        if ids is None:
            return None
        return [ids[i] for i in self._search_results._get_indices(self._row)]
    
    def get_scores(self):
        """The list of target scores, in the current ordering"""
        return self._search_results._get_scores(self._row)
        
    def get_ids_and_scores(self):
        """The list of (target identifier, target score) pairs, in the current ordering

        Raises a TypeError if the target IDs are not available.
        """
        ids = self._search_results.target_ids
        if ids is None:
            raise TypeError("target_ids are not available")
        return zip(self.get_ids(), self.get_scores())

    def get_indices_and_scores(self):
        """The list of (target index, score) pairs, in the current ordering"""
        return self._search_results._get_indices_and_scores(self._row)
            
    def reorder(self, ordering="decreasing-score"):
        """Reorder the hits based on the requested ordering.

        The available orderings are:
          increasing-score: sort by increasing score
          decreasing-score: sort by decreasing score
          increasing-index: sort by increasing target index
          decreasing-index: sort by decreasing target index
          move-closest-first: move the hit with the highest score to the first position
          reverse: reverse the current ordering

        :param ordering: the name of the ordering to use
        """
        self._search_results._reorder_row(self._row, ordering)

    def count(self, min_score=None, max_score=None, interval="[]"):
        """Count the number of hits with a score between `min_score` and `max_score`

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
        """
        return self._search_results._count_row(self._row, min_score, max_score, interval)

    def cumulative_score(self, min_score=None, max_score=None, interval="[]"):
        """The sum of the scores which are between `min_score` and `max_score`

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
        """
        return self._search_results._cumulative_score_row(self._row, min_score, max_score, interval)

    ## ??? What does this do?
    ## @property
    ## def target_id(self):
    ##     ids = self._search_results.target_ids
    ##     if ids is None:
    ##         return None
    ##     return ids[self._row]

class SearchResults(_chemfp.SearchResults):
    """Search results for a list of query fingerprints against a target arena

    This acts like a list of SearchResult elements, with the ability
    to iterate over each search results, look them up by index, and
    get the number of scores.

    In addition, there are helper methods to iterate over each hit and
    to get the hit indicies, scores, and identifiers directly as Python
    lists, sort the list contents, and more.
    
    """
    def __init__(self, n, arena_ids=None):
        """`n` is the number of SearchResult instances and `arena_ids` the target arena ids

        There is one SearchResult for each query fingerprint. The `arena_ids`
        are used to map the hit index back to the hit id.
        """
        super(SearchResults, self).__init__(n, arena_ids)
        self._results = [SearchResult(self, i) for i in xrange(n)]

    def __iter__(self):
        """Iterate over each SearchResult hit"""
        return iter(self._results)

    def __len__(self):
        """The number of rows in the SearchResults"""
        super(SearchResults, self).__len__(self)

    def __getitem__(self, i):
        """Get the `i`th SearchResult"""
        try:
            return self._results[i]
        except IndexError:
            raise IndexError("row index is out of range")

    def iter_indices(self):
        """For each hit, yield the list of target indices"""
        for i in xrange(len(self)):
            yield self._get_indices(i)

    def iter_ids(self):
        """For each hit, yield the list of target identifiers"""
        ids = self.target_ids
        for indicies in self.iter_indices():
            yield [ids[idx] for idx in indicies]

    def iter_scores(self):
        """For each hit, yield the list of target scores"""
        for i in xrange(len(self)):
            yield self._get_scores(i)

    def iter_indices_and_scores(self):
        """For each hit, yield the list of (target index, score) tuples"""
        for i in xrange(len(self)):
            yield zip(self._get_indices(i), self._get_scores(i))
    
    def iter_ids_and_scores(self):
        """For each hit, yield the list of (target id, score) tuples"""
        ids = self.target_ids
        for i in xrange(len(self)):
            yield [(ids[idx], score) for (idx, score) in self[i]]


    # I don't like how C-level doc strings can't report the call
    # signature even though keyword arguments are supported. I also
    # don't like maintaining the docstrings in C code.
    # Problem solved by interposing these Python methods
    def clear_all(self):
        """Remove all hits from all of the search results"""
        return super(SearchResults, self).clear_all()
            
    def count_all(self, min_score=None, max_score=None, interval="[]"):
        """Count the number of hits with a score between `min_score` and `max_score`
        
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
        """
        return super(SearchResults, self).count_all(min_score, max_score, interval)
        
    def cumulative_score_all(self, min_score=None, max_score=None, interval="[]"):
        """The sum of all scores in all rows which are between `min_score` and `max_score`

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
        """
        return super(SearchResults, self).cumulative_score_all(min_score, max_score, interval)

    def reorder_all(self, order="decreasing-score"):
        """Reorder the hits for all of the rows based on the requested `order`.
        
        The available orderings are:
          increasing-score: sort by increasing score
          decreasing-score: sort by decreasing score
          increasing-index: sort by increasing target index
          decreasing-index: sort by decreasing target index
          move-closest-first: move the hit with the highest score to the first position
          reverse: reverse the current ordering
        
        :param ordering: the name of the ordering to use
        """
        return super(SearchResults, self).reorder_all(order)


        
def _require_matching_fp_size(query_fp, target_arena):
    if len(query_fp) != target_arena.metadata.num_bytes:
        raise ValueError("query_fp uses %d bytes while target_arena uses %d bytes" % (
            len(query_fp), target_arena.metadata.num_bytes))

def _require_matching_sizes(query_arena, target_arena):
    assert query_arena.metadata.num_bits is not None, "arenas must define num_bits"
    assert target_arena.metadata.num_bits is not None, "arenas must define num_bits"
    if query_arena.metadata.num_bits != target_arena.metadata.num_bits:
        raise ValueError("query_arena has %d bits while target_arena has %d bits" % (
            query_arena.metadata.num_bits, target_arena.metadata.num_bits))
    if query_arena.metadata.num_bytes != target_arena.metadata.num_bytes:
        raise ValueError("query_arena uses %d bytes while target_arena uses %d bytes" % (
            query_arena.metadata.num_bytes, target_arena.metadata.num_bytes))
    



def count_tanimoto_hits_fp(query_fp, target_arena, threshold=0.7):
    """Count the number of hits in `target_arena` at least `threshold` similar to the `query_fp`

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
    """
    _require_matching_fp_size(query_fp, target_arena)
    # Improve the alignment so the faster algorithms can be used
    query_start_padding, query_end_padding, query_fp = _chemfp.align_fingerprint(
        query_fp, target_arena.alignment, target_arena.storage_size)
                                                 
    counts = array.array("i", (0 for i in xrange(len(query_fp))))
    _chemfp.count_tanimoto_arena(threshold, target_arena.num_bits,
                                 query_start_padding, query_end_padding,
                                 target_arena.storage_size, query_fp, 0, 1,
                                 target_arena.start_padding, target_arena.end_padding,
                                 target_arena.storage_size, target_arena.arena,
                                 target_arena.start, target_arena.end,
                                 target_arena.popcount_indices,
                                 counts)
    return counts[0]


def count_tanimoto_hits_arena(query_arena, target_arena, threshold=0.7):
    """For each fingerprint in `query_arena`, count the number of hits in `target_arena` at least `threshold` similar to it
    
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
    """
    _require_matching_sizes(query_arena, target_arena)

    counts = (ctypes.c_int*len(query_arena))()
    _chemfp.count_tanimoto_arena(threshold, target_arena.num_bits,
                                 query_arena.start_padding, query_arena.end_padding,
                                 query_arena.storage_size,
                                 query_arena.arena, query_arena.start, query_arena.end,
                                 target_arena.start_padding, target_arena.end_padding,
                                 target_arena.storage_size,
                                 target_arena.arena, target_arena.start, target_arena.end,
                                 target_arena.popcount_indices,
                                 counts)
    return counts    

def count_tanimoto_hits_symmetric(arena, threshold=0.7, batch_size=100):
    """For each fingerprint in the `arena`, count the number of other fingerprints at least `threshold` similar to it

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
    """
    N = len(arena)
    counts = (ctypes.c_int * N)()

    # This spends the entire time in C, which means ^C won't work until it finishes.
    # While it's theoretically slightly higher performance, I can't measure the
    # difference, and it's much better to let people be able to interrupt the program.
    #    _chemfp.count_tanimoto_hits_arena_symmetric(
    #        threshold, arena.num_bits,
    #        arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
    #        0, N, 0, N,
    #        arena.popcount_indices,
    #        counts)
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    # Process N rows at a time, which lets Python handle ^C at times.
    # Since the code processes a triangle, this means that early
    # on there will be more time between ^C checks than later.
    # I'm not able to detect the Python overhead, so I'm not going
    # to make it more "efficient".
    for query_start in xrange(0, N, batch_size):
        query_end = min(query_start + batch_size, N)
        _chemfp.count_tanimoto_hits_arena_symmetric(
            threshold, arena.num_bits,
            arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
            query_start, query_end, 0, N,
            arena.popcount_indices,
            counts)

    return counts


def partial_count_tanimoto_hits_symmetric(counts, arena, threshold=0.7,
                                          query_start=0, query_end=None,
                                          target_start=0, target_end=None):
    """Compute a portion of the symmetric Tanimoto counts

    For most cases, use count_tanimoto_hits_symmetric instead of this
    function!
    
    This function is only useful for thread-pool implementations. In
    that case, set the number of OpenMP threads to 1.

    `counts` is a contiguous array of integers. It should be
    initialized to zeros, and reused for successive calls.

    The function adds counts for counts[query_start:query_end] based
    on computing the upper-triangle portion contained in the rectangle
    query_start:query_end and target_start:target_end and using
    symmetry to fill in the lower half.

    You know, this is pretty complicated. Here's the bare minimum
    example of how to use it correctly to process 10 rows at a time
    using up to 4 threads::

        import chemfp
        import chemfp.search
        from chemfp import futures
        import array
        
        chemfp.set_num_threads(1)  # Globally disable OpenMP
        
        arena = chemfp.load_fingerprints("targets.fps")  # Load the fingerprints
        n = len(arena)
        counts = array.array("i", [0]*n)
        
        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            for row in xrange(0, n, 10):
                executor.submit(chemfp.search.partial_count_tanimoto_hits_symmetric,
                                counts, arena, threshold=0.2,
                                query_start=row, query_end=min(row+10, n))
        
        print counts

    :param counts: the accumulated Tanimoto counts
    :type counts: a contiguous block of integer
    :param arena: the fingerprints.
    :type arena: a FingerprintArena
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param query_start: the query start row
    :type query_start: an integer
    :param query_end: the query end row
    :type query_end: an integer, or None to mean the last query row
    :param target_start: the target start row
    :type target_start: an integer
    :param target_end: the target end row
    :type target_end: an integer, or None to mean the last target row
    :returns: nothing
    """
    N = len(arena)
    
    if query_end is None:
        query_end = N
    elif query_end > N:
        query_end = N
        
    if target_end is None:
        target_end = N
    elif target_end > N:
        target_end = N

    if query_end > len(counts):
        raise ValueError("counts array is too small for the given query range")
    if target_end > len(counts):
        raise ValueError("counts array is too small for the given target range")

    _chemfp.count_tanimoto_hits_arena_symmetric(
        threshold, arena.num_bits,
        arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
        query_start, query_end, target_start, target_end,
        arena.popcount_indices,
        counts)


# These all return indices into the arena!

def threshold_tanimoto_search_fp(query_fp, target_arena, threshold=0.7):
    """Search for fingerprint hits in `target_arena` which are at least `threshold` similar to `query_fp`

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
    """
    _require_matching_fp_size(query_fp, target_arena)

    # Improve the alignment so the faster algorithms can be used
    query_start_padding, query_end_padding, query_fp = _chemfp.align_fingerprint(
        query_fp, target_arena.alignment, target_arena.storage_size)


    results = SearchResults(1, target_arena.arena_ids)
    _chemfp.threshold_tanimoto_arena(
        threshold, target_arena.num_bits,
        query_start_padding, query_end_padding, target_arena.storage_size, query_fp, 0, 1,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena,
        target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        results, 0)
    return results[0]


def threshold_tanimoto_search_arena(query_arena, target_arena, threshold=0.7):
    """Search for the hits in the `target_arena` at least `threshold` similar to the fingerprints in `query_arena`

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
    """
    _require_matching_sizes(query_arena, target_arena)

    num_queries = len(query_arena)

    results = SearchResults(num_queries, target_arena.arena_ids)
    if num_queries:
        _chemfp.threshold_tanimoto_arena(
            threshold, target_arena.num_bits,
            query_arena.start_padding, query_arena.end_padding,
            query_arena.storage_size, query_arena.arena, query_arena.start, query_arena.end,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            results, 0)
    
    return results

def threshold_tanimoto_search_symmetric(arena, threshold=0.7, include_lower_triangle=True, batch_size=100):
    """Search for the hits in the `arena` at least `threshold` similar to the fingerprints in the arena

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
    """
    
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    N = len(arena)
    results = SearchResults(N, arena.arena_ids)

    if N:
        # Break it up into batch_size groups in order to let Python's
        # interrupt handler check for a ^C, which is otherwise
        # suppressed until the function finishes.
        for query_start in xrange(0, N, batch_size):
            query_end = min(query_start + batch_size, N)
            _chemfp.threshold_tanimoto_arena_symmetric(
                threshold, arena.num_bits,
                arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
                query_start, query_end, 0, N,
                arena.popcount_indices,
                results)

        if include_lower_triangle:
            _chemfp.fill_lower_triangle(results, N)
        
    return results



#def XXXpartial_threshold_tanimoto_search(results, query_arena, target_arena, threshold,
#                                      results_offsets=0):
#    pass

def partial_threshold_tanimoto_search_symmetric(results, arena, threshold=0.7,
                                                query_start=0, query_end=None,
                                                target_start=0, target_end=None):
    """Compute a portion of the symmetric Tanimoto search results

    For most cases, use threshold_tanimoto_arena_symmetric instead of this
    function!
    
    This function is only useful for thread-pool implementations. In
    that case, set the number of OpenMP threads to 1.

    `results` is a SearchResults instance which is at least as large
    as the arena. It should be reused for successive updates.

    The function adds hits to results[query_start:query_end] based
    on computing the upper-triangle portion contained in the rectangle
    query_start:query_end and target_start:target_end.

    It does not fill in the lower triangle. To get the full matrix,
    call `fill_lower_triangle`.

    You know, this is pretty complicated. Here's the bare minimum
    example of how to use it correctly to process 10 rows at a time
    using up to 4 threads::

        import chemfp
        import chemfp.search
        from chemfp import futures
        import array

        chemfp.set_num_threads(1)

        arena = chemfp.load_fingerprints("targets.fps")
        n = len(arena)
        results = chemfp.search.SearchResults(n, arena.ids)

        with futures.ThreadPoolExecutor(max_workers=4) as executor:
            for row in xrange(0, n, 10):
                executor.submit(chemfp.search.partial_threshold_tanimoto_search_symmetric,
                                results, arena, threshold=0.2,
                                query_start=row, query_end=min(row+10, n))

        chemfp.search.fill_lower_triangle(results)

    The hits in the `SearchResults` are in arbitrary order.

    :param counts: the intermediate search results
    :type counts: a SearchResults instance
    :param arena: the fingerprints.
    :type arena: a FingerprintArena
    :param threshold: The minimum score threshold.
    :type threshold: float between 0.0 and 1.0, inclusive
    :param query_start: the query start row
    :type query_start: an integer
    :param query_end: the query end row
    :type query_end: an integer, or None to mean the last query row
    :param target_start: the target start row
    :type target_start: an integer
    :param target_end: the target end row
    :type target_end: an integer, or None to mean the last target row
    :returns: nothing
    """
    assert arena.popcount_indices
    N = len(arena)
    
    if query_end is None:
        query_end = N
    elif query_end > N:
        query_end = N
        
    if target_end is None:
        target_end = N
    elif target_end > N:
        target_end = N

    if query_end > N:
        raise ValueError("counts array is too small for the given query range")
    if target_end > N:
        raise ValueError("counts array is too small for the given target range")

    if N:
        _chemfp.threshold_tanimoto_arena_symmetric(
            threshold, arena.num_bits,
            arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
            query_start, query_end, target_start, target_end,
            arena.popcount_indices,
            results)


def fill_lower_triangle(results):
    """Duplicate each entry of `results` to its transpose

    This is used after the symmetric threshold search to turn the
    upper-triangle results into a full matrix.
    """
    _chemfp.fill_lower_triangle(results, len(results))



# These all return indices into the arena!

def knearest_tanimoto_search_fp(query_fp, target_arena, k=3, threshold=0.7):
    """Search for `k`-nearest hits in `target_arena` which are at least `threshold` similar to `query_fp`

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
    """
    _require_matching_fp_size(query_fp, target_arena)
    query_start_padding, query_end_padding, query_fp = _chemfp.align_fingerprint(
        query_fp, target_arena.alignment, target_arena.storage_size)
    
    if k < 0:
        raise ValueError("k must be non-negative")

    results = SearchResults(1, target_arena.arena_ids)
    _chemfp.knearest_tanimoto_arena(
        k, threshold, target_arena.num_bits,
        query_start_padding, query_end_padding, target_arena.storage_size, query_fp, 0, 1,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        results, 0)
    _chemfp.knearest_results_finalize(results, 0, 1)

    return results[0]

def knearest_tanimoto_search_arena(query_arena, target_arena, k=3, threshold=0.7):
    """Search for the `k` nearest hits in the `target_arena` at least `threshold` similar to the fingerprints in `query_arena`

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
    """
    _require_matching_sizes(query_arena, target_arena)

    num_queries = len(query_arena)

    results = SearchResults(num_queries, target_arena.arena_ids)

    _chemfp.knearest_tanimoto_arena(
        k, threshold, target_arena.num_bits,
        query_arena.start_padding, query_arena.end_padding,
        query_arena.storage_size, query_arena.arena, query_arena.start, query_arena.end,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        results, 0)
    
    _chemfp.knearest_results_finalize(results, 0, num_queries)
    
    return results


def knearest_tanimoto_search_symmetric(arena, k=3, threshold=0.7, batch_size=100):
    """Search for the `k`-nearest hits in the `arena` at least `threshold` similar to the fingerprints in the arena

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
    """
    N = len(arena)
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

    results = SearchResults(N, arena.arena_ids)

    if N:
        # Break it up into batch_size groups in order to let Python's
        # interrupt handler check for a ^C, which is otherwise
        # suppressed until the function finishes.
        for query_start in xrange(0, N, batch_size):
            query_end = min(query_start + batch_size, N)            
            _chemfp.knearest_tanimoto_arena_symmetric(
                k, threshold, arena.num_bits,
                arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
                query_start, query_end, 0, N,
                arena.popcount_indices,
                results)
        _chemfp.knearest_results_finalize(results, 0, N)
    
    return results
