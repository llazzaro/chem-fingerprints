"Ways to search an arena"

import _chemfp
import ctypes
import array


class SearchResultsHandle(object):
    def __init__(self, size):
        self.size = size
        self.handle = 0
        self.handle = _chemfp.alloc_threshold_results(size)

    def __del__(self, free=_chemfp.free_threshold_results):
        if self.handle:
            free(self.handle, self.size)

    def __int__(self):
        return self.handle
    def __long__(self):
        return self.handle



class SearchResults(object):
    def __init__(self, num_results, target_ids):
        self.num_results = num_results
        if num_results:
            self._handle = SearchResultsHandle(num_results)
        else:
            self._handle = None
        self.target_ids = target_ids

    def __len__(self):
        return self.num_results

    def size(self, i):
        i = xrange(self.num_results)[i]  # Use this trick to support negative index lookups
        return _chemfp.get_num_threshold_hits(self._handle, i)

    def __getitem__(self, i):
        i = xrange(self.num_results)[i]  # Use this trick to support negative index lookups
        ids = self.target_ids
        return [(ids[idx], score) for (idx, score) in
                    _chemfp.threshold_result_get_hits(self._handle, i)]
    
    def __iter__(self): # XXXX fixme
        ids = self.target_ids
        for i in range(0, self.num_results):
            yield [(ids[idx], score) for (idx, score) in
                         _chemfp.threshold_result_get_hits(self._handle, i)]

    def iter_hits(self):
        for i in range(0, self.num_results):
            yield _chemfp.threshold_result_get_hits(self._handle, i)

    def iter_indices(self):
        # This can be optimized with more C code
        for i in range(0, self.num_results):
            yield [idx for (idx, score) in
                        _chemfp.threshold_result_get_hits(self._handle, i)]

    def iter_ids_and_scores(self):
        ids = self.target_ids
        for i in range(0, self.num_results):
            yield [(ids[idx], score) for (idx, score) in
                         _chemfp.threshold_result_get_hits(self._handle, i)]

        
def require_matching_fp_size(query_fp, target_arena):
    if len(query_fp) != target_arena.metadata.num_bytes:
        raise ValueError("query_fp uses %d bytes while target_arena uses %d bytes" % (
            len(query_fp), target_arena.metadata.num_bytes))

def require_matching_sizes(query_arena, target_arena):
    assert query_arena.metadata.num_bits is not None, "arenas must define num_bits"
    assert target_arena.metadata.num_bits is not None, "arenas must define num_bits"
    if query_arena.metadata.num_bits != target_arena.metadata.num_bits:
        raise ValueError("query_arena has %d bits while target_arena has %d bits" % (
            query_arena.metadata.num_bits, target_arena.metadata.num_bits))
    if query_arena.metadata.num_bytes != target_arena.metadata.num_bytes:
        raise ValueError("query_arena uses %d bytes while target_arena uses %d bytes" % (
            query_arena.metadata.num_bytes, target_arena.metadata.num_bytes))
    



def count_tanimoto_hits_fp(query_fp, target_arena, threshold):
    require_matching_fp_size(query_fp, target_arena)
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


def count_tanimoto_hits(query_arena, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)

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

def count_tanimoto_hits_symmetric(arena, threshold):
    N = len(arena)
    counts = (ctypes.c_int * N)()

    _chemfp.count_tanimoto_hits_arena_symmetric(
        threshold, arena.num_bits,
        arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
        0, N, 0, N,
        arena.popcount_indices,
        counts)

    return counts


def partial_count_tanimoto_hits_symmetric(counts, arena, threshold,
                                          query_start=0, query_end=None,
                                          target_start=0, target_end=None):
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

def threshold_tanimoto_search_fp(query_fp, target_arena, threshold):
    require_matching_fp_size(query_fp, target_arena)

    # Improve the alignment so the faster algorithms can be used
    query_start_padding, query_end_padding, query_fp = _chemfp.align_fingerprint(
        query_fp, target_arena.alignment, target_arena.storage_size)


    results = _chemfp.alloc_threshold_results(1)
    try:
        _chemfp.threshold_tanimoto_arena(
            threshold, target_arena.num_bits,
            query_start_padding, query_end_padding, target_arena.storage_size, query_fp, 0, 1,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena,
            target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            results, 0)
        return _chemfp.threshold_result_get_hits(results, 0)
    finally:
        _chemfp.free_threshold_results(results, 1)


def threshold_tanimoto_search(query_arena, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)

    num_queries = len(query_arena)

    results = SearchResults(num_queries, target_arena.ids)
    if num_queries:
        _chemfp.threshold_tanimoto_arena(
            threshold, target_arena.num_bits,
            query_arena.start_padding, query_arena.end_padding,
            query_arena.storage_size, query_arena.arena, query_arena.start, query_arena.end,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            results._handle, 0)
    
    return results

def threshold_tanimoto_search_symmetric(arena, threshold, include_lower_triangle=True):
    assert arena.popcount_indices
    N = len(arena)
    results = SearchResults(N, arena.ids)

    if N:
        _chemfp.threshold_tanimoto_arena_symmetric(
            threshold, arena.num_bits,
            arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
            0, N, 0, N,
            arena.popcount_indices,
            results._handle)

        if include_lower_triangle:
            _chemfp.fill_lower_triangle(results._handle, N)
        
    return results



#def XXXpartial_threshold_tanimoto_search(results, query_arena, target_arena, threshold,
#                                      results_offsets=0):
#    pass

def partial_threshold_tanimoto_search_symmetric(results, arena, threshold,
                                                query_start=0, query_end=None,
                                                target_start=0, target_end=None):
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
            results._handle)


def fill_lower_triangle(results):
    _chemfp.fill_lower_triangle(results._handle, len(results))




# These all return indices into the arena!

def knearest_tanimoto_search_fp(query_fp, target_arena, k, threshold):
    require_matching_fp_size(query_fp, target_arena)
    query_start_padding, query_end_padding, query_fp = _chemfp.align_fingerprint(
        query_fp, target_arena.alignment, target_arena.storage_size)
    
    if k < 0:
        raise ValueError("k must be non-negative")

    results = _chemfp.alloc_threshold_results(1)
    try:
        _chemfp.knearest_tanimoto_arena(
            k, threshold, target_arena.num_bits,
            query_start_padding, query_end_padding, target_arena.storage_size, query_fp, 0, 1,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            results, 0)
        _chemfp.knearest_results_finalize(results, 0, 1)
        return _chemfp.threshold_result_get_hits(results, 0)
    finally:
        _chemfp.free_threshold_results(results, 1)
        pass


def knearest_tanimoto_search(query_arena, target_arena, k, threshold):
    require_matching_sizes(query_arena, target_arena)

    num_queries = len(query_arena)

    results = SearchResults(num_queries, target_arena.ids)

    _chemfp.knearest_tanimoto_arena(
        k, threshold, target_arena.num_bits,
        query_arena.start_padding, query_arena.end_padding,
        query_arena.storage_size, query_arena.arena, query_arena.start, query_arena.end,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        results._handle, 0)
    
    _chemfp.knearest_results_finalize(results._handle, 0, num_queries)
    
    return results


def knearest_tanimoto_search_symmetric(arena, k, threshold):
    N = len(arena)

    results = SearchResults(N, arena.ids)

    _chemfp.knearest_tanimoto_arena_symmetric(
        k, threshold, arena.num_bits,
        arena.start_padding, arena.end_padding, arena.storage_size, arena.arena,
        0, N, 0, N,
        arena.popcount_indices,
        results._handle)
    _chemfp.knearest_results_finalize(results._handle, 0, N)
    
    return results



def partial_knearest_tanimoto_search(results, arena, k, threshold, results_offset=0):
    pass

def partial_knearest_tanimoto_search_symmetric(results, arena, k, threshold,
                                               query_start=0, query_end=None,
                                               target_start=0, target_end=None):
    pass
