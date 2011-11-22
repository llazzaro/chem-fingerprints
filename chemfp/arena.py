"""Algorithms and data structure for working with a FingerprintArena.

NOTE: This module should not be used directly.

A FingerprintArena stores the fingerprints as a contiguous byte
string, called the `arena`. Each fingerprint takes `storage_size`
bytes, which may be larger than `num_bytes` if the fingerprints have a
specific memory alignment. The bytes for fingerprint i are
  arena[i*storage_size:i*storage_size+num_bytes]
Additional bytes must contain NUL bytes.

The lookup for `ids[i]` contains the id for fingerprint `i`.

A FingerprintArena has an optional `indices` attribute. When
available, it means that the arena fingerprints and corresponding ids
are ordered by population count, and the fingerprints with popcount
`p` start at index indices[p] and end just before indices[p+1].

"""

from __future__ import absolute_import

import ctypes
from cStringIO import StringIO
import array

from chemfp import FingerprintReader, check_fp_problems, check_metadata_problems
import _chemfp
from chemfp import bitops

__all__ = []

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
    counts = array.array("i", (0 for i in xrange(len(query_fp))))
    _chemfp.count_tanimoto_arena(threshold, target_arena.num_bits,
                                 0, 0, len(query_fp), query_fp, 0, -1,
                                 target_arena.start_padding, target_arena.end_padding,
                                 target_arena.storage_size, target_arena.arena,
                                 target_arena.start, target_arena.end,
                                 target_arena.popcount_indices,
                                 counts)
    return counts[0]


def count_tanimoto_hits_arena(query_arena, target_arena, threshold):
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


# Search results stored in a compressed sparse row form

class SearchResults(object):
    """Contains the result of a Tanimoto threshold or k-nearest search

    Each result contains a list of hits, where the hit is a
    two-element tuple. If you iterate over the SearchResult then
    you'll get the hits as (target_id, target_score) pairs.
    tuples. If you iterate using the method `iter_hits()` then you'll
    get the hits as (target_index, target_score) pairs.

    """
    def __init__(self, offsets, indices, scores, query_ids, target_ids):
        assert len(offsets) > 0
        self.offsets = offsets
        self.indices = indices
        self.scores = scores
        self.query_ids = query_ids
        self.target_ids = target_ids
        
    def __len__(self):
        """Number of search results"""
        return len(self.offsets)-1

    def size(self, i):
        """The number of hits for result at position i

        :param i: index into the search results
        :type i: int
        :returns: int
        
        """
        i = xrange(len(self.offsets)-1)[i]  # Use this trick to support negative index lookups
        return self.offsets[i+1]-self.offsets[i]
    
    def __getitem__(self, i):
        """The list of hits for result at position i

        Each hit contains a (id, score) tuple.
        """
        i = xrange(len(self.offsets)-1)[i]  # Use this trick to support negative index lookups
        start, end = self.offsets[i:i+2]
        ids = self.target_ids
        return zip((ids[idx] for idx in self.indices[start:end]), self.scores[start:end])

    def __iter__(self):
        """Iterate over the named hits for each result

        Each term is a list of hits. A hit contains (id, score) tuples.
        The order of the hits depends on the search algorithm.
        """
        target_ids = self.target_ids
        indices = self.indices
        scores = self.scores
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield zip((target_ids[index] for index in indices[start:end]),
                      scores[start:end])
            start = end

    def iter_hits(self):
        """Iterate over the indexed hits for each result

        Each term is a list of hits. A hit contains (index, score) tuples.
        The order of the hits depends on the search algorithm.
        """
        indices = self.indices
        scores = self.scores
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield zip(indices[start:end], scores[start:end])
            start = end

    def iter_indices(self):
        """Iterate over the indices for each result

        Each term is a list of indices. The order of the indices in a list
        depends on the search algorithm.
        """
        indices = self.indices
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield indices[start:end]
            start = end
        

def threshold_tanimoto_search_fp_indices(query_fp, target_arena, threshold):
    require_matching_fp_size(query_fp, target_arena)

    offsets = (ctypes.c_int * 2)()
    offsets[0] = 0
    
    num_cells = len(target_arena)
    indices = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    num_added = _chemfp.threshold_tanimoto_arena(
        threshold, target_arena.num_bits,
        0, 0, len(query_fp), query_fp, 0, -1,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena,
        target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        offsets, 0,
        indices, scores)

    assert num_added == 1

    end = offsets[1]
    return [(indices[i], scores[i]) for i in xrange(end)]

def threshold_tanimoto_search_fp(query_fp, target_arena, threshold):
    require_matching_fp_size(query_fp, target_arena)
    result = threshold_tanimoto_search_fp_indices(query_fp, target_arena, threshold)
    return [(target_arena.ids[index], score) for (index, score) in result]


def threshold_tanimoto_search_arena(query_arena, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)
    num_bits = target_arena.num_bits

    num_queries = len(query_arena)

    offsets = (ctypes.c_int * (num_queries+1))()

    offsets[0] = 0
    
    product = num_queries*len(target_arena)
    if product < 100:
        min_rows = num_queries
    else:
        max_cells = min(10000, product // 4)
        min_rows = max(2, max_cells // len(target_arena))

    num_cells = min_rows * len(target_arena)
    indices = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()
    
    query_start = query_arena.start
    query_end = query_arena.end


    def add_rows(query_start, offset_start):
        return _chemfp.threshold_tanimoto_arena(
            threshold, num_bits,
            query_arena.start_padding, query_arena.end_padding,
            query_arena.storage_size, query_arena.arena, query_start, query_end,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            offsets, offset_start, # XXX should query_start=0?
            indices, scores)

    return _search(query_start, query_end, offsets, indices, scores, add_rows,
                   query_arena.ids, target_arena.ids)

def knearest_tanimoto_search_fp(query_fp, target_arena, k, threshold):
    result = knearest_tanimoto_search_fp_indices(query_fp, target_arena, k, threshold)
    return [(target_arena.ids[index], score) for (index, score) in result]

def knearest_tanimoto_search_fp_indices(query_fp, target_arena, k, threshold):
    require_matching_fp_size(query_fp, target_arena)
    if k < 0:
        raise ValueError("k must be non-negative")

    offsets = (ctypes.c_int * 2)()
    offsets[0] = 0
    indices = (ctypes.c_int * k)()
    scores = (ctypes.c_double * k)()

    num_added = _chemfp.knearest_tanimoto_arena(
        k, threshold, target_arena.num_bits,
        0, 0, len(query_fp), query_fp, 0, 1,
        target_arena.start_padding, target_arena.end_padding,
        target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
        target_arena.popcount_indices,
        offsets, 0,
        indices, scores)
    assert num_added > 0, num_added
    end = offsets[1]
    return [(indices[i], scores[i]) for i in xrange(end)]

def knearest_tanimoto_search_arena(query_arena, target_arena, k, threshold):
    require_matching_sizes(query_arena, target_arena)
    num_bits = query_arena.metadata.num_bits

    num_queries = len(query_arena)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(query_arena))*k

    indices = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    query_start = query_arena.start
    query_end = query_arena.end

    def add_rows(query_start, offset_start):
        return _chemfp.knearest_tanimoto_arena(
            k, threshold, num_bits,
            query_arena.start_padding, query_arena.end_padding,
            query_arena.storage_size, query_arena.arena, query_start, query_end,
            target_arena.start_padding, target_arena.end_padding,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indices,
            offsets, offset_start,
            indices, scores)

    return _search(query_start, query_end, offsets, indices, scores, add_rows,
                   query_arena.ids, target_arena.ids)


# Core of the Tanimoto search routine

def _search(query_start, query_end, offsets, indices, scores,
            add_rows, query_ids, target_ids):
    num_added = add_rows(query_start, 0)
    if num_added == query_end:
        return SearchResults(offsets, indices, scores, query_ids, target_ids)

    query_start = query_start + num_added
    offset_start = num_added

    last = offsets[num_added]
    all_indices = indices[:last]
    all_scores = scores[:last]

    while query_start < query_end:
        num_added = add_rows(query_start, offset_start)
        assert num_added > 0

        prev_last = offsets[query_start]
        all_indices[prev_last:] = indices
        all_scores[prev_last:] = scores

        offset_start += num_added
        query_start += num_added

    return SearchResults(offsets, all_indices, all_scores, query_ids, target_ids)




class FingerprintLookup(object):
    "This is an unpublished API and may be removed in the future"
    def __init__(self, fp_size, start_padding, end_padding, storage_size, arena):
        self._fp_size = fp_size
        self._start_padding = start_padding
        self._end_padding = end_padding
        self._storage_size = storage_size
        self._arena = arena
        self._range_check = xrange(len(self))

    def __len__(self):
        if not self._storage_size:
            return 0
        return (len(self._arena)-self._start_padding-self._end_padding)/self._storage_size

    def __iter__(self):
        fp_size = self._fp_size
        arena = self._arena
        for id, start_offset in zip(self.ids,
                                    xrange(self._start_padding,
                                           len(arena)-self.end_padding, storage_size)):
            yield id, arena[start_offset:start_offset+target_fp_size]
        
        
    def __getitem__(self, i):
        try:
            start_offset = self._range_check[i] * self._storage_size + self._start_padding
        except IndexError:
            raise IndexError("arena fingerprint index out of range")
            
        return self._arena[start_offset:start_offset+self._fp_size]

class FingerprintArena(FingerprintReader):
    """Stores fingerprints in a contiguous block of memory

    The public attributes are:
       metadata
           `Metadata` about the fingerprints
       ids
           list of identifiers, ordered by position
    """
    def __init__(self, metadata,
                 start_padding, end_padding, storage_size, arena,
                 popcount_indices, ids, start=0, end=None):
        if metadata.num_bits is None:
            raise TypeError("Missing metadata num_bits information")
        if metadata.num_bytes is None:
            raise TypeError("Missing metadata num_bytes information")
        self.metadata = metadata
        self.num_bits = metadata.num_bits
        self.start_padding = start_padding
        self.end_padding = end_padding
        self.storage_size = storage_size
        self.arena = arena
        self.popcount_indices = popcount_indices
        self.ids = ids
        self.fingerprints = FingerprintLookup(metadata.num_bytes, self.start_padding,
                                              self.end_padding, storage_size, arena)
        self.start = start
        if end is None:
            if self.metadata.num_bytes:
                end = (len(arena) - start_padding - end_padding) // self.storage_size
            else:
                end = 0
        self.end = end
        assert end >= start
        self._range_check = xrange(end-start)

    def __len__(self):
        """Number of fingerprint records in the FingerprintArena"""
        return self.end - self.start

    def __getitem__(self, i):
        """Return the (id, fingerprint) at position i"""
        if isinstance(i, slice):
            start, end, step = i.indices(self.end - self.start)
            if start >= end:
                return FingerprintArena(self.metadata, 0, 0, self.storage_size, "",
                                        "", [], 0, 0)
            if step != 1:
                raise IndexError("arena slice step size must be 1")
            return FingerprintArena(self.metadata, self.start_padding, self.end_padding,
                                    self.storage_size, self.arena,
                                    self.popcount_indices, self.ids[start:end],
                                    self.start+start, self.start+end)
        try:
            i = self._range_check[i]
        except IndexError:
            raise IndexError("arena fingerprint index out of range")
        arena_i = i + self.start
        start_offset = arena_i * self.storage_size + self.start_padding
        end_offset = start_offset + self.metadata.num_bytes
        return self.ids[i], self.arena[start_offset:end_offset]


    def save(self, destination):
        """Save the arena contents to the given filename or file object"""
        from . import io
        need_close = False
        if isinstance(destination, basestring):
            need_close = True
            output = io.open_output(destination)
        else:
            output = destination

        try:
            io.write_fps1_magic(output)
            io.write_fps1_header(output, self.metadata)
            for id, fp in self:
                io.write_fps1_fingerprint(output, fp, id)
        finally:
            if need_close:
                output.close()
                
    def reset(self):
        """This method is not documented"""
        pass

    def __iter__(self):
        """Iterate over the (id, fingerprint) contents of the arena"""
        storage_size = self.storage_size
        if not storage_size:
            return
        target_fp_size = self.metadata.num_bytes
        arena = self.arena
        for id, start_offset in zip(self.ids,
                                    xrange(self.start*storage_size+self.start_padding,
                                           self.end*storage_size+self.start_padding,
                                           storage_size)):
            yield id, arena[start_offset:start_offset+target_fp_size]

    def iter_arenas(self, arena_size = 1000):
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
            yield self
            return
        
        storage_size = self.storage_size
        start = self.start
        for i in xrange(0, len(self), arena_size):
            ids = self.ids[i:i+arena_size]
            end = start+len(ids)
            yield FingerprintArena(self.metadata, self.start_padding, self.end_padding,
                                   self.storage_size, self.arena,
                                   self.popcount_indices, ids, start, end)
            start = end

    def count_tanimoto_hits_fp(self, query_fp, threshold=0.7):
        """Count the fingerprints which are similar enough to the query fingerprint

        Return the number of fingerprints in this arena which are
        at least `threshold` similar to the query fingerprint `query_fp`.

        :param query_fp: query fingerprint
        :type query_fp: byte string
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: integer count
        """
        return count_tanimoto_hits_fp(query_fp, self, threshold)

    def count_tanimoto_hits_arena(self, query_arena, threshold=0.7):
        """Count the fingerprints which are similar enough to each query fingerprint

        For each fingerprint in the `query_arena`, count the number of
        fingerprints in this arena with Tanimoto similarity of at
        least `threshold`. The resulting list order is the same as the
        query fingerprint order.
        
        :param query_fp: query arena
        :type query_fp: FingerprintArena
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: list of integer counts
        """
        return count_tanimoto_hits_arena(query_arena, self, threshold)

    def threshold_tanimoto_search_fp(self, query_fp, threshold=0.7):
        """Find the fingerprints which are similar enough to the query fingerprint

        Find all of the fingerprints in this arena which are at least
        `threshold` similar to the query fingerprint `query_fp`.
        The hits are returned as a list containing (id, score) tuples
        in arbitrary order.
        
        :param query_fp: query fingerprint
        :type query_fp: byte string
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: list of (int, score) tuples
        """
        return threshold_tanimoto_search_fp(query_fp, self, threshold)

    def threshold_tanimoto_search_arena(self, query_arena, threshold=0.7):
        """Find the fingerprints which are similar to each of the query fingerprints

        For each fingerprint in the `query_arena`, find all of the
        fingerprints in this arena which are at least `threshold`
        similar. The hits are returned as a `SearchResults` instance.
        
        :param query_arena: query arena
        :type query_arena: FingerprintArena
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: SearchResults
        """
        return threshold_tanimoto_search_arena(query_arena, self, threshold)

    def knearest_tanimoto_search_fp(self, query_fp, k=3, threshold=0.7):
        """Find the k-nearest fingerprints which are similar to the query fingerprint

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
        """
        return knearest_tanimoto_search_fp(query_fp, self, k, threshold)

    def knearest_tanimoto_search_arena(self, query_arena, k=3, threshold=0.7):
        """Find the k-nearest fingerprint which are similar to each of the query fingerprints

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
        """
        return knearest_tanimoto_search_arena(query_arena, self, k, threshold)


class ChemFPOrderedPopcount(ctypes.Structure):
    _fields_ = [("popcount", ctypes.c_int),
                ("index", ctypes.c_int)]


_methods = bitops.get_methods()
_has_popcnt = "POPCNT" in _methods
_has_ssse3 = "ssse3" in _methods

def get_optimal_alignment(num_bits):
    if num_bits <= 32:
        # Just in case!
        if num_bits <= 8:
            return 1
        return 4

    # Since the ssse3 method must examine at least 512 bits while the
    # Gillies method doesn't, this puts the time tradeoff around 210 bits.
    # I decided to save a bit of space and round that up to 224 bits.
    # (Experience will tell us if 256 is a better boundary.)
    if num_bits <= 224:
        return 8

    # If you have POPCNT (and you're using it) then there's no reason
    # to use a larger alignment
    if _has_popcnt:
        if num_bits >= 768:
            if bitops.get_alignment_method("align8-large") == "POPCNT":
                return 8
        else:
            if bitops.get_alignment_method("align8-small") == "POPCNT":
                return 8

    # If you don't have SSSE3 or you aren't using it, then use 8
    if not _has_ssse3 or bitops.get_alignment_method("align-ssse3") != "ssse3":
        return 8

    # In my timing tests:
    #    Lauradoux takes 12.6s
    #    ssse3 takes in 9.0s
    #    Gillies takes 22s


    # Otherwise, go ahead and pad up to 64 bytes
    # (Even at 768 bits/96 bytes, the SSSE3 method is faster.)
    return 64


def fps_to_arena(fps_reader, metadata=None, reorder=True, alignment=None):
    if metadata is None:
        metadata = fps_reader.metadata
    num_bits = metadata.num_bits
    if not num_bits:
        num_bits = metadata.num_bytes * 8
    #assert num_bits

    if alignment is None:
        alignment = get_optimal_alignment(num_bits)

    storage_size = metadata.num_bytes
    if storage_size % alignment != 0:
        n = alignment - storage_size % alignment
        end_padding = "\0" * n
        storage_size += n
    else:
        end_padding = None

    ids = []
    unsorted_fps = StringIO()
    for (id, fp) in fps_reader:
        unsorted_fps.write(fp)
        if end_padding:
            unsorted_fps.write(end_padding)
        ids.append(id)

    unsorted_arena = unsorted_fps.getvalue()
    unsorted_fps.close()
    unsorted_fps = None


    if not reorder or not metadata.num_bits:
        start_padding, end_padding, unsorted_arena = _chemfp.make_unsorted_aligned_arena(
            unsorted_arena, alignment)
        return FingerprintArena(metadata, start_padding, end_padding, storage_size,
                                unsorted_arena, "", ids)

    # Reorder
        
    ordering = (ChemFPOrderedPopcount*len(ids))()
    popcounts = array.array("i", (0,)*(metadata.num_bits+2))

    start_padding, end_padding, unsorted_arena = _chemfp.make_sorted_aligned_arena(
        num_bits, storage_size, unsorted_arena, len(ids),
        ordering, popcounts, alignment)

    new_ids = [ids[item.index] for item in ordering]
    return FingerprintArena(metadata,
                            start_padding, end_padding, storage_size,
                            unsorted_arena, popcounts.tostring(), new_ids)
