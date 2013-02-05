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

from chemfp import FingerprintReader
import _chemfp
from chemfp import bitops, search

__all__ = []
    
class FingerprintArena(FingerprintReader):
    """Stores fingerprints in a contiguous block of memory

    The public attributes are:
       metadata
           `Metadata` about the fingerprints
       ids
           list of identifiers, ordered by position
    """
    def __init__(self, metadata, alignment,
                 start_padding, end_padding, storage_size, arena,
                 popcount_indices, arena_ids, start=0, end=None,
                 id_lookup=None,
                 ):
        if metadata.num_bits is None:
            raise TypeError("Missing metadata num_bits information")
        if metadata.num_bytes is None:
            raise TypeError("Missing metadata num_bytes information")
        self.metadata = metadata
        self.alignment = alignment
        self.num_bits = metadata.num_bits
        self.start_padding = start_padding
        self.end_padding = end_padding
        self.storage_size = storage_size
        self.arena = arena
        self.popcount_indices = popcount_indices
        self.arena_ids = arena_ids
        self.start = start   # the starting index in the arena (not byte position!)
        if end is None:      # the ending index in the arena (not byte position!)
            if self.metadata.num_bytes:
                end = (len(arena) - start_padding - end_padding) // self.storage_size
            else:
                end = 0
        self.end = end
        if self.start == 0 and self.end == len(arena_ids):
            self._ids = arena_ids
        else:
            self._ids = None
        self._id_lookup = id_lookup
        assert end >= start
        self._range_check = xrange(end-start)

    def __len__(self):
        """Number of fingerprint records in the FingerprintArena"""
        return self.end - self.start

    @property
    def ids(self):
        ids = self._ids
        if ids is None:
            ids = self.arena_ids[self.start:self.end]
            self._ids = ids
        return ids

    def __getitem__(self, i):
        """Return the (id, fingerprint) at position i"""
        if isinstance(i, slice):
            start, end, step = i.indices(self.end - self.start)
            if step != 1:
                raise IndexError("arena slice step size must be 1")
            if start >= end:
                return FingerprintArena(self.metadata, self.alignment,
                                        0, 0, self.storage_size, "",
                                        "", [], 0, 0)
            return FingerprintArena(self.metadata, self.alignment,
                                    self.start_padding, self.end_padding,
                                    self.storage_size, self.arena,
                                    self.popcount_indices, self.arena_ids,
                                    self.start+start, self.start+end)
        try:
            i = self._range_check[i]
        except IndexError:
            raise IndexError("arena fingerprint index out of range")
        arena_i = i + self.start
        start_offset = arena_i * self.storage_size + self.start_padding
        end_offset = start_offset + self.metadata.num_bytes
        return self.arena_ids[arena_i], self.arena[start_offset:end_offset]

    def _make_id_lookup(self):
        d = dict((id, i) for (i, id) in enumerate(self.ids))
        self._id_lookup = d.get
        return self._id_lookup
        
    def get_by_id(self, id):
        """Given the record identifier, return the (id, fingerprint) tuple or None if not present"""
        id_lookup = self._id_lookup
        if id_lookup is None:
            id_lookup = self._make_id_lookup()
        i = id_lookup(id)
        if i is None:
            return None
        arena_i = i + self.start
        start_offset = arena_i * self.storage_size + self.start_padding
        end_offset = start_offset + self.metadata.num_bytes
        return self.arena_ids[arena_i], self.arena[start_offset:end_offset]

    def get_index_by_id(self, id):
        """Given the record identifier, return the record index or None if not present"""
        id_lookup = self._id_lookup
        if id_lookup is None:
            id_lookup = self._make_id_lookup()
        return id_lookup(id)

    def get_fingerprint_by_id(self, id):
        """Given the record identifier, return its fingerprint or None if not present"""
        id_lookup = self._id_lookup
        if id_lookup is None:
            id_lookup = self._make_id_lookup()
        i = id_lookup(id)
        if i is None:
            return None
        arena_i = i + self.start
        start_offset = arena_i * self.storage_size + self.start_padding
        end_offset = start_offset + self.metadata.num_bytes
        return self.arena[start_offset:end_offset]

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
            try:
                for i, (id, fp) in enumerate(self):
                    io.write_fps1_fingerprint(output, fp, id)
            except ValueError, err:
                raise ValueError("%s in record %i" % (err, i+1))
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
        for id, start_offset in zip(self.arena_ids[self.start:self.end],
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
            end = start+arena_size
            if end > self.end:
                end = self.end
            yield FingerprintArena(self.metadata, self.alignment,
                                   self.start_padding, self.end_padding,
                                   storage_size, self.arena,
                                   self.popcount_indices, self.arena_ids, start, end)
            start = end

    def count_tanimoto_hits_fp(self, query_fp, threshold=0.7):
        """Count the fingerprints which are similar enough to the query fingerprint

        DEPRECATED: Use chemfp.search.count_tanimoto_hits_fp
        
        Return the number of fingerprints in this arena which are
        at least `threshold` similar to the query fingerprint `query_fp`.

        :param query_fp: query fingerprint
        :type query_fp: byte string
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: integer count
        """
        return search.count_tanimoto_hits_fp(query_fp, self, threshold)

    def count_tanimoto_hits_arena(self, queries, threshold=0.7):
        """Count the fingerprints which are similar enough to each query fingerprint

        DEPRECATED: Use chemfp.search.count_tanimoto_hits_arena
        
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
        """
        return search.count_tanimoto_hits_arena(queries, self, threshold)

    def threshold_tanimoto_search_fp(self, query_fp, threshold=0.7):
        """Find the fingerprints which are similar enough to the query fingerprint

        DEPRECATED: Use chemfp.search.threshold_tanimoto_search_fp

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
        return search.threshold_tanimoto_search_fp(query_fp, self, threshold)

    def threshold_tanimoto_search_arena(self, queries, threshold=0.7, arena_size=100):
        """Find the fingerprints which are similar to each of the query fingerprints

        DEPRECATED: Use chemfp.search.threshold_tanimoto_search_arena

        For each fingerprint in the `query_arena`, find all of the
        fingerprints in this arena which are at least `threshold`
        similar. The hits are returned as a `SearchResults` instance.
        
        :param query_arena: query arena
        :type query_arena: FingerprintArena
        :param threshold: minimum similarity threshold (default: 0.7)
        :type threshold: float between 0.0 and 1.0, inclusive
        :returns: SearchResults
        """
        return search.threshold_tanimoto_search_arena(queries, self, threshold)

    def knearest_tanimoto_search_fp(self, query_fp, k=3, threshold=0.7):
        """Find the k-nearest fingerprints which are similar to the query fingerprint

        DEPRECATED: Use chemfp.search.knearest_tanimoto_search_fp
        
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
        return search.knearest_tanimoto_search_fp(query_fp, self, k, threshold)

    def knearest_tanimoto_search_arena(self, queries, k=3, threshold=0.7):
        """Find the k-nearest fingerprint which are similar to each of the query fingerprints

        DEPRECATED: Use chemfp.search.knearest_tanimoto_search_arena or chemfp.search.knearest_tanimoto_search_symmetric

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
        return search.knearest_tanimoto_search_arena(queries, self, k, threshold)

    def copy(self, indices=None, reorder=None):
        if reorder is None:
            if indices is None:
                # This is a pure copy. Reorder only if there are popcount indices.
                reorder = (self.popcount_indices != "")
            else:
                # The default is to go fast. If you want to preserve index order
                # then you'll need to set reorder=False
                reorder = True
            
        if indices is None:
            # Make a completely new arena
            # Handle the trivial case where I don't need to do anything.
            if (self.start == 0 and
                (self.end*self.storage_size + self.start_padding + self.end_padding == len(self.arena)) and
                (not reorder or self.popcount_indices)):
                return FingerprintArena(self.metadata, self.alignment,
                                        self.start_padding, self.end_padding, self.storage_size, self.arena,
                                        self.popcount_indices, self.arena_ids,
                                        start = 0, end = self.end,
                                        id_lookup = self._id_lookup)
            
            # Otherwise I need to do some work
            # Make a copy of the actual fingerprints. (Which could be a subarena.)
            start = self.start_padding + self.start*self.storage_size
            end = self.start_padding + self.end*self.storage_size
            arena = self.arena[start:end]

            # If we don't have popcount_indices and don't want them ordered
            # then just do the alignment and we're done.
            if not reorder and not self.popcount_indices:
                # Don't reorder the unordered fingerprints
                start_padding, end_padding, unsorted_arena = (
                    _chemfp.make_unsorted_aligned_arena(arena, self.alignment))
                return FingerprintArena(self.metadata, self.alignment, start_padding, end_padding,
                                        self.storage_size, unsorted_arena, "", self.ids,
                                        id_lookup = self._id_lookup)

            # Either we're already sorted or we should become sorted.
            # If we're sorted then make_sorted_aligned_arena will detect
            # that and keep the old arena. Otherwise it sorts first and
            # makes a new arena block.
            current_ids = self.ids
            ordering = (ChemFPOrderedPopcount*len(current_ids))()
            popcounts = array.array("i", (0,)*(self.metadata.num_bits+2))
            start_padding, end_padding, arena = _chemfp.make_sorted_aligned_arena(
                self.metadata.num_bits, self.storage_size, arena, len(current_ids),
                ordering, popcounts, self.alignment)

            reordered_ids = [current_ids[item.index] for item in ordering]
            return FingerprintArena(self.metadata, self.alignment,
                                    start_padding, end_padding, self.storage_size,
                                    arena, popcounts.tostring(), reordered_ids)

        # On this pathway, we want to make a new arena which contains
        # selected fingerprints given indices into the old arena.
        
        arena = self.arena
        storage_size = self.storage_size
        start = self.start
        start_padding = self.start_padding
        arena_ids = self.arena_ids
        
        # First make sure that all of the indices are in range.
        # This will also convert negative indices into positive ones.
        new_indices = []
        range_check = self._range_check
        try:
            for i in indices:
                new_indices.append(range_check[i])
        except IndexError:
            raise IndexError("arena fingerprint index %d is out of range" % (i,))

        if reorder and self.popcount_indices:
            # There's a slight performance benefit because
            # make_sorted_aligned_arena will see that the fingerprints
            # are already in sorted order and not resort.
            # XXX Is that true? Why do a Python sort instead of a C sort?
            # Perhaps because then I don't need to copy fingerprints?
            new_indices.sort()

        # Copy the fingerprints over to a new arena block
        unsorted_fps = []
        new_ids = []
        for new_i in new_indices:
            start_offset = start_padding + new_i*storage_size
            end_offset = start_offset + storage_size
            unsorted_fps.append(arena[start_offset:end_offset])
            new_ids.append(arena_ids[new_i])
                
        unsorted_arena = "".join(unsorted_fps)
        unsorted_fps = None   # regain some memory

        # If the caller doesn't want ordered data, then leave it unsorted
        if not reorder:
            start_padding, end_padding, unsorted_arena = _chemfp.make_unsorted_aligned_arena(
                unsorted_arena, self.alignment)
            return FingerprintArena(self.metadata, self.alignment, start_padding, end_padding, storage_size,
                                    unsorted_arena, "", new_ids)

        # Otherwise, reorder and align the area, along with popcount information
        ordering = (ChemFPOrderedPopcount*len(new_ids))()
        popcounts = array.array("i", (0,)*(self.metadata.num_bits+2))

        start_padding, end_padding, sorted_arena = _chemfp.make_sorted_aligned_arena(
            self.metadata.num_bits, storage_size, unsorted_arena, len(new_ids),
            ordering, popcounts, self.alignment)

        reordered_ids = [new_ids[item.index] for item in ordering]
        return FingerprintArena(self.metadata, self.alignment,
                                start_padding, end_padding, storage_size,
                                sorted_arena, popcounts.tostring(), reordered_ids)
        

# TODO: push more of this malloc-management down into C
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
        if metadata.num_bytes is None:
            raise ValueError("metadata must contain at least one of num_bits or num_bytes")
        num_bits = metadata.num_bytes * 8
    #assert num_bits

    if alignment is None:
        alignment = get_optimal_alignment(num_bits)

    num_bytes = metadata.num_bytes

    storage_size = num_bytes
    if storage_size % alignment != 0:
        n = alignment - storage_size % alignment
        end_padding = "\0" * n
        storage_size += n
    else:
        end_padding = None

    ids = []
    unsorted_fps = StringIO()
    for (id, fp) in fps_reader:
        if len(fp) != num_bytes:
            raise ValueError("Fingerprint for id %r has %d bytes while the metadata says it should have %d"
                             % (id, len(fp), num_bytes))
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
        return FingerprintArena(metadata, alignment, start_padding, end_padding, storage_size,
                                unsorted_arena, "", ids)

    # Reorder
        
    ordering = (ChemFPOrderedPopcount*len(ids))()
    popcounts = array.array("i", (0,)*(metadata.num_bits+2))

    start_padding, end_padding, unsorted_arena = _chemfp.make_sorted_aligned_arena(
        num_bits, storage_size, unsorted_arena, len(ids),
        ordering, popcounts, alignment)

    new_ids = [ids[item.index] for item in ordering]
    return FingerprintArena(metadata, alignment,
                            start_padding, end_padding, storage_size,
                            unsorted_arena, popcounts.tostring(), new_ids)
