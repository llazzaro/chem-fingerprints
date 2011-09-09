from __future__ import absolute_import

import ctypes
from cStringIO import StringIO

from chemfp import _THRESHOLD, _K, FingerprintReader, check_fp_problems, check_metadata_problems
import _chemfp

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
                                 len(query_fp), query_fp, 0, -1,
                                 target_arena.storage_size, target_arena.arena,
                                 target_arena.start, target_arena.end,
                                 target_arena.popcount_indicies,
                                 counts)
    return counts[0]


def count_tanimoto_hits_arena(query_arena, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)

    counts = (ctypes.c_int*len(query_arena))()
    _chemfp.count_tanimoto_arena(threshold, target_arena.num_bits,
                                 query_arena.storage_size,
                                 query_arena.arena, query_arena.start, query_arena.end,
                                 target_arena.storage_size,
                                 target_arena.arena, target_arena.start, target_arena.end,
                                 target_arena.popcount_indicies,
                                 counts)
    return counts


# Search results stored in a compressed sparse row form

class SearchHits(object):
    def __init__(self, offsets, indicies, scores, query_ids, target_ids):
        assert len(offsets) > 0
        self.offsets = offsets
        self.indicies = indicies
        self.scores = scores
        self.query_ids = query_ids
        self.target_ids = target_ids
        
    def __len__(self):
        return len(self.offsets)-1

    def size(self, i):
        return self.offsets[i+1]-self.offsets[i]
    
    def __getitem__(self, i):
        start, end = self.offsets[i:i+2]
        return zip(self.indicies[start:end], self.scores[start:end])

    def __iter__(self):
        target_ids = self.target_ids
        indicies = self.indicies
        scores = self.scores
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield zip((target_ids[index] for index in indicies[start:end]),
                      scores[start:end])
            start = end

    def iterhits():
        indicies = self.indicies
        scores = self.scores
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield target_id, zip(indicies[start:end], scores[start:end])
            start = end

def threshold_tanimoto_search_fp_indicies(query_fp, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)
    num_bits = target_arena.num_bits

    offsets = (ctypes.c_int * 2)()
    offsets[0] = 0
    
    num_cells = len(target_arena)
    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    num_added = _chemfp.threshold_tanimoto_arena(
        threshold, num_bits,
        len(query_fp), query_fp, 0, -1,
        target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
        target_arena.popcount_indicies,
        offsets, 0,
        indicies, scores)

    assert num_added == 1

    end = offsets[1]
    return [(indicies[i], scores[i]) for i in xrange(end)]

def threshold_tanimoto_search_fp(query_fp, target_arena, threshold):
    require_matching_fp_size(query_fp, target_arena)
    result = threshold_tanimoto_search_fp_indicies(query_fp, target_arena, threshold)
    return [(target_arena.ids[index], score) for (index, score) in result]


def threshold_tanimoto_search_arena(query_arena, target_arena, threshold):
    require_matching_sizes(query_arena, target_arena)
    num_bits = target_arena.num_bits

    num_queries = len(query_arena)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, num_queries) * len(target_arena)
    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()
    
    query_start = query_arena.start
    query_end = query_arena.end


    def add_rows(query_start, offset_start):
        return _chemfp.threshold_tanimoto_arena(
            threshold, num_bits,
            query_arena.storage_size, query_arena.arena, query_start, query_end,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indicies,
            offsets, offset_start, # XXX should query_start=0?
            indicies, scores)

    return _search(query_start, query_end, offsets, indicies, scores, add_rows,
                   query_arena.ids, target_arena.ids)

def knearest_tanimoto_search_arena(query_arena, target_arena, k, threshold):
    require_matching_sizes(query_arena, target_arena)
    num_bits = query_arena.metadata.num_bits

    num_queries = len(query_arena)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(query_arena))*k

    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    query_start = query_arena.start
    query_end = query_arena.end

    def add_rows(query_start, offset_start):
        return _chemfp.knearest_tanimoto_arena(
            k, threshold, num_bits,
            query_arena.storage_size, query_arena.arena, query_start, query_end,
            target_arena.storage_size, target_arena.arena, target_arena.start, target_arena.end,
            target_arena.popcount_indicies,
            offsets, offset_start,
            indicies, scores)

    return _search(query_start, query_end, offsets, indicies, scores, add_rows,
                   query_arena.ids, target_arena.ids)


# Core of the Tanimoto search routine

def _search(query_start, query_end, offsets, indicies, scores,
            add_rows, query_ids, target_ids):
    num_added = add_rows(query_start, 0)
    if num_added == query_end:
        return SearchHits(offsets, indicies, scores, query_ids, target_ids)

    query_start = query_start + num_added
    offset_start = num_added

    last = offsets[num_added]
    all_indicies = indicies[:last]
    all_scores = scores[:last]

    while query_start < query_end:
        num_added = add_rows(query_start, offset_start)
        assert num_added > 0

        prev_last = offsets[query_start]
        all_indicies[prev_last:] = indicies
        all_scores[prev_last:] = scores

        query_start += num_added

    return SearchHits(offsets, all_indicies, all_scores, query_ids, target_ids)




class FingerprintLookup(object):
    def __init__(self, fp_size, storage_size, arena):
        self._fp_size = fp_size
        self._storage_size = storage_size
        self._arena = arena
        self._range_check = xrange(len(self))

    def __len__(self):
        return len(self._arena) / self._storage_size

    def __iter__(self):
        fp_size = self._fp_size
        arena = self._arena
        for id, start_offset in zip(self.ids, xrange(0, len(arena), storage_size)):
            yield id, arena[start_offset:start_offset+target_fp_size]
        
        
    def __getitem__(self, i):
        start_offset = self._range_check[i] * self._storage_size
        return self._arena[start_offset:start_offset+self._fp_size]

class FingerprintArena(FingerprintReader):
    def __init__(self, metadata, storage_size, arena, popcount_indicies, ids,
                 start=0, end=None):
        if metadata.num_bits is None:
            raise TypeError("Missing metadata num_bits information")
        if metadata.num_bytes is None:
            raise TypeError("Missing metadata num_bytes information")
        self.metadata = metadata
        self.num_bits = metadata.num_bits
        self.storage_size = storage_size
        self.arena = arena
        self.popcount_indicies = popcount_indicies
        self.ids = ids
        self.fingerprints = FingerprintLookup(metadata.num_bytes, storage_size, arena)
        self.start = start
        if end is None:
            end = len(arena) // self.metadata.num_bytes
        self.end = end
        assert end >= start
        self._range_check = xrange(end-start)

    def __len__(self):
        return self.end - self.start

    def __getitem__(self, i):
        i = self._range_check[i]
        arena_i = i + self.start
        start_offset = arena_i * self._storage_size
        end_offset = start_offset + self.metadata.num_bytes
        return self.ids[i], self.arena[start_offset:end_offset]


    def save(self, destination):
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
        pass

    def __iter__(self):
        storage_size = self.storage_size
        target_fp_size = self.metadata.num_bytes
        arena = self.arena
        for id, start_offset in zip(self.ids, xrange(self.start*storage_size,
                                                     self.end*storage_size, storage_size)):
            yield id, arena[start_offset:start_offset+target_fp_size]

    def iter_arenas(self, arena_size):
        if arena_size is None:
            yield self
            return
        
        storage_size = self.storage_size
        start = self.start
        for i in xrange(0, len(self), arena_size):
            ids = self.ids[i:i+arena_size]
            end = start+len(ids)
            yield FingerprintArena(self.metadata, self.storage_size, self.arena,
                                   self.popcount_indicies, ids, start, end)
            start = end

    def count_tanimoto_hits_fp(self, query_fp, threshold=_THRESHOLD):
        return count_tanimoto_hits_fp(query_fp, self, threshold)

    def count_tanimoto_hits_arena(self, query_arena, threshold=_THRESHOLD):
        return count_tanimoto_hits_arena(query_arena, self, threshold)

    def threshold_tanimoto_search_fp(self, query_fp, threshold=_THRESHOLD):
        return threshold_tanimoto_search_fp(query_fp, self, threshold)

    def threshold_tanimoto_search_arena(self, query_arena, threshold=_THRESHOLD):
        return threshold_tanimoto_search_arena(query_arena, self, threshold)

    def knearest_tanimoto_search_fp(self, query_fp, k=_K, threshold=_THRESHOLD):
        return knearest_tanimoto_search_fp(query_fp, self, k, threshold)

    def knearest_tanimoto_search_arena(self, query_arena, k=_K, threshold=_THRESHOLD):
        return knearest_tanimoto_search_arena(query_arena, self, k, threshold)


class ChemFPOrderedPopcount(ctypes.Structure):
    _fields_ = [("popcount", ctypes.c_int),
                ("index", ctypes.c_int)]

import array
def reorder_fingerprints(fingerprints):
    ordering = (ChemFPOrderedPopcount*len(fingerprints))()
    popcounts = array.array("i", (0,)*(fingerprints.metadata.num_bits+1))

    new_arena = _chemfp.reorder_by_popcount(
        fingerprints.metadata.num_bits, fingerprints.storage_size,
        fingerprints.arena, fingerprints.start, fingerprints.end, ordering, popcounts)

    new_ids = [fingerprints.ids[item.index] for item in ordering]
    return FingerprintArena(fingerprints.metadata, fingerprints.storage_size,
                            new_arena, popcounts.tostring(), new_ids)
                                


def fps_to_arena(fps_reader, metadata=None, reorder=True):
    if metadata is None:
        metadata = fps_reader.metadata
    num_bits = metadata.num_bits
    if not num_bits:
        num_bits = metadata.num_bytes * 8
    assert num_bits

    ids = []
    unsorted_fps = StringIO()
    for (id, fp) in fps_reader:
        unsorted_fps.write(fp)
        ids.append(id)

    unsorted_arena = unsorted_fps.getvalue()
    unsorted_fps.close()
    unsorted_fps = None

    fingerprints = FingerprintArena(metadata, metadata.num_bytes,
                                    unsorted_arena, "", ids)

    if reorder:
        return reorder_fingerprints(fingerprints)
    else:
        return fingerprints
