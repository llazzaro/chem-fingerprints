import ctypes
from cStringIO import StringIO

import _chemfp

def check_fp_compatibility(query_fp, targets):
    if len(query_fp) != targets.header.num_bytes_per_fp:
        raise TypeError("Query fingerprint size does not match target library")

def check_compatibility(queries, targets):
    if queries.header.num_bits != targets.header.num_bits:
        raise TypeError("Incompatible fingerprint sizes: queries=%d and targets=%d" %
                        (queries.header.num_bits, targets.header.num_bits))



def tanimoto_count_fp(query_fp, targets, threshold):
    check_fp_compatibility(query_fp, targets)

    counts = array.array("i", (0 for i in xrange(len(query_fp))))
    _chemfp.count_tanimoto_arena(threshold, targets.num_bits,
                                 len(query_fp), query_fp, 0, -1,
                                 targets.storage_size, targets.arena, 0, -1,
                                 targets.popcount_indicies,
                                 counts)
    return counts[0]


def tanimoto_count(queries, targets, threshold):
    check_fp_compatibility(query_fp, targets)

    counts = (ctypes.c_int*len(queries))()
    result = ctypes.count_tanimoto_arena(threshold, targets.num_bits,
                                         queries.storage_size, queries.arena, 0, -1,
                                         target.storage_size, target.arena, 0, -1,
                                         target.popcount_indicies,
                                         counts)
    assert result == len(queries), result
    return counts[0]


# Search results stored in a compressed sparse row form

class SearchResultRows(object):
    def __init__(self, offsets, indicies, scores, query_ids, target_ids):
        assert len(offsets) > 0
        if offsets:
            assert indicies[offsets[-2]]
            assert scores[offsets[-2]]
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
        for target_id, end in zip(self.query_ids, self.offsets[1:]):
            yield target_id, zip((target_ids[index] for index in indicies[start:end]),
                                 scores[start:end])
            start = end

    def iterhits():
        indicies = self.indicies
        scores = self.scores
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield target_id, zip(indicies[start:end], scores[start:end])
            start = end

def threshold_tanimoto_search_fp_indicies(query_fp, targets, threshold):
    check_fp_compatibility(query_fp, targets)
    num_bits = targets.num_bits

    offsets = (ctypes.c_int * 2)()
    offsets[0] = 0
    
    num_cells = len(targets)
    indicies = (ctypes.c_int * num_cells)()
    offsets = (ctypes.c_double * num_cells)()

    num_added = _chemfp.threshold_tanimoto_arena(
        threshold, num_bits,
        len(query_fp), query_fp, 0, -1,
        targets.storage_size, targets.arena, 0, -1,
        targets.popcount_indicies,
        offsets, 0,
        indicies, scores)

    assert num_added == 1

    end = offsets[1]
    return [(indicies[i], scores[i]) for i in xrange(end)]

def threshold_tanimoto_search_fp(query_fp, targets, threshold):
    result = threshold_tanimoto_search_fp_indicies(query_fp, targets, threshold)
    return [(targets.ids[i], score) for (index, score) in result]


def threshold_tanimoto_search(queries, targets, threshold):
    check_compatibility(queries, targets)
    num_bits = targets.num_bits

    num_queries = len(queries)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(queries)) * len(targets)
    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()
    
    query_start = 0
    query_end = len(queries)


    def add_rows(query_start):
        return _chemfp.threshold_tanimoto_arena(
            threshold, num_bits,
            queries.storage_size, queries.arena, query_start, query_end,
            targets.storage_size, targets.arena, 0, -1,
            targets.popcount_indicies,
            offsets, query_start,
            indicies, scores)

    return _search(query_end, offsets, indicies, scores, add_rows,
                   queries.ids, targets.ids)

def knearest_tanimoto_search(queries, targets, k, threshold):
    check_compatibility(queries, targets)
    num_bits = queries.header.num_bits

    num_queries = len(queries)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(queries))*k

    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    query_start = 0
    query_end = len(queries)

    def add_rows(query_start):
        return _chemfp.klargest_tanimoto_arena(
            k, threshold, num_bits,
            queries.storage_size, queries.arena, query_start, query_end,
            targets.storage_size, targets.arena, 0, -1,
            targets.popcount_indicies,
            offsets, query_start,
            indicies, scores)

    return _search(query_end, offsets, indicies, scores, add_rows,
                   queries.ids, targets.ids)


# Core of the Tanimoto search routine

def _search(query_end, offsets, indicies, scores, add_rows, query_ids, target_ids):
    num_added = add_rows(0)
    if num_added == query_end:
        return SearchResultRows(offsets, indicies, scores, query_ids, target_ids)

    query_start = num_added

    last = offsets[num_added]
    all_indicies = indicies[:last]
    all_scores = scores[:last]

    while query_start < query_end:
        num_added = add_rows(query_start)
        assert num_added > 0

        prev_last = offsets[query_start]
        all_indicies[prev_last:] = indicies
        all_scores[prev_last:] = scores

        query_start += num_added

    return SearchResultRows(offsets, all_indicies, all_scores, query_ids, target_ids)




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

class Library(object):
    def __init__(self, header, storage_size, arena, popcount_indicies, ids):
        self.header = header
        self.num_bits = header.num_bits
        self.storage_size = storage_size
        self.arena = arena
        self.popcount_indicies = popcount_indicies
        self.ids = ids
        self.fingerprints = FingerprintLookup(header.num_bytes_per_fp, storage_size, arena)
        self._range_check = xrange(len(self))

    def __len__(self):
        return len(self.arena) / self.header.num_bytes_per_fp

    def __getitem__(self):
        start_offset = self._range_check[i] * self._storage_size
        self.arena[start_offset:start_offset+self._fp_size]
        return self.ids[i], arena[start_offset:start_offset+self.header.num_bytes_per_fp]
        
    def reset(self):
        pass

    def __iter__(self):
        storage_size = self.storage_size
        target_fp_size = self.header.num_bytes_per_fp
        arena = self.arena
        for id, start_offset in zip(self.ids, xrange(0, len(self.arena), storage_size)):
            yield id, arena[start_offset:start_offset+target_fp_size]

    _tanimoto_count_fp_ = staticmethod(tanimoto_count_fp)
    _tanimoto_count_ = staticmethod(tanimoto_count)
    _tanimoto_count_once_ = staticmethod(tanimoto_count)

    _threshold_tanimoto_search_fp_ = staticmethod(threshold_tanimoto_search_fp)
    _threshold_tanimoto_search_ = staticmethod(threshold_tanimoto_search)
    _threshold_tanimoto_search_once_ = staticmethod(threshold_tanimoto_search)

#    _knearest_tanimoto_search_fp = staticmethod(knearest_tanimoto_search_fp)
    _knearest_tanimoto_search_ = staticmethod(knearest_tanimoto_search)
    _knearest_tanimoto_search_once_ = staticmethod(knearest_tanimoto_search)


class ChemFPOrderedPopcount(ctypes.Structure):
    _fields_ = [("popcount", ctypes.c_int),
                ("index", ctypes.c_int)]

import array
def reorder_fingerprints(fingerprints):
    ordering = (ChemFPOrderedPopcount*len(fingerprints))()
    popcounts = array.array("i", (0,)*(fingerprints.header.num_bits+1))
    #popcounts = (ctypes.c_int*(fingerprints.header.num_bits+1))()

    new_arena = _chemfp.reorder_by_popcount(
        fingerprints.header.num_bits, fingerprints.storage_size,
        fingerprints.arena, 0, -1, ordering, popcounts)

    new_ids = [fingerprints.ids[item.index] for item in ordering]
    return Library(fingerprints.header, fingerprints.storage_size,
                   new_arena, popcounts.tostring(), new_ids)
                                

def knearest_tanimoto_search_fp(query_fp, targets, k, threshold):
    if not isinstance(query_fp, Library):
        raise Spam
    



def fps_to_library(fps_reader, header=None, sort=True):
    if header is None:
        header = fps_reader.header
    num_bits = header.num_bits
    assert num_bits

    ids = []
    unsorted_fps = StringIO()
    for (id, fp) in fps_reader:
        unsorted_fps.write(fp)
        ids.append(id)

    unsorted_arena = unsorted_fps.getvalue()
    unsorted_fps.close()
    unsorted_fps = None

    fingerprints = Library(header, header.num_bytes_per_fp,
                           unsorted_arena, "", ids)

    if sort:
        return reorder_fingerprints(fingerprints)
    else:
        return fingerprints
