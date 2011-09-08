# Internal module to help with FPS-based searches
from __future__ import absolute_import

import ctypes
import itertools
import array

import _chemfp
from . import check_fp_problems, check_metadata_problems

def report_errors(problem_report):
    for (severity, error, msg_template) in problem_report:
        if severity == "error":
            raise TypeError(msg_template.format(metadata1 = "query",
                                                metadata2 = "target"))

######## count Tanimoto search #########

def _fp_to_arena(query_fp, metadata):
    assert len(query_fp) == metadata.num_bytes
    from . import arena
    return arena.Library(metadata, len(query_fp), query_fp, "", [None])

def count_tanimoto_hits_fp(query_fp, target_reader, threshold):
    return count_tanimoto_hits_arena(_fp_to_arena(query_fp, target_reader.metadata), threshold)[0]

def count_tanimoto_hits_arena(query_arena, target_reader, threshold):
    counts = array.array("i", (0 for i in xrange(len(query_arena))))

    lineno = target_reader._first_fp_lineno
    num_bits = target_reader.metadata.num_bits

    for block in target_reader.iter_blocks():
        err, num_lines = _chemfp.fps_count_tanimoto_hits(
            num_bits, query_arena.storage_size, query_arena.arena, 0, -1,
            block, 0, -1,
            threshold, counts)
        lineno += num_lines
        if err:
            ERROR

    return list(counts)
    



######## threshold Tanimoto search #########

class TanimotoCell(ctypes.Structure):
    _fields_ = [("score", ctypes.c_double),
                ("query_index", ctypes.c_int),
                ("id_start", ctypes.c_int),
                ("id_end", ctypes.c_int)]


def threshold_tanimoto_search_fp(query_fp, target_reader, threshold):
    hits = []

    fp_size = len(query_fp)
    num_bits = fp_size * 8
        
    NUM_CELLS = 1000
    cells = (TanimotoCell*NUM_CELLS)()

    lineno = target_reader._first_fp_lineno
    
    for block in target_reader.iter_blocks():
        start = 0
        end = len(block)
        while 1:
            err, start, num_lines, num_cells = _chemfp.fps_threshold_tanimoto_search(
                num_bits, fp_size, query_fp, 0, -1,
                block, start, end,
                threshold, cells)
            lineno += num_lines
            if err:
                raise TypeError("Error at or after line %r of %r: %s" %
                                (lineno, target_reader, _chemfp.strerror(err)))
                
            for cell in itertools.islice(cells, 0, num_cells):
                    id = block[cell.id_start:cell.id_end]
                    hits.append( (id, cell.score) )
            if start == end:
                break
    return hits

def threshold_tanimoto_search_all(query_arena, target_reader, threshold):
    all_hits = [[] for i in xrange(len(query_arena))]
    if not all_hits:
        return []

    
    # Compute at least 100 tanimotos per query, but at most 10,000 at a time
    # (That's about 200K of memory)
    NUM_CELLS = max(10000, len(query_arena) * 100)
    cells = (TanimotoCell*NUM_CELLS)()

    lineno = target_reader._first_fp_lineno

    for block in target_reader.iter_blocks():
        start = 0
        end = len(block)
        while 1:
            err, start, num_lines, num_cells = _chemfp.fps_threshold_tanimoto_search(
                query_arena.metadata.num_bits, query_arena.storage_size,
                query_arena.arena, 0, -1,
                block, start, end,
                threshold, cells)
            lineno += num_lines
            if err:
                ERROR
                
            for cell in itertools.islice(cells, 0, num_cells):
                id = block[cell.id_start:cell.id_end]
                all_hits[cell.query_index].append( (id, cell.score) )

            if start == end:
                break

    return all_hits
            
######### k-nearest Tanimoto search, with threshold


# Support for peering into the chemfp_fps_heap data structure

def _make_knearest_search(num_queries, k):
    class TanimotoHeap(ctypes.Structure):
        _fields_ = [("size", ctypes.c_int),
                    ("heap_state", ctypes.c_int),
                    ("indicies", ctypes.POINTER(ctypes.c_int*k)),
                    ("ids", ctypes.POINTER(ctypes.c_char_p*k)),
                    ("scores", ctypes.POINTER(ctypes.c_double*k))]

    class KNearestSearch(ctypes.Structure):
        _fields_ = [("queries_start", ctypes.c_char_p),
                    ("num_queries", ctypes.c_int),
                    ("query_fp_size", ctypes.c_int),
                    ("query_storage_size", ctypes.c_int),
                    ("k", ctypes.c_int),
                    ("search_state", ctypes.c_int),
                    ("threshold", ctypes.c_double),
                    ("heaps", ctypes.POINTER(TanimotoHeap*num_queries)),
                    ("num_targets_processed", ctypes.c_int),
                    ("_all_ids", ctypes.c_void_p),
                    ("_all_scores", ctypes.c_void_p)]

    return KNearestSearch()


def check_num_bits_compatibility(query_arena, target_reader):
    num_bits = query_arena.metadata.num_bits
    if num_bits != target_reader.metadata.num_bits:
        raise TypeError("The query fingerprints have %d bits while the targets have %d" %
                        (num_bits, target_reader.metadata.num_bits))
    return num_bits

def knearest_tanimoto_search_fp(query_fp, target_reader, k, threshold):
    query_arena = _fp_to_arena(query_fp, target_reader.metadata)
    return knearest_tanimoto_search_all(query_arena, target_reader, k, threshold)[0]

def knearest_tanimoto_search_all(query_arena, target_reader, k, threshold):
    report_errors(query_arena.metadata, target_reader.metadata)
    num_bits = query_arena.metadata.num_bits

    num_queries = len(query_arena)
    search = _make_knearest_search(num_queries, k)

    _chemfp.fps_knearest_search_init(
        search,
        num_bits, query_arena.storage_size, query_arena.arena, 0, -1,
        k, threshold)

    try:
        for block in target_reader.iter_blocks():
            err = _chemfp.fps_knearest_tanimoto_search_feed(search, block)
            if err:
                ERR

        _chemfp.fps_knearest_search_finish(search)

        all_results = []
        for query_index in xrange(num_queries):
            heap = search.heaps[0][query_index]
            results = []
            for i in xrange(heap.size):
                id = ctypes.string_at(heap.ids[0][i])
                score = heap.scores[0][i]
                results.append( (id, score) )
            all_results.append(results)
        return all_results

    finally:
        _chemfp.fps_knearest_search_free(search)
