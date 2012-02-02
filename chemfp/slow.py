import heapq
import operator

from . import bitops
from . import Fingerprints
from . import fps_search

def count_tanimoto_hits_fp(query_fp, targets, threshold):
    return sum(1 for target in targets
                   if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

## def iter_count_tanimoto_hits(queries, targets, threshold):
##     for query_id, query_fp in queries:
##         yield query_id, sum(1 for target in targets
##                                if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

def count_tanimoto_hits_arena(queries, targets, threshold):
    counts = []
    for query_id, query_fp in queries:
        counts.append(sum(1 for target in targets
                          if bitops.byte_tanimoto(query_fp, target[1]) >= threshold))
    return counts

def tanimoto_count_once(queries, targets, threshold):
    # Only go through the list of queries once
    results = [0] * len(queries)
    query_ids = []
    query_fps = []
    for (id, fp) in queries:
        query_ids.append(id)
        query_fps.append(fp)
    for target_id, target_fp in targets:
        for i, query_fp in enumerate(query_fps):
            if bitops.byte_tanimoto(query_fp, target_fp) >= threshold:
                results[i] += 1
    return zip(query_ids, results)

##########

def threshold_tanimoto_search_fp(query_fp, targets, threshold):
    ids = []
    scores = []
    for target_id, target_fp in targets:
        score = bitops.byte_tanimoto(query_fp, target_fp)
        if score >= threshold:
            ids.append(target_id)
            scores.append(score)
    return fps_search.FPSSearchResult(ids, scores)

def _iter_threshold_tanimoto_search(queries, targets, threshold):
    for query_id, query_fp in queries:
        yield threshold_tanimoto_search_fp(query_fp, targets, threshold)

def threshold_tanimoto_search_arena(queries, targets, threshold):
    results = []
    for query_id, query_fp in queries:
        results.append(threshold_tanimoto_search_fp(query_fp, targets, threshold))
    return fps_search.FPSSearchResults(results)

#def tanimoto_search_all(queries, targets, threshold):
#    results = [[] for i in xrange(len(queries))]
#    query_fps = [query[1] for query in queries]
#    for target_id, target_fp in targets:
#        for i, query_fp in enumerate(query_fps):
#            score = bitops.byte_tanimoto(query_fp, target_fp)
#            if score >= threshold:
#                results[i].append( (target_id, score) )
#    return results

##########


def knearest_tanimoto_search_fp(query_fp, targets, k, threshold):
    hits = heapq.nlargest(k, threshold_tanimoto_search_fp(query_fp, targets, threshold),
                          key = operator.itemgetter(1))
    if hits:
        ids, scores = zip(*hits)
    else:
        ids, scores = [], []
    return fps_search.FPSSearchResult(ids, scores)

def _iter_knearest_tanimoto_search(queries, targets, k, threshold):
    for hits in _iter_threshold_tanimoto_search(queries, targets, threshold):
        hits = heapq.nlargest(k, hits, key = operator.itemgetter(1))
        ids, scores = zip(*hits)
        yield fps_search.FPSSearchResult(ids, scores)

def knearest_tanimoto_search_arena(queries, targets, k, threshold):
    results = []
    for row in _iter_knearest_tanimoto_search(queries, targets, k, threshold):
        results.append(row)
    return fps_search.FPSSearchResults(results)

# I am not going to optimize this.
#def knearest_tanimoto_search_all(queries, targets, k, threshold):
#    results = tanimoto_search(queries, targets, threshold)
#    for i, (id, hits) in enumerate(results):
#        results[i] = (id, heapq.nlargest(k, hits, key = operator.itemgetter(1)))
#    return results

def _check_threshold(threshold):
    if not (0.0 <= threshold <= 1.0):
        raise ValueError("threshold must between 0.0 and 1.0, inclusive")

class SlowFingerprints(Fingerprints):
    def count_tanimoto_hits_fp(self, query_fp, threshold=0.7):
        if not (0.0 <= threshold <= 1.0):
            raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return count_tanimoto_hits_fp(query_fp, self._id_fp_pairs, threshold)

    def count_tanimoto_hits_arena(self, queries, threshold=0.7):
        _check_threshold(threshold)
        return count_tanimoto_hits_arena(queries, self._id_fp_pairs, threshold)

    def threshold_tanimoto_search_fp(self, fp, threshold=0.7):
        _check_threshold(threshold)
        return threshold_tanimoto_search_fp(fp, self._id_fp_pairs, threshold)

    def threshold_tanimoto_search_arena(self, queries, threshold=0.7):
        _check_threshold(threshold)
        return threshold_tanimoto_search_arena(queries, self._id_fp_pairs, threshold)

    def knearest_tanimoto_search_fp(self, fp, k=3, threshold=0.7):
        if k < 0:
            raise ValueError("k must be non-negative")
        _check_threshold(threshold)
        return knearest_tanimoto_search_fp(fp, self._id_fp_pairs, k, threshold)

    def knearest_tanimoto_search_arena(self, queries, k=3, threshold=0.7, batch_size=100):
        if k < 0:
            raise ValueError("k must be non-negative")
        _check_threshold(threshold)
        return knearest_tanimoto_search_arena(queries, self._id_fp_pairs, k, threshold)
