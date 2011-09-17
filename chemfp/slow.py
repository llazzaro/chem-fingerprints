import heapq
import operator

from . import bitops
from . import Fingerprints

def count_tanimoto_hits_fp(query_fp, targets, threshold):
    return sum(1 for target in targets
                   if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

## def iter_count_tanimoto_hits(queries, targets, threshold):
##     for query_id, query_fp in queries:
##         yield query_id, sum(1 for target in targets
##                                if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

def count_tanimoto_hits(queries, targets, threshold):
    results = []
    for query_id, query_fp in queries:
        results.append(sum(1 for target in targets
                           if bitops.byte_tanimoto(query_fp, target[1]) >= threshold))
    return results

def tanimoto_count_once(queries, targets, threshold):
    results = [0] * len(queries)
    query_fps = [query[1] for query in queries]
    for target_id, target_fp in targets:
        for i, query_fp in enumerate(query_fps):
            if bitops.byte_tanimoto(query_fp, target_fp) >= threshold:
                results[i] += 1
    return results

##########

def threshold_tanimoto_search_fp(query_fp, targets, threshold):
    results = []
    for target_id, target_fp in targets:
        score = bitops.byte_tanimoto(query_fp, target_fp)
        if score >= threshold:
            results.append( (target_id, score) )
    return results

def iter_threshold_tanimoto_search(queries, targets, threshold):
    for query_id, query_fp in queries:
        yield threshold_tanimoto_search_fp(query_fp, targets, threshold)

def threshold_tanimoto_search(queries, targets, threshold):
    return list(iter_threshold_tanimoto_search(queries, targets, threshold))

def tanimoto_search_all(queries, targets, threshold):
    results = [[] for i in xrange(len(queries))]
    query_fps = [query[1] for query in queries]
    for target_id, target_fp in targets:
        for i, query_fp in enumerate(query_fps):
            score = bitops.byte_tanimoto(query_fp, target_fp)
            if score >= threshold:
                results[i].append( (target_id, score) )
    return results

##########


def knearest_tanimoto_search_fp(query_fp, targets, k, threshold):
    return heapq.nlargest(k, threshold_tanimoto_search_fp(query_fp, targets, threshold),
                          key = operator.itemgetter(1))

def iter_knearest_tanimoto_search(queries, targets, k, threshold):
    for hits in iter_threshold_tanimoto_search(queries, targets, threshold):
        yield heapq.nlargest(k, hits, key = operator.itemgetter(1))

def knearest_tanimoto_search(queries, targets, k, threshold):
    return list(iter_knearest_tanimoto_search(queries, targets, k, threshold))

# I am not going to optimize this.
def knearest_tanimoto_search_all(queries, targets, k, threshold):
    results = tanimoto_search(queries, targets, threshold)
    for i, (id, hits) in enumerate(results):
        results[i] = (id, heapq.nlargest(k, hits, key = operator.itemgetter(1)))
    return results

class SlowFingerprints(Fingerprints):
    def count_tanimoto_hits_fp(self, query_fp, threshold=0.7):
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return count_tanimoto_hits_fp(query_fp, self._id_fp_pairs, threshold)

    def count_tanimoto_hits_arena(self, query_arena, threshold=0.7):
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return count_tanimoto_hits(query_arena, self._id_fp_pairs, threshold)

    def threshold_tanimoto_search_fp(self, fp, threshold=0.7):
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return threshold_tanimoto_search_fp(fp, self._id_fp_pairs, threshold)

    def threshold_tanimoto_search_arena(self, query_arena, threshold=0.7):
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return threshold_tanimoto_search(query_arena, self._id_fp_pairs, threshold)

    def knearest_tanimoto_search_fp(self, fp, k=3, threshold=0.7):
        if k < 0: raise ValueError("k must be non-negative")
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return knearest_tanimoto_search_fp(fp, self._id_fp_pairs, k, threshold)

    def knearest_tanimoto_search_arena(self, query_arena, k=3, threshold=0.7):
        if k < 0: raise ValueError("k must be non-negative")
        if not (0.0 <= threshold <= 1.0): raise ValueError("threshold must between 0.0 and 1.0, inclusive")
        return knearest_tanimoto_search(query_arena, self._id_fp_pairs, k, threshold)
    
