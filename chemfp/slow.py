import heapq
from . import bitops

def tanimoto_count_fp(query_fp, targets, threshold):
    return sum(1 for target in targets
                   if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

def iter_tanimoto_count(queries, targets, threshold):
    for query_id, query_fp in queries:
        yield query_id, sum(1 for target in targets
                               if bitops.byte_tanimoto(query_fp, target[1]) >= threshold)

def tanimoto_count(queries, targets, threshold):
    return list(iter_tanimoto_count(queries, targets, threshold))

def tanimoto_count_once(queries, targets, threshold):
    results = [0] * len(queries)
    query_fps = [query[1] for query in queries]
    for target_id, target_fp in targets:
        for i, query_fp in enumerate(query_fps):
            if bitops.byte_tanimoto(query_fp, target_fp) >= threshold:
                results[i] += 1
    return results

##########

def tanimoto_search_fp(query_fp, targets, threshold):
    results = []
    for target_id, target_fp in targets:
        score = bitops.byte_tanimoto(query_fp, target_fp)
        if score >= threshold:
            results.append( (target_id, score) )
    return results

def tanimoto_search_iter(queries, targets, threshold):
    for query_id, query_fp in queries:
        yield (query_id, tanimoto_search_fp(query_fp, targets, threshold))

def tanimoto_search(queries, targets, threshold):
    return list(tanimoto_search_iter(queries, targets, threshold))

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
    return heapq.nlargest(k, tanimoto_search_fp(query_fp, targets, threshold),
                          key = operator.itemgetter(1))

def knearest_tanimoto_search_iter(queries, targets, k, threshold):
    for id, hits in tanimoto_search_iter(queries, targets, threshold):
        print id, hits
        yield id, heapq.nlargest(k, hits, key = operator.itemgetter(1))

def knearest_tanimoto_search(queries, targets, k, threshold):
    return list(knearest_tanimoto_search_iter(queries, targets, k, threshold))

# I am not going to optimize this.
def knearest_tanimoto_search_all(queries, targets, k, threshold):
    results = tanimoto_search(queries, targets, threshold)
    for i, (id, hits) in enumerate(results):
        results[i] = (id, heapq.nlargest(k, hits, key = operator.itemgetter(1)))
    return results
