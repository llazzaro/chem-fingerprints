"""This is an internal module. Do not call these functions directly"""
import operator
import heapq

from chemfp import bitops

# This code is so much simpler than the C extension code.
# And about 6 times slower, even with bitops.byte_tanimoto in C!

def generic_tanimoto(query, targets, threshold):
    for fp, id in targets:
        score = bitops.byte_tanimoto(query, fp)
        if score == -1:
            raise AssertionError("should not happen")
        if score >= threshold:
            yield id, score
    
def generic_tanimoto_count(query, targets, threshold):
    return sum(1 for (fp, id) in targets if bitops.byte_tanimoto(query, fp) >= threshold)

def generic_tanimoto_knearest(query, targets, k, threshold):
    results = heapq.nlargest(k, generic_tanimoto(query, targets, threshold),
                             key=operator.itemgetter(1))
    return results
