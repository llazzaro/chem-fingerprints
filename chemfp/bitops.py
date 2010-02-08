import array
import heapq
import functools
import operator

from _chemfp import (hex_isvalid, hex_popcount, hex_intersect_popcount,
                     hex_tanimoto, hex_contains)

from _chemfp import (byte_popcount, byte_intersect_popcount,
                     byte_tanimoto, byte_contains)

from _chemfp import nlargest_tanimoto_block as _nlargest_tanimoto_block

def _score_iterator(calculate_score, target_fps):
    for i, (target_fp, target_id) in enumerate(target_fps):
        print "Calc for", repr(target_fp), target_fps
        yield calculate_score(target_fp), -i, target_id


def _nlargest_tanimoto(n, calculate_score, target_fps):
    result = heapq.nlargest(n, _score_iterator(calculate_score, target_fps))
    return [(-i, target_id, target_score) for (target_score, i, target_id) in result]

import itertools
c = itertools.count()
def hex_nlargest_tanimoto(n, query_fp, target_fps):
    def calculate_hex_score(target_fp, hex_tanimoto=hex_tanimoto, query_fp=query_fp):
        print repr(query_fp)
        print repr(target_fp)
        score = hex_tanimoto(query_fp, target_fp)
        if score == -1.0:
            raise TypeError("Fingerprint index %i is not a hex fingerprint" % (index,))
        return score
    return _nlargest_tanimoto(n, calculate_hex_score, target_fps)

def byte_nlargest_tanimoto(n, query_fp, target_fps):
    return _nlargest_tanimoto(n, functools.partial(byte_tanimoto, query_fp), target_fps)


def byte_nlargest_tanimoto_block(n, query_fp, target_block, offset=0,
                                 storage_len=None, threshold=0.0):
    if storage_len is None:
        storage_len = len(query_fp)
    indicies = array.array('i', (0,)*n)
    scores = array.array('d', (0.0,)*n)
    k = _nlargest_tanimoto_block(n, query_fp, target_block, offset,
                                 storage_len, threshold, indicies, scores)
    return [(indicies[i], scores[i]) for i in range(k)]

