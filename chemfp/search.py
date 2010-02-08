import operator
import heapq

from chemfp import bitops

def tanimoto_search(query_fp, targets, threshold=0.0, k=None):
    if not (0.0 <= threshold <= 1.0):
        raise TypeError("threshold must be between 0.0 and 1.1, inclusive")
    if not (k is None or k > 0):
        raise TypeError("k must be None or a positive value")

    f = getattr(targets, "tanimoto_search", None)
    if f is not None:
        return f(query_fp, threshold, k)

    preferred_form = getattr(targets, "preferred_form", "byte")
    if preferred_form == "hex":
        tsearch = _hex_tanimoto_search(query_fp.encode("hex"), targets, threshold)
    else:
        tsearch = _byte_tanimoto_serach(query_fp, targets, threshold)

    if k is not None:
        tsearch = heapq.nlargest(k, tsearch, key=operator.itemgetter(1))
    return iter(tsearch)
        
# This exists to make searching FPS files faster
def _hex_tanimoto_search(hex_query_fp, targets, threshold):
    for rec in targets.hex_iter():
        score = bitops.hex_tanimoto(hex_query_fp, rec.hex_fp)
        if score == -1:
            # only happens when one of the fingerprints wasn't hex
            # but by construction, that's not possible here
            raise AssertionError("there should be no way to return -1 here")
        if score >= threshold:
            yield rec.id, score

def _byte_tanimoto_search(query_fp, targets, threshold):
    for rec in targets:
        score = bitopts.byte_tanimoto(query_fp, rec.byte_fp)
        if score == -1:
            raise AssertionError("there should be no way to return -1 here")
        if score >= threshold:
            yield rec.id, score


HEX = "0022010310401020310004000000000800000010000000000000000000320100180000022303000000060d300001004009024000500000000102000010000004000040000000008001200c0012f0010e0030075000002000e041020a0000200008200020040020000000000000000000400c3004004014200034000100000100"
FP = HEX.decode("hex")

def main():
    from chemfp import fps_reader
    reader = fps_reader.open("nci.fps")
    for id, score in tanimoto_search(FP, reader, k=10):
        print id, score

if __name__ == "__main__":
    main()
