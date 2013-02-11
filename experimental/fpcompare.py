# Compare multiple fingerprints

# fpcompare a.fps b.fps c.fps
##fpcounts
# Says that bit 120 differs, with a=0, b=0, c=1
#1 120:4 <name>

import chemfp
import itertools
import sys
from collections import defaultdict

def different_bits(reference_id, rows):
    reference_fp, reference_id = rows[0]
    for fp, fp_id in rows[1:]:
        if fp_id != reference_id:
            raise AssertionError
        if fp != reference_fp:
            break
    else:
        # They are all the same. Simple to skip
        return []

    # Something is different.
    bits = []
    terms = []
    for byteno, chars in enumerate(itertools.izip(*[fp for (fp, id) in rows])):
        if set(chars) == 1:
            continue
        words = map(ord, chars)
        for bitno in range(8):
            b = 0
            for i, w in enumerate(words):
                if w & (1<<bitno):
                    b += (1<<i)
            if b not in (0, 2**len(words)-1):
                terms.append( (byteno*8+bitno, b) )
                
    assert len(terms)
    return terms
    

filenames = sys.argv[1:]
readers = map(chemfp.open, filenames)
mismatches = defaultdict(int)

# Make sure the num_bits are the same
for rows in itertools.izip(*readers):
    reference_fp, reference_id = rows[0]
    terms = different_bits(reference_id, rows)
    for bitno, pattern in terms:
        mismatches[bitno] += 1
    print " ".join([str(len(terms))] + ["%s:%s" % term for term in terms] + [reference_id])
for k, v in sorted(mismatches.items()):
    print k, v
