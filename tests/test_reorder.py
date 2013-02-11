# Test the fingerprint reordering implementation

# Note: this is the ordering by popcount, and NOT the ordering of the results!

from __future__ import absolute_import, with_statement
import unittest2
from cStringIO import StringIO
import struct

import chemfp
from chemfp.bitops import byte_popcount

def _load(fingerprints, reorder):
    if len(fingerprints) == 0:
        num_bits = 16
    else:
        num_bits = len(fingerprints[0])*8
    id_fps = ((str(i), fp) for (i, fp) in enumerate(fingerprints))
    return chemfp.load_fingerprints(id_fps,
                                    metadata=chemfp.Metadata(num_bits=num_bits),
                                    reorder=reorder, alignment=1)

def verify_popcount_indices(arena):
    assert len(arena.popcount_indices) % 4 == 0
    format = "i"*(len(arena.popcount_indices)//4)
    values = struct.unpack(format, arena.popcount_indices)
    assert len(values) == arena.metadata.num_bits+2, (len(values), arena.metadata.num_bits)
    assert values[-1] == len(arena), (values[-1], len(arena), values)
    for i in range(len(values)-1):
        start = values[i]
        end = values[i+1]

        for j in range(start, end):
            fp = arena[j][1]
            assert byte_popcount(fp) == i, (byte_popcount(fp), i, start, end)
    
    

class TestReorder(unittest2.TestCase):
    def test_empty(self):
        arena = _load([], True)
        self.assertEquals(arena.arena, "")
        verify_popcount_indices(arena)
        
        arena = _load([], False)
        self.assertEquals(arena.arena, "")
        self.assertEquals(arena.popcount_indices, "")

    def test_single(self):
        arena = _load(["1234"], True)
        self.assertEquals(arena.arena, "1234")
        verify_popcount_indices(arena)
        
        arena = _load(["1234"], False)
        self.assertEquals(arena.arena, "1234")

    def test_two(self):
        arena = _load(["AA", "CC"], True)
        self.assertEquals(arena.arena, "AACC")
        verify_popcount_indices(arena)
        
        arena = _load(["AA", "CC"], False)
        self.assertEquals(arena.arena, "AACC")
        
        arena = _load(["CC", "AA"], True)
        self.assertEquals(arena.arena, "AACC")
        verify_popcount_indices(arena)
        
        arena = _load(["CC", "AA"], False)
        self.assertEquals(arena.arena, "CCAA")

    def test_every_bit_unsorted(self):
        arena = _load([chr(i) for i in range(256)], False)
        self.assertEquals(arena.arena, "".join(chr(i) for i in range(256)))

    def test_every_bit_sorted(self):
        # This is a bit tricker since there's no guaranteed order
        # of the contents of the arena
        
        all_bytes = [chr(i) for i in range(256)]
        expected = sorted(byte_popcount(fp) for fp in all_bytes)

        arena = _load(all_bytes, True)
        popcounts = map(byte_popcount, arena.arena)
        self.assertEquals(popcounts, expected)
        verify_popcount_indices(arena)

        arena = _load(all_bytes[::-1], True)
        popcounts = map(byte_popcount, arena.arena)
        self.assertEquals(popcounts, expected)
        verify_popcount_indices(arena)

    def test_every_bit_sorted_tripled(self):
        # This is a bit tricker since there's no guaranteed order
        # of the contents of the arena
        
        all_bytes = [chr(i) for i in range(256)]
        all_bytes *= 3
        expected = sorted(byte_popcount(fp) for fp in all_bytes)

        arena = _load(all_bytes, True)
        popcounts = map(byte_popcount, arena.arena)
        self.assertEquals(popcounts, expected)
        verify_popcount_indices(arena)

        arena = _load(all_bytes[::-1], True)
        popcounts = map(byte_popcount, arena.arena)
        self.assertEquals(popcounts, expected)
        verify_popcount_indices(arena)


    def test_all_ones(self):
        arena = _load([chr(255), chr(255), chr(255), chr(255)], True)
        self.assertEquals(arena.arena, "\xff\xff\xff\xff")
        verify_popcount_indices(arena)

if __name__ == "__main__":
    unittest2.main()
