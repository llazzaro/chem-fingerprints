from cStringIO import StringIO

import unittest2

import chemfp
from chemfp import arena, bitops

zeros = ("0000\tfirst\n"
         "0010\tsecond\n"
         "0000\tthird\n")

ordered_zeros = ("0000\tfirst\n"
                 "0000\tthird\n"
                 "0010\tsecond\n")

class TestUnsortedAlignment(unittest2.TestCase):

    # Python strings are 1-, 2-, and 4-byte aligned
    
    def test_1_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=False, alignment=1)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x10\x00\x00")
        self.assertEquals(a.storage_size, 2)
    def test_2_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=False, alignment=2)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x10\x00\x00")
        self.assertEquals(a.storage_size, 2)
    def test_4_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=False, alignment=4)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00")
        self.assertEquals(a.storage_size, 4)

    # A Python string might be 8-aligned, but no guarantee
    def test_8_alignment(self):
        arenas = [chemfp.load_fingerprints(StringIO(zeros), reorder=False, alignment=8)
                  for i in range(10)]
        for a in arenas:
            if a.start_padding == a.end_padding == 0:
                s = a.arena
            else:
                self.assertEquals(a.arena[:a.start_padding], "\x00" * a.start_padding)
                self.assertEquals(a.arena[-a.end_padding:], "\x00" * a.end_padding)
                s = a.arena[a.start_padding:-a.end_padding]
            
            self.assertEquals(s,
                              "\x00\x00\x00\x00\x00\x00\x00\x00"
                              "\x00\x10\x00\x00\x00\x00\x00\x00"
                              "\x00\x00\x00\x00\x00\x00\x00\x00")
            self.assertEquals(a.storage_size, 8)

class TestReorderedAlignment(unittest2.TestCase):
    def test_1_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=True, alignment=1)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x10")

    def test_2_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=True, alignment=2)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 1)
        # The code overallocates one byte
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x10\x00")

    def test_4_alignment(self):
        a = chemfp.load_fingerprints(StringIO(zeros), reorder=True, alignment=4)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 3)
        self.assertEquals(a.arena,
                          "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00")


class TestAlreadyOrderedAlignment(unittest2.TestCase):
    def test_1_alignment(self):
        a = chemfp.load_fingerprints(StringIO(ordered_zeros), reorder=True, alignment=1)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x10")

    def test_2_alignment(self):
        a = chemfp.load_fingerprints(StringIO(ordered_zeros), reorder=True, alignment=2)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x10")

    def test_4_alignment(self):
        a = chemfp.load_fingerprints(StringIO(ordered_zeros), reorder=True, alignment=4)
        self.assertEquals(a.start_padding, 0)
        self.assertEquals(a.end_padding, 0)
        self.assertEquals(a.arena, "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00")

    def test_16_alignment(self):
        arenas = [chemfp.load_fingerprints(StringIO(ordered_zeros), reorder=True, alignment=16)
                  for i in range(10)]
        for a in arenas:
            if a.start_padding == 0 and a.end_padding == 0:
                s = a.arena
            else:
                self.assertEquals(a.start_padding + a.end_padding + 1, 16)
                self.assertEquals(a.arena[:a.start_padding], "\x00" * a.start_padding)
                self.assertEquals(a.arena[-a.end_padding:], "\x00" * a.end_padding)
                s = a.arena[a.start_padding:-a.end_padding]

            self.assertEquals(s, ("\x00"*16 +
                                  "\x00"*16 +
                                  "\x00\x10" + "\x00"*14))

class TestOptimalAlignment(unittest2.TestCase):
    def setUp(self):
        self._has_popcnt = arena._has_popcnt
        self._has_ssse3 = arena._has_ssse3
        self.get_alignment_method = bitops.get_alignment_method
        bitops.get_alignment_method = self.hook

    def tearDown(self):
        arena._has_popcnt = self._has_popcnt
        arena._has_ssse3 = self._has_ssse3
        bitops.get_alignment_method = self.get_alignment_method

    def hook(self, name):
        return self.data[name]

    def test_tiny(self):
        for i in range(1, 9):
            self.assertEquals(arena.get_optimal_alignment(i), 1)

    def test_small(self):
        for i in range(9, 33):
            self.assertEquals(arena.get_optimal_alignment(i), 4)


    def test_medium(self):
        arena._has_popcnt = False
        arena._has_ssse3 = False
        for i in range(35, 225):
            self.assertEquals(arena.get_optimal_alignment(i), 8)

    def test_popcnt_ssse3_combinations(self):
        arena._has_popcnt = True
        arena._has_ssse3 = True
        
        self.data = {"align8-large": "POPCNT",
                     "align8-small": "POPCNT",
                     "align-ssse3": "ssse3"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 8)
        
        self.data = {"align8-large": "POPCNT",
                     "align8-small": "LUT8-1",
                     "align-ssse3": "ssse3"}
        self.assertEquals(arena.get_optimal_alignment(300), 64)
        self.assertEquals(arena.get_optimal_alignment(800), 8)

        self.data = {"align8-large": "LUT8-1",
                     "align8-small": "POPCNT",
                     "align-ssse3": "ssse3"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 64)

        self.data = {"align8-large": "LUT8-1",
                     "align8-small": "LUT8-1",
                     "align-ssse3": "ssse3"}
        self.assertEquals(arena.get_optimal_alignment(300), 64)
        self.assertEquals(arena.get_optimal_alignment(800), 64)

        # And now, with SSSE3 disabled

        self.data = {"align8-large": "POPCNT",
                     "align8-small": "POPCNT",
                     "align-ssse3": "POPCNT"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 8)
        
        self.data = {"align8-large": "POPCNT",
                     "align8-small": "LUT8-1",
                     "align-ssse3": "POPCNT"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 8)

        self.data = {"align8-large": "LUT8-1",
                     "align8-small": "POPCNT",
                     "align-ssse3": "POPCNT"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 8)

        self.data = {"align8-large": "LUT8-1",
                     "align8-small": "LUT8-1",
                     "align-ssse3": "POPCNT"}
        self.assertEquals(arena.get_optimal_alignment(300), 8)
        self.assertEquals(arena.get_optimal_alignment(800), 8)


        
        

if __name__ == "__main__":
    unittest2.main()
