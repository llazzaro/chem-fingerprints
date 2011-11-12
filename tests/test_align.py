from cStringIO import StringIO

import unittest2

import chemfp

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

if __name__ == "__main__":
    unittest2.main()
