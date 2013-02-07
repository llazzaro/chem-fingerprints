# This tests check for memory leaks.

import unittest2

import resource
import sys
import array

import _chemfp
import chemfp

def get_memory():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

def memory_growth(f, *args, **kwargs):
    x = f(*args, **kwargs)
    m1 = get_memory()
    count = 0
    
    for i in xrange(1000):
        x = f(*args, **kwargs)
        m2 = get_memory()
        if m2 > m1:
            count += 1
            m1 = m2
    #print "count", count
    if count > 900:
        raise AssertionError("Memory growth in %r (%d/%d)!" % (f, count, 1000))
    return x

def fps_parse_id_fp():
    memory_growth(_chemfp.fps_parse_id_fp, 4, "ABCD\tSpam\n")


# I mistakenly used BuildValue("O") instead of BuildValue("N").
# This added an extra refcount to some of the returned strings,
# which meant they would never become freed.
class TestBuildValueMemoryLeaks(unittest2.TestCase):
    def test_make_unsorted_aligned_arena_when_aligned(self):
        s = "SpamAndEggs!"
        n1 = sys.getrefcount(s)
        memory_growth(_chemfp.make_unsorted_aligned_arena, s, 4)
        n2 = sys.getrefcount(s)
        if n2 > n1 + 10:
            raise AssertionError("Aligned growth has an extra refcount")

    def test_make_unsorted_aligned_arena_when_unaligned(self):
        SIZE = 1024*8
        s = "S" * SIZE
        n1 = sys.getrefcount(s)
        memory_growth(_chemfp.make_unsorted_aligned_arena, s, SIZE)
        n2 = sys.getrefcount(s)
        if n2 > n1 + 10:
            raise AssertionError("Aligned growth has an extra refcount: %d, %d" % (n1, n2))

    def test_align_fingerprint_same_size(self):
        s = "abcd"
        n1 = sys.getrefcount(s)
        memory_growth(_chemfp.align_fingerprint, s, 4, 4)
        n2 = sys.getrefcount(s)
        if n2 > n1 + 10:
            raise AssertionError("align_fingerprint with alignment has extra refcount: %d, %d" %
                                 (n1, n2))

    def test_align_fingerprint_needs_new_string(self):
        x = memory_growth(_chemfp.align_fingerprint, "ab", 16, 1024)
        # I don't know why the above doesn't find the leak.
        # I can probe the count directly.
        i = sys.getrefcount(x[2])
        if i != 2:
            raise AssertionError("Unexpected refcount: %d" % (i,))
    
    def test_make_sorted_aligned_arena_trivial(self):
        ordering = array.array("c", "\0\0\0\0"*16)  # space for a ChemFPOrderedPopcount (and extra)
        popcounts = array.array("c", "\0\0\0\0"*36)  # num_fingerprints + 2 (and extra)
        arena = "BLAH"
        n1 = sys.getrefcount(arena)
        memory_growth(_chemfp.make_sorted_aligned_arena, 32, 4, arena, 0, ordering, popcounts, 4)
        n2 = sys.getrefcount(arena)
        if n1 != n2:
            # This one doesn't need "N" because it borrowed the input refcount.
            raise AssertionError("Borrowed refcount should have been okay.")

    def test_make_sorted_aligned_arena_already_sorted_and_aligned(self):
        ordering = array.array("c", "\0\0\0\0"*16)  # space for a ChemFPOrderedPopcount (and extra)
        popcounts = array.array("c", "\0\0\0\0"*36)  # num_fingerprints + 2 (and extra)
        arena = "BDFL"
        n1 = sys.getrefcount(arena)
        memory_growth(_chemfp.make_sorted_aligned_arena, 32, 4, arena, 1, ordering, popcounts, 4)
        n2 = sys.getrefcount(arena)
        if n1 != n2:
            # This one doesn't need "N" because it borrowed the input refcount.
            raise AssertionError("Borrowed refcount should have been okay. %d != %d" % (n1, n2))

    def test_make_sorted_aligned_arena_already_sorted_but_not_aligned(self):
        ordering = array.array("c", "\0\0\0\0"*16)  # space for a ChemFPOrderedPopcount (and extra)
        popcounts = array.array("c", "\0\0\0\0"*40)  # num_fingerprints + 2 (and extra)
        arena = "QWOP" + "\0" * (4096-4)
        n1 = sys.getrefcount(arena)
        memory_growth(_chemfp.make_sorted_aligned_arena, 32, 4096, arena, 1, ordering, popcounts, 1024)
        n2 = sys.getrefcount(arena)
        if n1 != n2:
            # This one shouldn't touch the input arena
            raise AssertionError("Borrowed refcount should have been okay. %d != %d" % (n1, n2))

    def test_make_sorted_aligned_arena_when_not_sorted(self):
        ordering = array.array("c", "\0\0\0\0"*16)  # space for a ChemFPOrderedPopcount (and extra)
        popcounts = array.array("c", "\0\0\0\0"*40)  # num_fingerprints + 2 (and extra)
        arena = "QWON" + "\0" * (4096-4)
        n1 = sys.getrefcount(arena)
        memory_growth(_chemfp.make_sorted_aligned_arena, 32, 4, arena, 2, ordering, popcounts, 1024)
        n2 = sys.getrefcount(arena)
        if n1 != n2:
            # This one shouldn't touch the input arena
            raise AssertionError("Borrowed refcount should have been okay. %d != %d" % (n1, n2))

        
    def test_arena_copy(self):
        data = "ABCD\t0\n" * 200
        from cStringIO import StringIO
        arena = chemfp.load_fingerprints(StringIO(data))
        def make_subarena_copy():
            arena[1:].copy()
        memory_growth(make_subarena_copy)
    
if __name__ == "__main__":
    unittest2.main()
