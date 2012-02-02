from __future__ import absolute_import, with_statement
import unittest2

from support import fullpath

import chemfp
import chemfp.bitops
import _chemfp

set_alignment_method = chemfp.bitops.set_alignment_method
get_alignment_method = chemfp.bitops.get_alignment_method

CHEBI_TARGETS = fullpath("chebi_rdmaccs.fps")
CHEBI_QUERIES = fullpath("chebi_queries.fps.gz")

targets = chemfp.load_fingerprints(CHEBI_TARGETS, alignment=8)
targets_64 = chemfp.load_fingerprints(CHEBI_TARGETS, alignment=64)

available_methods = chemfp.bitops.get_methods()
alignment_methods = chemfp.bitops.get_alignment_methods()


all_methods = dict.fromkeys("LUT8-1 LUT8-4 LUT16-4 Lauradoux POPCNT Gillies ssse3".split())

class TestMethods(unittest2.TestCase):
    def test_no_duplicates(self):
        methods = chemfp.bitops.get_methods()
        self.assertEquals(len(methods), len(set(methods)))

    def test_for_unknown_methods(self):
        for method in chemfp.bitops.get_methods():
            self.assertIn(method, all_methods)

    def test_for_possible_missing_popcnt(self):
        if len(all_methods) == 4:
            self.assertNotIn("POPCNT", chemfp.get_methods())
        
    def test_internal_bad_args(self):
        with self.assertRaisesRegexp(IndexError, "method index is out of range"):
            _chemfp.get_method_name(-1)
        with self.assertRaisesRegexp(IndexError, "method index is out of range"):
            _chemfp.get_method_name(_chemfp.get_num_methods())

all_alignments = dict.fromkeys("align1 align4 align8-small align8-large align-ssse3".split())

class TestAlignments(unittest2.TestCase):
    def test_no_duplicates(self):
        alignments = chemfp.bitops.get_alignments()
        self.assertEquals(len(alignments), len(set(alignments)))
        self.assertEquals(len(alignments), len(all_alignments))

    def test_for_unknown_alignments(self):
        for alignment in chemfp.bitops.get_alignments():
            self.assertIn(alignment, all_alignments)

    def test_get_set_alignment_method(self):
        for alignment in chemfp.bitops.get_alignments():
            method = get_alignment_method(alignment)
            self.assertIn(method, all_methods)
            set_alignment_method(alignment, "LUT8-1")
            self.assertEqual(get_alignment_method(alignment), "LUT8-1")
            set_alignment_method(alignment, method)
            self.assertEqual(get_alignment_method(alignment), method)
            
        
    def test_internal_bad_args(self):
        with self.assertRaisesRegexp(IndexError, "alignment index is out of range"):
            _chemfp.get_alignment_name(-1)
        with self.assertRaisesRegexp(IndexError, "alignment index is out of range"):
            _chemfp.get_alignment_name(_chemfp.get_num_methods())


        # I didn't want a better error code for this
        with self.assertRaisesRegexp(ValueError, "Bad argument"):
            _chemfp.get_alignment_name(_chemfp.get_alignment_method(-1))
        with self.assertRaisesRegexp(ValueError, "Bad argument"):
            _chemfp.get_alignment_name(_chemfp.get_alignment_method(100))

        with self.assertRaisesRegexp(ValueError, "Bad argument"):
            _chemfp.get_alignment_name(_chemfp.set_alignment_method(-1, 0))
        with self.assertRaisesRegexp(ValueError, "Bad argument"):
            _chemfp.get_alignment_name(_chemfp.set_alignment_method(100, 0))


    def test_cannot_use_64_bit_method_for_shorter_bit_alignment(self):
        msg = "Mismatch between popcount method and alignment type"
        available_methods = chemfp.bitops.get_methods()
        for method in ("Lauradoux", "Gillies", "POPCNT"):
            if (method == "POPCNT") and ("POPCNT" not in available_methods):
                continue
            with self.assertRaisesRegexp(ValueError, msg):
                set_alignment_method("align1", method)
            with self.assertRaisesRegexp(ValueError, msg):
                set_alignment_method("align4", method)


    @unittest2.skipIf("ssse3" not in available_methods, "CPU does not implement SSSE3")
    def test_ssse3(self):
        method = get_alignment_method("align-ssse3")
        # This disables SSSE3 support
        set_alignment_method("align-ssse3", "LUT8-1")
        self.assertEquals(get_alignment_method("align-ssse3"), "LUT8-1")
        set_alignment_method("align-ssse3", "ssse3")
        self.assertEquals(get_alignment_method("align-ssse3"), "ssse3")
        set_alignment_method("align-ssse3", method)



class TestAlign8SmallMethods(unittest2.TestCase):
    def setUp(self):
        self.small_method = get_alignment_method("align8-small")
        self.large_method = get_alignment_method("align8-large")
        self.ssse3_method = get_alignment_method("align-ssse3")
    def tearDown(self):
        set_alignment_method("align8-small", self.small_method)
        set_alignment_method("align8-large", self.large_method)
        set_alignment_method("align-ssse3", self.ssse3_method)
        
    def _doit(self, method):
        for alignment in ("align8-small", "align8-large", "align-ssse3"):
            set_alignment_method(alignment, method)
            self.assertEquals(get_alignment_method(alignment), method)
        
        hits = targets.knearest_tanimoto_search_fp("00000000100410200290000b03a29241846163ee1f".decode("hex"), k=12, threshold=0.2).get_ids_and_scores()
        self.assertEqual(hits, [('CHEBI:8069', 1.0),
                                ('CHEBI:6758', 0.78723404255319152),
                                ('CHEBI:7983', 0.73999999999999999),
                                ('CHEBI:8107', 0.6956521739130435),
                                ('CHEBI:17568', 0.6904761904761905),
                                ('CHEBI:16294', 0.6818181818181818),
                                ('CHEBI:16964', 0.673469387755102),
                                ('CHEBI:17477', 0.6458333333333334),
                                ('CHEBI:17025', 0.62),
                                ('CHEBI:15901', 0.6122448979591837),
                                ('CHEBI:16742', 0.6122448979591837),
                                ('CHEBI:4888', 0.6078431372549019)])

    def test_lut8_1(self):
        self._doit("LUT8-1")
        
    def test_lut8_4(self):
        self._doit("LUT8-4")
        
    def test_lut16_4(self):
        self._doit("LUT16-4")

    def test_lauradoux(self):
        with self.assertRaisesRegexp(ValueError, "Mismatch between popcount method and alignment type"):
            set_alignment_method("align8-small", "Lauradoux")

    @unittest2.skipIf("POPCNT" not in alignment_methods, "CPU does not implement POPCNT")
    def test_popcnt(self):
        self._doit("POPCNT")

class TestAlign8LargeMethods(unittest2.TestCase):
    def setUp(self):
        self.large_method = get_alignment_method("align8-large")
        self.ssse3_method = get_alignment_method("align-ssse3")
    def tearDown(self):
        set_alignment_method("align8-large", self.large_method)
        set_alignment_method("align-ssse3", self.ssse3_method)
        
    def _doit(self, method, use_ssse3=False):
        set_alignment_method("align8-large", method)
        self.assertEquals(get_alignment_method("align8-large"), method)
        if use_ssse3:
            set_alignment_method("align-ssse3", "ssse3")
            self.assertEquals(get_alignment_method("align-ssse3"), "ssse3")
        else:
            set_alignment_method("align-ssse3", "LUT8-1")
            self.assertEquals(get_alignment_method("align-ssse3"), "LUT8-1")
        
        hits = targets_64.knearest_tanimoto_search_fp("00000000100410200290000b03a29241846163ee1f".decode("hex"), k=12, threshold=0.2).get_ids_and_scores()
        self.assertEqual(hits, [('CHEBI:8069', 1.0),
                                ('CHEBI:6758', 0.78723404255319152),
                                ('CHEBI:7983', 0.73999999999999999),
                                ('CHEBI:8107', 0.6956521739130435),
                                ('CHEBI:17568', 0.6904761904761905),
                                ('CHEBI:16294', 0.6818181818181818),
                                ('CHEBI:16964', 0.673469387755102),
                                ('CHEBI:17477', 0.6458333333333334),
                                ('CHEBI:17025', 0.62),
                                ('CHEBI:15901', 0.6122448979591837),
                                ('CHEBI:16742', 0.6122448979591837),
                                ('CHEBI:4888', 0.6078431372549019)])

    def test_lut8_1(self):
        self._doit("LUT8-1")
        
    def test_lut8_4(self):
        self._doit("LUT8-4")
        
    def test_lut16_4(self):
        self._doit("LUT16-4")
        
    def test_lauradoux(self):
        self._doit("Lauradoux")
        
    def test_gillies(self):
        self._doit("Gillies")

    @unittest2.skipIf("ssse3" not in available_methods, "CPU does not implement SSSE3")
    def test_ssse3(self):
        self._doit("Lauradoux", use_ssse3=True)

    @unittest2.skipIf("POPCNT" not in available_methods, "CPU does not implement POPCNT")
    def test_popcnt(self):
        self._doit("POPCNT")
        

class TestSelectFastestMethod(unittest2.TestCase):
    def setUp(self):
        self._alignment_methods = chemfp.bitops.get_alignment_methods()
    def tearDown(self):
        for k,v in self._alignment_methods.items():
            set_alignment_method(k, v)
            
    def test_select_fastest(self):
        for alignment in all_alignments:
            set_alignment_method(alignment, "LUT8-1")
            self.assertEquals(get_alignment_method(alignment), "LUT8-1")

        chemfp.bitops.select_fastest_method()

        best_methods1 = chemfp.bitops.get_alignment_methods()

        for alignment in all_alignments:
            set_alignment_method(alignment, "LUT8-1")
            self.assertEquals(get_alignment_method(alignment), "LUT8-1")

        chemfp.bitops.select_fastest_method()

        best_methods2 = chemfp.bitops.get_alignment_methods()
        self.assertEquals(best_methods1, best_methods2) # This might fail if two methods have nearly identical timings

        chemfp.bitops.select_fastest_method(repeat=-1000)

if __name__ == "__main__":
    chemfp.bitops.use_environment_variables()
    unittest2.main()
    
