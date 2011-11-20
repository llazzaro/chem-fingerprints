from __future__ import absolute_import, with_statement
import unittest2

from support import fullpath

import chemfp
import chemfp.bitops
import _chemfp

CHEBI_TARGETS = fullpath("chebi_rdmaccs.fps")
CHEBI_QUERIES = fullpath("chebi_queries.fps.gz")

targets = chemfp.load_fingerprints(CHEBI_TARGETS, alignment=8)

alignment_methods = chemfp.bitops.get_alignment_methods()


all_methods = dict.fromkeys("LUT8-1 LUT8-4 LUT16-4 Lauradoux POPCNT".split())

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

all_alignments = dict.fromkeys("align1 align4 align8-small align8-large".split())

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
            method = chemfp.bitops.get_alignment_method(alignment)
            self.assertIn(method, all_methods)
            chemfp.bitops.set_alignment_method(alignment, "LUT8-1")
            self.assertEqual(chemfp.bitops.get_alignment_method(alignment), "LUT8-1")
            chemfp.bitops.set_alignment_method(alignment, method)
            self.assertEqual(chemfp.bitops.get_alignment_method(alignment), method)
            
        
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



class TestAlign8SmallMethods(unittest2.TestCase):
    def setUp(self):
        self.method = chemfp.bitops.get_alignment_method("align8-small")
    def tearDown(self):
        chemfp.bitops.set_alignment_method("align8-small", self.method)
        
    def _doit(self, method):
        chemfp.bitops.set_alignment_method("align8-small", method)
        self.assertEquals(chemfp.bitops.get_alignment_method("align8-small"), method)
        
        hits = targets.knearest_tanimoto_search_fp("00000000100410200290000b03a29241846163ee1f".decode("hex"), k=12, threshold=0.2)
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
            chemfp.bitops.set_alignment_method("align8-small", "Lauradoux")

    @unittest2.skipIf("POPCNT" not in alignment_methods, "CPU does not implement POPCNT")
    def test_popcnt(self):
        chemfp.bitops.set_alignment_method("align8-small", "POPCNT")
        self._doit()
        

class TestSelectFastestMethod(unittest2.TestCase):
    def setUp(self):
        self._alignment_methods = chemfp.bitops.get_alignment_methods()
    def tearDown(self):
        for k,v in self._alignment_methods.items():
            chemfp.bitops.set_alignment_method(k, v)
            
    def test_select_fastest(self):
        for alignment in all_alignments:
            chemfp.bitops.set_alignment_method(alignment, "LUT8-1")
            self.assertEquals(chemfp.bitops.get_alignment_method(alignment), "LUT8-1")

        chemfp.bitops.select_fastest_method()

        best_methods1 = chemfp.bitops.get_alignment_methods()

        for alignment in all_alignments:
            chemfp.bitops.set_alignment_method(alignment, "LUT8-1")
            self.assertEquals(chemfp.bitops.get_alignment_method(alignment), "LUT8-1")

        chemfp.bitops.select_fastest_method()

        best_methods2 = chemfp.bitops.get_alignment_methods()
        self.assertEquals(best_methods1, best_methods2)

        chemfp.bitops.select_fastest_method(repeat=-1000)

if __name__ == "__main__":
    unittest2.main()
    
