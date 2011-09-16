from __future__ import absolute_import
import unittest2

import chemfp

from support import fullpath

CHEBI_TARGETS = fullpath("chebi_rdmaccs.fps")
CHEBI_QUERIES = fullpath("chebi_queries.fps.gz")

# Backwards compatibility for Python 2.5
try:
    next
except NameError:
    def next(it):
        return it.next()
    

class CommonReaderAPI(object):
    _open = None
    
    def _check_target_metadata(self, metadata):
        self.assertEquals(metadata.num_bits, 166)
        self.assertEquals(metadata.num_bytes, 21)
        self.assertEquals(metadata.software, "OEChem/1.7.4 (20100809)")
        self.assertEquals(metadata.type, "RDMACCS-OpenEye/1")
        self.assertEquals(metadata.sources, ["/Users/dalke/databases/ChEBI_lite.sdf.gz"])
        self.assertEquals(metadata.date, "2011-09-16T13:49:04")
        self.assertEquals(metadata.aromaticity, "mmff")

    def _check_query_metadata(self, metadata):
        self.assertEquals(metadata.num_bits, 166)
        self.assertEquals(metadata.num_bytes, 21)
        self.assertEquals(metadata.software, "OEChem/1.7.4 (20100809)")
        self.assertEquals(metadata.type, "RDMACCS-OpenEye/1")
        self.assertEquals(metadata.sources, ["/Users/dalke/databases/ChEBI_lite.sdf.gz"])
        self.assertEquals(metadata.date, "2011-09-16T13:28:43")
        self.assertEquals(metadata.aromaticity, "openeye")
        
    
    def test_uncompressed_open(self):
        reader = self._open(CHEBI_TARGETS)
        self._check_target_metadata(reader.metadata)
        num = sum(1 for x in reader)
        self.assertEquals(num, 2000)

    def test_compressed_open(self):
        reader = self._open(CHEBI_QUERIES)
        self._check_query_metadata(reader.metadata)
        num = sum(1 for x in reader)
        self.assertEquals(num, 154)

    def test_iteration(self):
        reader = iter(self._open(CHEBI_TARGETS))
        fields = [next(reader) for i in range(5)]
        self.assertEquals(fields, 
                          [("CHEBI:776", "00000000000000008200008490892dc00dc4a7d21e".decode("hex")),
                           ("CHEBI:1148", "000000000000200080000002800002040c0482d608".decode("hex")),
                           ("CHEBI:1734", "0000000000000221000800111601017000c1a3d21e".decode("hex")),
                           ("CHEBI:1895", "00000000000000000000020000100000000400951e".decode("hex")),
                           ("CHEBI:2303", "0000000002001021820a00011681015004cdb3d21e".decode("hex"))])

      
    def test_iter_arenas_default_sizes(self):
        reader = self._open(CHEBI_TARGETS)
        count = 0
        for arena in reader.iter_arenas():
            self._check_target_metadata(arena.metadata)
            if count == 0:
                self.assertEquals(arena.ids[-5:],
                                  ['CHEBI:16316', 'CHEBI:16317', 'CHEBI:16318', 'CHEBI:16319', 'CHEBI:16320'])
                
            self.assertEquals(len(arena), 1000)  # There should be two of these
            count += 1
        self.assertEquals(count, 2)
        self.assertEquals(arena.ids[-5:],
                          ['CHEBI:17578', 'CHEBI:17579', 'CHEBI:17580', 'CHEBI:17581', 'CHEBI:17582'])

    def test_iter_arenas_select_size(self):
        reader = self._open(CHEBI_TARGETS)
        count = 0
        for arena in reader.iter_arenas(100):
            self._check_target_metadata(arena.metadata)
            if count == 0:
                self.assertEquals(arena.ids[-5:],
                                  ['CHEBI:5280', 'CHEBI:5445', 'CHEBI:5706', 'CHEBI:5722', 'CHEBI:5864'])
            self.assertEquals(len(arena), 100)  # There should be two of these
            count += 1
        self.assertEquals(count, 20)
        self.assertEquals(arena.ids[:5],
                          ['CHEBI:17457', 'CHEBI:17458', 'CHEBI:17459', 'CHEBI:17460', 'CHEBI:17464'])
    
class TestFPSReader(unittest2.TestCase, CommonReaderAPI):
    _open = staticmethod(chemfp.open)

    def test_row_iteration(self):
        reader = chemfp.open(CHEBI_TARGETS)
        num = sum(1 for x in reader.iter_rows())
        self.assertEquals(num, 2000)
        
        row_reader = chemfp.open(CHEBI_TARGETS).iter_rows()
        fields = [next(row_reader) for i in range(5)]
        self.assertEquals(fields,  [
            ['00000000000000008200008490892dc00dc4a7d21e', 'CHEBI:776'],
            ['000000000000200080000002800002040c0482d608', 'CHEBI:1148'],
            ['0000000000000221000800111601017000c1a3d21e', 'CHEBI:1734'],
            ['00000000000000000000020000100000000400951e', 'CHEBI:1895'],
            ['0000000002001021820a00011681015004cdb3d21e', 'CHEBI:2303']])

    def test_iter_blocks(self):
        reader = chemfp.open(CHEBI_TARGETS)
        line_counts = 0
        has_776 = False
        has_17582 = False
        for block in reader.iter_blocks():
            line_counts += block.count("\n")
            if "00000000000000008200008490892dc00dc4a7d21e\tCHEBI:776" in block:
                has_776 = True
            if "00000000020012008008000104000064844ca2521c\tCHEBI:17582" in block:
                has_17582 = True

        self.assertEquals(line_counts, 2000)
        self.assertTrue(has_776, "Missing CHEBI:776")
        self.assertTrue(has_17582, "Missing CHEBI:17582")

    #
    # Count tanimoto hits using a fingerprint
    # 

    def test_count_tanimoto_hits_fp_default(self):
        reader = self._open(CHEBI_TARGETS)
        num_hits = reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"))
        self.assertEqual(num_hits, 176)

    def test_count_tanimoto_hits_fp_set_default(self):
        # This is set to the default value
        reader = self._open(CHEBI_TARGETS)
        num_hits = reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"),
                                                 threshold = 0.7)
        self.assertEqual(num_hits, 176)

    def test_count_tanimoto_hits_fp_set_threshold(self):
        reader = self._open(CHEBI_TARGETS)
        num_hits = reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"),
                                                 threshold = 0.8)
        self.assertEqual(num_hits, 108)

    def test_count_tanimoto_hits_fp_set_max_threshold(self):
        reader = self._open(CHEBI_TARGETS)
        num_hits = reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"),
                                                 threshold = 1.0)
        self.assertEqual(num_hits, 1)

    def test_count_tanimoto_hits_fp_threshold_range_error(self):
        reader = self._open(CHEBI_TARGETS)
        with self.assertRaisesRegexp(ValueError, "threshold must between 0.0 and 1.0, inclusive") as e:
            reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"),
                                          threshold = 1.1)
        reader = self._open(CHEBI_TARGETS)
        with self.assertRaisesRegexp(ValueError, "threshold must between 0.0 and 1.0, inclusive") as e:
            reader.count_tanimoto_hits_fp("000000102084322193de9fcfbffbbcfbdf7ffeff1f".decode("hex"),
                                          threshold = -0.00001)

    #
    # Count tanimoto hits using an arena
    #

    def test_count_tanimoto_arena_default(self):
        targets = self._open(CHEBI_TARGETS)
        query_arena = next(chemfp.open(CHEBI_QUERIES).iter_arenas(10))
        hits = targets.count_tanimoto_hits_arena(query_arena)
        self.assertEquals(hits, [4, 179, 40, 32, 1, 3, 28, 11, 46, 7])

    def test_count_tanimoto_arena_set_default(self):
        targets = self._open(CHEBI_TARGETS)
        query_arena = next(chemfp.open(CHEBI_QUERIES).iter_arenas(10))
        hits = targets.count_tanimoto_hits_arena(query_arena, threshold=0.7)
        self.assertEquals(hits, [4, 179, 40, 32, 1, 3, 28, 11, 46, 7])

    def test_count_tanimoto_arena_set_threshold(self):
        targets = self._open(CHEBI_TARGETS)
        query_arena = next(chemfp.open(CHEBI_QUERIES).iter_arenas(10))
        hits = targets.count_tanimoto_hits_arena(query_arena, threshold=0.9)
        self.assertEquals(hits, [0, 97, 7, 1, 0, 1, 1, 0, 1, 1])

    def test_count_tanimoto_hits_arena_threshold_range_error(self):
        reader = self._open(CHEBI_TARGETS)
        query_arena = next(chemfp.open(CHEBI_QUERIES).iter_arenas(10))
        with self.assertRaisesRegexp(ValueError, "threshold must between 0.0 and 1.0, inclusive") as e:
            reader.count_tanimoto_hits_arena(query_arena, threshold = 1.1)
                                          
        reader = self._open(CHEBI_TARGETS)
        with self.assertRaisesRegexp(ValueError, "threshold must between 0.0 and 1.0, inclusive") as e:
            reader.count_tanimoto_hits_arena(query_arena, threshold = -0.00001)
        
        
        


class TestLoadFingerprints(unittest2.TestCase, CommonReaderAPI):
    # Hook to handle the common API
    def _open(self, name):
        return chemfp.load_fingerprints(name, reorder=False)
    
            

if __name__ == "__main__":
    unittest2.main()

