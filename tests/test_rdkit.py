# Test the layer directly on top of RDKit

import unittest2
from chemfp import rdkit, decompressors

class TestNormalizeInput(unittest2.TestCase):
    def test_source_specified(self):
        E = self.assertEquals
        normalize_input = rdkit.normalize_input
        E(normalize_input(), ("smi", decompressors.Uncompressed))
        E(normalize_input("blah.smi"), ("smi", decompressors.Uncompressed))
        E(normalize_input("blah.smiles"), ("smi", decompressors.Uncompressed))
        E(normalize_input("blah.ism"), ("smi", decompressors.Uncompressed))
        E(normalize_input("blah.sdf"), ("sdf", decompressors.Uncompressed))
        E(normalize_input("blah.sdf2"), ("smi", decompressors.Uncompressed))

        E(normalize_input("blah.can.gz"), ("smi", decompressors.GzipDecompressor))
        E(normalize_input("blah.sdf.gz"), ("sdf", decompressors.GzipDecompressor))
        E(normalize_input("blah.sdf2.gz"), ("smi", decompressors.GzipDecompressor))
        E(normalize_input(source="blah.sdf2.gz"), ("smi", decompressors.GzipDecompressor))

        E(normalize_input("blah.smi.bz2"), ("smi", decompressors.Bzip2Decompressor))
        E(normalize_input("blah.sdf.bz2"), ("sdf", decompressors.Bzip2Decompressor))
        E(normalize_input("blah.sdf2.bz2"), ("smi", decompressors.Bzip2Decompressor))

        # Make sure it isn't just checking the last three characters
        E(normalize_input("blahsdf"), ("smi", decompressors.Uncompressed))

    def test_format_specified(self):
        E = self.assertEquals
        normalize_input = rdkit.normalize_input
        
        E(normalize_input(format="sdf"), ("sdf", decompressors.Uncompressed))

        # Check that the source doesn't affect anything
        E(normalize_input("blah.smi", format="sdf"), ("sdf", decompressors.Uncompressed))
        E(normalize_input("blah.sdf", format="sdf"), ("sdf", decompressors.Uncompressed))
        E(normalize_input("blah.sdf", "sdf"), ("sdf", decompressors.Uncompressed))

        # Various decompression
        E(normalize_input("blah.smi", format="sdf.gz"),
                                 ("sdf", decompressors.GzipDecompressor))
        E(normalize_input("blah.smi.bz2", format="sdf.gz"),
                                 ("sdf", decompressors.GzipDecompressor))


    def test_with_decompressor(self):
        E = self.assertEquals
        normalize_input = rdkit.normalize_input

        E(normalize_input("x.iso", "sdf", "bzip2"), ("sdf", decompressors.Bzip2Decompressor))
        E(normalize_input("x.iso", "smi", "bzip2"), ("smi", decompressors.Bzip2Decompressor))
        self.assertRaises(TypeError, normalize_input, "x.iso", "smi.gz", "bzip2")

        E(normalize_input("x.mdl", None, "bzip2"), ("sdf", decompressors.Bzip2Decompressor))
        E(normalize_input(None, None, "bzip2"), ("smi", decompressors.Bzip2Decompressor))
        E(normalize_input(None, None, "gzip"), ("smi", decompressors.GzipDecompressor))
        E(normalize_input(None, None, decompressors.GzipDecompressor),
                              ("smi", decompressors.GzipDecompressor))
        E(normalize_input(None, "sdf", decompressors.GzipDecompressor),
                              ("sdf", decompressors.GzipDecompressor))
        

    def test_unknown_format(self):
        self.assertRaises(TypeError, rdkit.normalize_input, "blah.smi", "mdl")
    
if __name__ == "__main__":
    unittest2.main()
