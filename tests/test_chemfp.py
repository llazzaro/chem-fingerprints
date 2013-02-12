import unittest2
import re
import os

import chemfp
import _chemfp

version_pattern = re.compile(r"\d+\.\d+(\.\d)?((a|b|p)\d+)?$")

class SystemTests(unittest2.TestCase):
    def test_version(self):
        m = version_pattern.match(_chemfp.version())
        self.assertNotEqual(m, None, "bad version: %s" % (_chemfp.version(),))


skip_omp = (chemfp.get_max_threads() == 1)

class OpenMPTests(unittest2.TestCase):
    def setUp(self):
        self._num_threads = chemfp.get_num_threads()
    def tearDown(self):
        chemfp.set_num_threads(self._num_threads)

    def test_num_threads_is_max_threads(self):
        self.assertEquals(chemfp.get_num_threads(), chemfp.get_max_threads())

    def test_set_to_zero(self):
        chemfp.set_num_threads(0)
        self.assertEquals(chemfp.get_num_threads(), 1)

    def test_set_to_one(self):
        chemfp.set_num_threads(1)
        self.assertEquals(chemfp.get_num_threads(), 1)

    def test_set_to_two(self):
        chemfp.set_num_threads(2)
        self.assertEquals(chemfp.get_num_threads(), 2)
    test_set_to_two = unittest2.skipIf(skip_omp, "Multiple OMP threads not available")(
        test_set_to_two)

    def test_set_to_max(self):
        chemfp.set_num_threads(chemfp.get_max_threads())
        self.assertEquals(chemfp.get_num_threads(), chemfp.get_max_threads())

    def test_set_beyond_max(self):
        chemfp.set_num_threads(chemfp.get_max_threads()+1)
        self.assertEquals(chemfp.get_num_threads(), chemfp.get_max_threads())


if __name__ == "__main__":
    unittest2.main()
