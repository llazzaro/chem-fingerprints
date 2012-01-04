from __future__ import with_statement

import unittest2
import random

from _chemfp import SearchResults


random_scores = [
  0.676, 0.384, 0.740, 0.970, 0.148, 0.361, 0.621, 0.715, 0.698,
  0.009, 0.667, 0.760, 0.743, 0.807, 0.772, 0.074, 0.622, 0.218,
  0.594, 0.247, 0.680, 0.214, 0.721, 0.590, 0.433, 0.725, 0.917,
  0.401, 0.818, 0.381, 0.039, 0.214, 0.133, 0.014, 0.072, 0.254,
  0.515, 0.965, 0.145, 0.548, 0.468, 0.205, 0.631, 0.132, 0.710,
  0.367, 0.313, 0.866, 0.611, 0.640, 0.727, 0.910, 0.057, 0.619,
  0.160, 0.390, 0.868, 0.101, 0.525, 0.689, 0.945, 0.473, 0.448,
  0.705, 0.399, 0.731, 0.214, 0.575, 0.721, 0.867, 0.514, 0.801,
  0.415, 0.742, 0.628, 0.686, 0.117, 0.016, 0.411, 0.336, 0.447,
  0.774, 0.028, 0.283, 0.937, 0.341, 0.348, 0.404, 0.956, 0.391,
  0.822, 0.976, 0.162, 0.422, 0.260, 0.688, 0.596, 0.298, 0.927,
  0.412, 0.979, 0.180, 0.258, 0.779, 0.893, 0.367, 0.219, 0.658,
  0.084, 0.966, 0.264, 0.024, 0.795, 0.703, 0.092, 0.007, 0.463,
  0.028, 0.567, 0.815, 0.403, 0.084, 0.760, 0.738, 0.125, 0.067,
  0.200, 0.044, 0.307, 0.696, 0.314, 0.244, 0.420, 0.135, 0.741,
  0.770, 0.047, 0.678, 0.186, 0.704, 0.732, 0.796, 0.017, 0.316,
  0.377, 0.256, 0.866, 0.964, 0.651, 0.046, 0.073, 0.787, 0.043,
  0.115, 0.269, 0.171, 0.374, 0.572, 0.448, 0.733, 0.845, 0.366,
  0.777, 0.057, 0.009, 0.145, 0.755, 0.939, 0.417, 0.813, 0.658,
  0.985, 0.866, 0.171, 0.425, 0.849, 0.316, 0.981, 0.689, 0.665,
  0.872, 0.342, 0.326, 0.057, 0.755, 0.903, 0.111, 0.618, 0.980,
  0.313, 0.829, 0.566, 0.876, 0.456, 0.954, 0.122, 0.598, 0.616,
  0.735, 0.319, 0.482, 0.680, 0.437, 0.025, 0.281, 0.688, 0.859,
  0.472, 0.321, 0.048, 0.601, 0.654, 0.991, 0.918, 0.754, 0.832,
  0.352, 0.036, 0.184, 0.089, 0.534, 0.875, 0.651, 0.482, 0.135,
  0.958, 0.805, 0.999, 0.655, 0.373, 0.092, 0.906, 0.919, 0.464,
  0.588, 0.752, 0.268, 0.907, 0.936, 0.215, 0.551, 0.848, 0.324,
  0.533, 0.787, 0.369, 0.695, 0.550, 0.594, 0.775, 0.731, 0.774,
  0.881, 0.445, 0.034, 0.946, 0.979, 0.051, 0.494, 0.247, 0.475,
  0.650, 0.452, 0.482, 0.691, 0.018, 0.191, 0.773, 0.859, 0.240,
  0.782, 0.321, 0.613, 0.507, 0.453, 0.438, 0.461, 0.707, 0.973,
  0.621, 0.882, 0.216, 0.232, 0.832, 0.456, 0.383, 0.507, 0.100,
  0.422, 0.807, 0.008, 0.525, 0.962, 0.912, 0.329, 0.929, 0.041,
  0.107, 0.897, 0.810, 0.666, 0.560, 0.593, 0.934, 0.223, 0.460,
  0.447, 0.555, 0.649, 0.210, 0.297, 0.056, 0.854, 0.007, 0.024,
  0.914, 0.547, 0.620, 0.981, 0.369, 0.877, 0.320, 0.564, 0.767,
  0.672, 0.731, 0.386, 0.136, 0.569, 0.870, 0.950, 0.905, 0.504,
  0.130, 0.477, 0.850, 0.823, 0.765, 0.378, 0.065, 0.694, 0.733,
  0.283, 0.391, 0.084, 0.123, 0.362, 0.720, 0.327, 0.969, 0.077,
  0.179, 0.732, 0.747, 0.621, 0.784, 0.116, 0.570, 0.964, 0.064,
  0.138, 0.696, 0.596, 0.308, 0.256, 0.491, 0.678, 0.775, 0.958,
  0.202, 0.342, 0.426, 0.946, 0.969, 0.236, 0.435, 0.231, 0.345,
  0.242, 0.096, 0.146, 0.142, 0.723, 0.035, 0.449, 0.208, 0.123,
  0.596, 0.499, 0.797, 0.146, 0.115, 0.547, 0.678, 0.989, 0.200,
  0.042, 0.971, 0.732, 0.775]


class TestBasicAPI(unittest2.TestCase):
    def test_len(self):
        for i in (0, 1, 4, 5):
            results = SearchResults(i)
            self.assertEquals(len(results), i)

    def test_size(self):
        results = SearchResults(5)
        results._add_hit(0, 1, 0.1)
        results._add_hit(1, 2, 0.2)
        results._add_hit(1, 3, 0.25)
        results._add_hit(2, 1, 0.15)
        results._add_hit(2, 5, 0.7)
        results._add_hit(2, 6, 0.8)
        results._add_hit(3, 8, 0.9)
        self.assertEquals(results.size(0), 1)
        self.assertEquals(results.size(1), 2)
        self.assertEquals(results.size(2), 3)
        self.assertEquals(results.size(3), 1)
        self.assertEquals(results.size(4), 0)

        self.assertEquals(results.size(-5), 1)
        self.assertEquals(results.size(-4), 2)
        self.assertEquals(results.size(-3), 3)
        self.assertEquals(results.size(-2), 1)
        self.assertEquals(results.size(-1), 0)
        
    def test_negative_index(self):
        results = SearchResults(3)
        results._add_hit(0, 1, 0.0)
        results._add_hit(1, 12, 1.0)
        self.assertEquals(results[1], [(12, 1.0)])
        self.assertEquals(results[-2], [(12, 1.0)])
        self.assertEquals(results[-3], [(1, 0.0)])

    def test_clear(self):
        results = SearchResults(3)
        results._add_hit(0, 1, 0.0)
        results._add_hit(1, 12, 1.0)
        self.assertNotEquals(results[0], [])
        self.assertNotEquals(results[1], [])
        results.clear()
        self.assertEquals(results[0], [])
        self.assertEquals(results[1], [])

    def test_clear_row(self):
        results = SearchResults(3)
        results._add_hit(0, 1, 0.0)
        results._add_hit(1, 12, 1.0)
        self.assertEquals(results[1], [(12, 1.0)])
        results.clear_row(1)
        self.assertEquals(results[0], [(1, 0.0)])
        self.assertEquals(results[1], [])
        results.clear_row(0)
        self.assertEquals(results[0], [])
        self.assertEquals(results[1], [])

    def test_clear_negative_row(self):
        results = SearchResults(3)
        results._add_hit(0, 1, 0.0)
        results._add_hit(1, 12, 1.0)
        self.assertEquals(results[1], [(12, 1.0)])
        results.clear_row(-2)
        self.assertEquals(results[0], [(1, 0.0)])
        self.assertEquals(results[1], [])
        results.clear_row(-3)
        self.assertEquals(results[0], [])
        self.assertEquals(results[1], [])




class TestErrors(unittest2.TestCase):
    def test_bad_sort_order(self):
        results = SearchResults(5)
        with self.assertRaisesRegexp(ValueError, "Unknown sort order"):
            results.sort("increasing")

    def test_index_out_of_range(self):
        results = SearchResults(5)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results[5]

    def test_illegal_negative_index(self):
        results = SearchResults(3)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results[-4]

    def test_illegal_clear(self):
        results = SearchResults(3)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.clear_row(3)

        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.clear_row(-4)

    def test_illegal_indices(self):
        results = SearchResults(4)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.get_indices(4)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.get_indices(-5)
        
    def test_illegal_scores(self):
        results = SearchResults(4)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.get_scores(4)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results.get_scores(-5)


class TestGetHitInfo(unittest2.TestCase):
    def setUp(self):
        results = SearchResults(4)
        results._add_hit(1, 1, 0.1)
        results._add_hit(2, 7, 1.0)
        results._add_hit(3, 3, 0.2)
        results._add_hit(3, 4, 0.5)
        results._add_hit(3, 5, 0.6)
        results._add_hit(3, 6, 0.7)
        self.results = results

    def test_index(self):
        self.assertEquals(self.results[0], [])
        self.assertEquals(self.results[1], [(1, 0.1)])
        self.assertEquals(self.results[3], [(3, 0.2), (4, 0.5), (5, 0.6), (6, 0.7)])
        self.assertEquals(self.results[-2], [(7, 1.0)])

    def test_get_indices(self):
        self.assertEquals(self.results.get_indices(0), [])
        self.assertEquals(self.results.get_indices(1), [1])
        self.assertEquals(self.results.get_indices(3), [3, 4, 5, 6])
        self.assertEquals(self.results.get_indices(-2), [7])

    def test_get_scores(self):
        self.assertEquals(self.results.get_scores(0), [])
        self.assertEquals(self.results.get_scores(1), [0.1])
        self.assertEquals(self.results.get_scores(3), [0.2, 0.5, 0.6, 0.7])
        self.assertEquals(self.results.get_scores(-2), [1.0])
        

_get_sort_key = {
    "increasing-score": lambda (index, score): (score, index),
    "decreasing-score": lambda (index, score): (-score, index),
    "increasing-index": lambda (index, score): index,
    "decreasing-index": lambda (index, score): -index,
}
class TestSortOrder(unittest2.TestCase):
    def test_size_0(self):
        results = SearchResults(5)
        self.assertEquals(results[0], [])
    def test_size_1(self):
        results = SearchResults(5)
        results._add_hit(1, 5, 0.2)
        self.assertEquals(results[1], [(5, 0.2)])
    def test_size_2(self):
        results = SearchResults(5)
        results._add_hit(1, 5, 0.2)
        results._add_hit(1, 6, 0.4)
        self.assertEquals(results[1], [(5, 0.2), (6, 0.4)])
        results.sort("increasing-score")
        self.assertEquals(results[1], [(5, 0.2), (6, 0.4)])
        results.sort("decreasing-score")
        self.assertEquals(results[1], [(6, 0.4), (5, 0.2)])

    def test_size_3(self):
        results = SearchResults(5)
        results._add_hit(1, 5, 0.2)
        results._add_hit(1, 6, 0.4)
        results._add_hit(1, 7, 0.2)
        self.assertEquals(results[1], [(5, 0.2), (6, 0.4), (7, 0.2)])
        results.sort("increasing-score")
        self.assertEquals(results[1], [(5, 0.2), (7, 0.2), (6, 0.4)])
        results.sort("decreasing-score")
        self.assertEquals(results[1], [(6, 0.4), (5, 0.2), (7, 0.2)])

    def test_random_values(self):
        # The underlying timsort does merge sorts of 64 element
        # blocks.  Hence some of the code is not exercised unless the
        # input is at 128 elements long.
        for size in (3, 5, 10, 20, 70, 100, 400):
            results = SearchResults(1)
            expected = []
            for i in range(size):
                score = random_scores[i]
                expected.append((i, score))
                results._add_hit(0, i, score)
            self.assertEquals(results[0], expected)
            for name in ("increasing-score", "decreasing-score",
                         "increasing-index", "decreasing-index"):
                results.sort(name)
                expected.sort(key = _get_sort_key[name])
                self.assertEquals(results[0], expected, "error in %s:%d" % (name, size))

    def test_index_as_secondary_sort(self):
        # Timsort preserves input order. test_random_values uses
        # sequentially ordered indicies so can't tell the difference
        # between input order and index order. Here I reverse the
        # order so I can really test tie-breaking.
        for name in ("increasing-score", "decreasing-score",
                     "increasing-index", "decreasing-index"):
            results = SearchResults(1)
            expected = []
            for i in range(300):
                score = random_scores[i]
                results._add_hit(0, 400-i, score)
                expected.append((400-i, score))
                results.sort(name)
                expected.sort(key = _get_sort_key[name])
                self.assertEquals(results[0], expected, "error in %s (300)" % (name,))

    def test_sort_row(self):
        results = SearchResults(2)
        results._add_hit(0, 1, 0.1)
        results._add_hit(0, 2, 0.8)
        results._add_hit(0, 3, 0.6)
        
        results._add_hit(1, 6, 0.1)
        results._add_hit(1, 7, 0.8)
        results._add_hit(1, 8, 0.6)

        self.assertEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertEquals(results[1], [(6, 0.1), (7, 0.8), (8, 0.6)])

        results.sort_row(1, "increasing-score")
        self.assertEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        results.sort_row(0, "decreasing-score")
        self.assertEquals(results[0], [(2, 0.8), (3, 0.6), (1, 0.1)])
        self.assertEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

    def test_sort_negative_row(self):
        results = SearchResults(2)
        results._add_hit(0, 1, 0.1)
        results._add_hit(0, 2, 0.8)
        results._add_hit(0, 3, 0.6)
        
        results._add_hit(1, 6, 0.1)
        results._add_hit(1, 7, 0.8)
        results._add_hit(1, 8, 0.6)

        self.assertEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertEquals(results[1], [(6, 0.1), (7, 0.8), (8, 0.6)])

        results.sort_row(-1, "increasing-score")
        self.assertEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        results.sort_row(-2, "decreasing-score")
        self.assertEquals(results[0], [(2, 0.8), (3, 0.6), (1, 0.1)])
        self.assertEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])
        

class TestMoveClosestFirst(unittest2.TestCase):
    def test_empty(self):
        results = SearchResults(2)
        results.sort("move-closest-first")
        self.assertEquals(len(results), 2)
        self.assertEquals(len(results[0]), 0)
        self.assertEquals(len(results[1]), 0)

    def test_one(self):
        results = SearchResults(2)
        results._add_hit(0, 9, 0.1)
        results._add_hit(1, 8, 0.8)
        results.sort("move-closest-first")

        results.sort("move-closest-first")

        self.assertEquals(results[0], [(9, 0.1)])
        self.assertEquals(results[1], [(8, 0.8)])

    def test_two(self):
        results = SearchResults(2)
        results._add_hit(0, 1, 0.1)
        results._add_hit(0, 2, 0.8)
        results._add_hit(1, 2, 0.8)
        results._add_hit(1, 3, 0.6)

        results.sort("move-closest-first")

        self.assertEquals(results[0], [(2, 0.8), (1, 0.1)])
        self.assertEquals(results[1], [(2, 0.8), (3, 0.6)])

    def test_three(self):
        results = SearchResults(3)
        results._add_hit(0, 1, 0.1)
        results._add_hit(0, 2, 0.8)
        results._add_hit(0, 3, 0.6)

        results._add_hit(1, 12, 0.8)
        results._add_hit(1, 22, 0.1)
        results._add_hit(1, 32, 0.6)

        results._add_hit(2, 12, 0.6)
        results._add_hit(2, 22, 0.1)
        results._add_hit(2, 32, 0.8)
        
        results.sort("move-closest-first")
        
        self.assertEquals(results[0], [(2, 0.8), (1, 0.1), (3, 0.6)])
        self.assertEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertEquals(results[2], [(32, 0.8), (22, 0.1), (12, 0.6)])
        

if __name__ == "__main__":
    unittest2.main()
