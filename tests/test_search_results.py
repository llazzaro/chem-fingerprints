from __future__ import with_statement

import unittest2
import random

from chemfp.search import SearchResults
from chemfp.fps_search import FPSSearchResults, FPSSearchResult

try:
    next
except NameError:
    # Compatibility with Python 2.5
    def next(it):
        return it.next()
  

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

class TestCase(unittest2.TestCase):
    def assertListEquals(self, left, right, msg=None):
        self.assertEquals(list(left), right, msg)

class CreateSearchResults(object):
    def _create(self, i, values):
        results = SearchResults(i)
        for value in values:
            results._add_hit(*value)
        return results

class CreateFPSSearchResults(object):
    def _create(self, i, values):
        all_ids = [[] for x in range(i)]
        all_scores = [[] for x in range(i)]
        results = FPSSearchResults(i)
        for (row, id, score) in values:
            all_ids[row].append(id)
            all_scores[row].append(score)
        return FPSSearchResults([FPSSearchResult(rows, scores)
                                      for (rows, scores) in zip(all_ids, all_scores)])

class TestBasicAPI(object):
    def test_len(self):
        for i in (0, 1, 4, 5):
            results = SearchResults(i)
            self.assertEquals(len(results), i)

    def test_row_len(self):
        results = self._create(5, [
            (0, 1, 0.1),
            (1, 2, 0.2),
            (1, 3, 0.25),
            (2, 1, 0.15),
            (2, 5, 0.7),
            (2, 6, 0.8),
            (3, 8, 0.9)])
        self.assertEquals(len(results[0]), 1)
        self.assertEquals(len(results[1]), 2)
        self.assertEquals(len(results[2]), 3)
        self.assertEquals(len(results[3]), 1)
        self.assertEquals(len(results[4]), 0)

        self.assertEquals(len(results[-5]), 1)
        self.assertEquals(len(results[-4]), 2)
        self.assertEquals(len(results[-3]), 3)
        self.assertEquals(len(results[-2]), 1)
        self.assertEquals(len(results[-1]), 0)
        
    def test_negative_index(self):
        results = self._create(3, [
            (0, 1, 0.0),
            (1, 12, 1.0)])
        self.assertListEquals(results[1], [(12, 1.0)])
        self.assertListEquals(results[-2], [(12, 1.0)])
        self.assertListEquals(results[-3], [(1, 0.0)])

    def test_clear(self):
        results = self._create(3, [
            (0, 1, 0.0),
            (1, 12, 1.0)])
        self.assertTrue(results[0])
        self.assertTrue(results[1])
        results.clear_all()
        self.assertFalse(results[0])
        self.assertListEquals(results[0], [])
        self.assertFalse(results[1])

    def test_clear_row(self):
        results = self._create(3, [
            (0, 1, 0.0),
            (1, 12, 1.0)])
        self.assertListEquals(results[1], [(12, 1.0)])
        results[1].clear()
        self.assertListEquals(results[0], [(1, 0.0)])
        self.assertTrue(results[0])
        self.assertEquals(len(results[0]), 1)
        self.assertFalse(results[1])
        results[0].clear()
        self.assertFalse(results[0])
        self.assertEquals(len(results[0]), 0)
        self.assertFalse(results[1])

    def test_clear_negative_row(self):
        results = self._create(3, [
            (0, 1, 0.0),
            (1, 12, 1.0)])
        self.assertListEquals(results[1], [(12, 1.0)])
        results[-2].clear()
        self.assertListEquals(results[0], [(1, 0.0)])
        self.assertListEquals(results[1], [])
        results[-3].clear()
        self.assertListEquals(results[0], [])
        self.assertListEquals(results[1], [])

    def test_unknown_ordering(self):
        results = self._create(1, [(0, 1, 0.9)])
        with self.assertRaisesRegexp(ValueError, "Unknown sort order"):
            results.reorder_all("blah")

class TestArenaBasicAPI(TestCase, CreateSearchResults, TestBasicAPI):
    pass

class TestFPSBasicAPI(TestCase, CreateFPSSearchResults, TestBasicAPI):
    pass



class TestIterAPI(TestCase):
    def setUp(self):
        results = SearchResults(2)
        results.target_ids = map(str, range(40))
        all_scores = [random_scores[:10],
                      random_scores[30:35]]
        # First row has columns 0, 1, 2, ..., 9
        # second row has columns 0, 2, 4, 6, 8
        for row, scores in enumerate(all_scores):
            for i, score in enumerate(scores):
                results._add_hit(row, i*(row+1), score)
        self.results = results
        
    def test_iter_indices(self):
        results = self.results
        it = results.iter_indices()
        self.assertListEquals(next(it), range(10))
        self.assertListEquals(next(it), [0, 2, 4, 6, 8])
        with self.assertRaisesRegexp(StopIteration, ""):
            next(it)

    def test_iter_ids(self):
        results = self.results
        it = results.iter_ids()
        self.assertListEquals(next(it), map(str, range(10)))
        self.assertListEquals(next(it), map(str, [0, 2, 4, 6, 8]))
        with self.assertRaisesRegexp(StopIteration, ""):
            next(it)

    def test_iter_scores(self):
        results = self.results
        it = results.iter_scores()
        self.assertListEquals(next(it), random_scores[:10])
        self.assertListEquals(next(it), random_scores[30:35])
        with self.assertRaisesRegexp(StopIteration, ""):
            next(it)

    def test_iter_indices_and_scores(self):
        results = self.results
        it = results.iter_indices_and_scores()
        self.assertListEquals(next(it), zip(range(10), random_scores[:10]))
        self.assertListEquals(next(it), zip([0, 2, 4, 6, 8], random_scores[30:35]))
        with self.assertRaisesRegexp(StopIteration, ""):
            next(it)

    def test_iter_ids_and_scores(self):
        results = self.results
        it = results.iter_ids_and_scores()
        self.assertListEquals(next(it), zip(map(str, range(10)), random_scores[:10]))
        self.assertListEquals(next(it), zip(map(str, [0, 2, 4, 6, 8]), random_scores[30:35]))
        with self.assertRaisesRegexp(StopIteration, ""):
            next(it)


class TestErrors(TestCase):
    def test_bad_order(self):
        results = SearchResults(5)
        with self.assertRaisesRegexp(ValueError, "Unknown sort order"):
            results.reorder_all("xyzzy")

    def test_bad_row_order(self):
        results = SearchResults(5)
        with self.assertRaisesRegexp(ValueError, "Unknown sort order"):
            results[0].reorder("xyzzy")

    def test_index_out_of_range(self):
        results = SearchResults(5)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results[5]
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results[98]

    def test_illegal_negative_index(self):
        results = SearchResults(3)
        with self.assertRaisesRegexp(IndexError, "row index is out of range"):
            results[-4]


class TestGetHitInfo(TestCase):
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
        self.assertListEquals(self.results[0], [])
        self.assertListEquals(self.results[1], [(1, 0.1)])
        self.assertListEquals(self.results[3], [(3, 0.2), (4, 0.5), (5, 0.6), (6, 0.7)])
        self.assertListEquals(self.results[-2], [(7, 1.0)])

    def test_get_indices(self):
        self.assertListEquals(self.results[0].get_indices(), [])
        self.assertListEquals(self.results[1].get_indices(), [1])
        self.assertListEquals(self.results[3].get_indices(), [3, 4, 5, 6])
        self.assertListEquals(self.results[-2].get_indices(), [7])

    def test_get_scores(self):
        self.assertListEquals(self.results[0].get_scores(), [])
        self.assertListEquals(self.results[1].get_scores(), [0.1])
        self.assertListEquals(self.results[3].get_scores(), [0.2, 0.5, 0.6, 0.7])
        self.assertListEquals(self.results[-2].get_scores(), [1.0])


    def test_get_indices_and_scores(self):
        self.assertListEquals(self.results[0].get_indices_and_scores(), [])
        self.assertListEquals(self.results[1].get_indices_and_scores(), [(1, 0.1)])
        self.assertListEquals(self.results[3].get_indices_and_scores(),
                              [(3, 0.2), (4, 0.5), (5, 0.6), (6, 0.7)])
        self.assertListEquals(self.results[-2].get_indices_and_scores(), [(7, 1.0)])

    def test_get_ids_and_scores(self):
        self.results.target_ids = map(str, range(8))
        self.assertListEquals(self.results[0].get_ids_and_scores(), [])
        self.assertListEquals(self.results[1].get_ids_and_scores(), [("1", 0.1)])
        self.assertListEquals(self.results[3].get_ids_and_scores(),
                              [("3", 0.2), ("4", 0.5), ("5", 0.6), ("6", 0.7)])
        self.assertListEquals(self.results[-2].get_ids_and_scores(), [("7", 1.0)])

    def test_get_ids_and_scores_missing_target_ids(self):
        with self.assertRaisesRegexp(TypeError, "target_ids are not available"):
            self.assertListEquals(self.results[-2].get_ids_and_scores(), [("7", 1.0)])            

class TestGetHitInfoRow(TestCase):
    def setUp(self):
        results = SearchResults(4)
        results._add_hit(1, 1, 0.1)
        results._add_hit(2, 7, 1.0)
        results._add_hit(3, 3, 0.2)
        results._add_hit(3, 4, 0.5)
        results._add_hit(3, 5, 0.6)
        results._add_hit(3, 6, 0.7)
        self.results = results

    def test_get_indices(self):
        self.assertListEquals(self.results[0].get_indices(), [])
        self.assertListEquals(self.results[1].get_indices(), [1])
        self.assertListEquals(self.results[3].get_indices(), [3, 4, 5, 6])
        self.assertListEquals(self.results[-2].get_indices(), [7])

    def test_get_scores(self):
        self.assertListEquals(self.results[0].get_scores(), [])
        self.assertListEquals(self.results[1].get_scores(), [0.1])
        self.assertListEquals(self.results[3].get_scores(), [0.2, 0.5, 0.6, 0.7])
        self.assertListEquals(self.results[-2].get_scores(), [1.0])


    def test_get_indices_and_scores(self):
        self.assertListEquals(self.results[0].get_indices_and_scores(), [])
        self.assertListEquals(self.results[1].get_indices_and_scores(), [(1, 0.1)])
        self.assertListEquals(self.results[3].get_indices_and_scores(),
                              [(3, 0.2), (4, 0.5), (5, 0.6), (6, 0.7)])
        self.assertListEquals(self.results[-2].get_indices_and_scores(), [(7, 1.0)])

    def test_get_ids_and_scores(self):
        self.results.target_ids = map(str, range(8))
        self.assertListEquals(self.results[0].get_ids_and_scores(), [])
        self.assertListEquals(self.results[1].get_ids_and_scores(), [("1", 0.1)])
        self.assertListEquals(self.results[3].get_ids_and_scores(),
                              [("3", 0.2), ("4", 0.5), ("5", 0.6), ("6", 0.7)])
        self.assertListEquals(self.results[-2].get_ids_and_scores(), [("7", 1.0)])

    def test_get_ids_and_scores_missing_target_ids(self):
        with self.assertRaisesRegexp(TypeError, "target_ids are not available"):
            self.assertListEquals(self.results[-2].get_ids_and_scores(), [("7", 1.0)])            


_get_sort_key = {
    "increasing-score": lambda (index, score): (score, index),
    "decreasing-score": lambda (index, score): (-score, index),
    "increasing-index": lambda (index, score): index,
    "decreasing-index": lambda (index, score): -index,
}
class TestSortOrder(object):
    def test_size_0(self):
        results = self._create(5, [])
        results.reorder_all()
        self.assertListEquals(results[0], [])
    def test_size_1(self):
        results = self._create(5, [
            (1, 5, 0.2),])
        results.reorder_all()
        self.assertListEquals(results[1], [(5, 0.2)])
    def test_size_2(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),])
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4)])
        results.reorder_all("increasing-score")
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4)])
        results.reorder_all("decreasing-score")
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2)])

    def test_default_ordering_2(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),])
        results.reorder_all()
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2)])

    def test_size_3(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),
            (1, 7, 0.2),])
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4), (7, 0.2)])
        results.reorder_all("increasing-score")
        self.assertListEquals(results[1], [(5, 0.2), (7, 0.2), (6, 0.4)])
        results.reorder_all("decreasing-score")
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2), (7, 0.2)])

class TestArenaTestSortOrder(TestCase, CreateSearchResults, TestSortOrder):
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
                results.reorder_all(name)
                expected.sort(key = _get_sort_key[name])
                self.assertListEquals(results[0], expected, "error in %s (300)" % (name,))

    def test_random_values(self):
        # The underlying timsort does merge sorts of 64 element
        # blocks.  Hence some of the code is not exercised unless the
        # input is at least 128 elements long.
        for size in (3, 5, 10, 20, 70, 100, 400):
            results = SearchResults(1)
            expected = []
            for i in range(size):
                score = random_scores[i]
                expected.append((i, score))
                results._add_hit(0, i, score)

            self.assertListEquals(results[0], expected)
            for name in ("increasing-score", "decreasing-score",
                         "increasing-index", "decreasing-index"):
                results.reorder_all(name)
                expected.sort(key = _get_sort_key[name])
                self.assertListEquals(results[0], expected, "error in %s:%d" % (name, size))
    
class TestFPSSortOrder(TestCase, CreateFPSSearchResults, TestSortOrder):
    def test_order_by_id(self):
        results = self._create(2, [
            (0, "one", 0.4),
            (0, "two", 0.5),
            (0, "three", 0.6),

            (1, "ett", 0.3),
            (1, "tvaa", 0.2),
            (1, "tre", 0.1),])
        results.reorder_all("increasing-id")
        self.assertListEquals(results[0], [("one", 0.4), ("three", 0.6), ("two", 0.5)])
        self.assertListEquals(results[1], [("ett", 0.3), ("tre", 0.1), ("tvaa", 0.2)])
        
        results.reorder_all("decreasing-id")
        self.assertListEquals(results[0], [("one", 0.4), ("three", 0.6), ("two", 0.5)][::-1])
        self.assertListEquals(results[1], [("ett", 0.3), ("tre", 0.1), ("tvaa", 0.2)][::-1])


class TestSortOrderRow(object):
    def test_size_2(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),])
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4)])
        for result in results:
            result.reorder("increasing-score")
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4)])
        for result in results:
            result.reorder("decreasing-score")
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2)])

    def test_default_ordering_2(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),])
        for result in results:
            result.reorder()
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2)])

    def test_size_3(self):
        results = self._create(5, [
            (1, 5, 0.2),
            (1, 6, 0.4),
            (1, 7, 0.2),])
        self.assertListEquals(results[1], [(5, 0.2), (6, 0.4), (7, 0.2)])
        for result in results:
            result.reorder("increasing-score")
        self.assertListEquals(results[1], [(5, 0.2), (7, 0.2), (6, 0.4)])
        for result in results:
            result.reorder("decreasing-score")
        self.assertListEquals(results[1], [(6, 0.4), (5, 0.2), (7, 0.2)])

    def test_reorder_row(self):
        results = self._create(2, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),
        
            (1, 6, 0.1),
            (1, 7, 0.8),
            (1, 8, 0.6),])

        self.assertListEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertListEquals(results[1], [(6, 0.1), (7, 0.8), (8, 0.6)])

        results[1].reorder("increasing-score")
        self.assertListEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertListEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        results[0].reorder("decreasing-score")
        self.assertListEquals(results[0], [(2, 0.8), (3, 0.6), (1, 0.1)])
        self.assertListEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        # Check that the default works
        results[0].reorder("increasing-score")  # ensure the default only affects one row
        results[1].reorder()
        self.assertListEquals(results[0], [(1, 0.1), (3, 0.6), (2, 0.8)])
        self.assertListEquals(results[1], [(7, 0.8), (8, 0.6), (6, 0.1)])

    def test_sort_negative_row(self):
        results = self._create(2, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),
        
            (1, 6, 0.1),
            (1, 7, 0.8),
            (1, 8, 0.6),])

        self.assertListEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertListEquals(results[1], [(6, 0.1), (7, 0.8), (8, 0.6)])

        results[-1].reorder("increasing-score")
        self.assertListEquals(results[0], [(1, 0.1), (2, 0.8), (3, 0.6)])
        self.assertListEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        results[-2].reorder("decreasing-score")
        self.assertListEquals(results[0], [(2, 0.8), (3, 0.6), (1, 0.1)])
        self.assertListEquals(results[1], [(6, 0.1), (8, 0.6), (7, 0.8)])

        results[1].reorder()  # default is decreasing score
        self.assertListEquals(results[0], [(2, 0.8), (3, 0.6), (1, 0.1)])
        self.assertListEquals(results[1], [(7, 0.8), (8, 0.6), (6, 0.1)])
        
class TestArenaTestSortOrderRow(TestCase, CreateSearchResults, TestSortOrderRow):
    pass
    
class TestFPSSortOrderRow(TestCase, CreateFPSSearchResults, TestSortOrderRow):
    pass


class TestMoveClosestFirst(object):
    def test_empty(self):
        results = self._create(2, [])
        results.reorder_all("move-closest-first")
        self.assertEquals(len(results), 2)
        self.assertEquals(len(results[0]), 0)
        self.assertEquals(len(results[1]), 0)

    def test_one(self):
        results = self._create(2, [
            (0, 9, 0.1),
            (1, 8, 0.8)])
        results.reorder_all("move-closest-first")

        self.assertListEquals(results[0], [(9, 0.1)])
        self.assertListEquals(results[1], [(8, 0.8)])

    def test_two(self):
        results = self._create(2, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (1, 2, 0.8),
            (1, 3, 0.6),])

        results.reorder_all("move-closest-first")

        self.assertListEquals(results[0], [(2, 0.8), (1, 0.1)])
        self.assertListEquals(results[1], [(2, 0.8), (3, 0.6)])

    def test_three(self):
        results = self._create(3, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),

            (1, 12, 0.8),
            (1, 22, 0.1),
            (1, 32, 0.6),

            (2, 12, 0.6),
            (2, 22, 0.1),
            (2, 32, 0.8),])
        
        results.reorder_all("move-closest-first")
        
        self.assertListEquals(results[0], [(2, 0.8), (1, 0.1), (3, 0.6)])
        self.assertListEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertListEquals(results[2], [(32, 0.8), (22, 0.1), (12, 0.6)])

    def test_row(self):
        results = self._create(3, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),

            (1, 12, 0.8),
            (1, 22, 0.1),
            (1, 32, 0.6),
        
            (2, 12, 0.6),
            (2, 22, 0.1),
            (2, 32, 0.8),])
        
        results[0].reorder("move-closest-first")
        self.assertListEquals(results[0], [(2, 0.8), (1, 0.1), (3, 0.6)])
        self.assertListEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertListEquals(results[2], [(12, 0.6), (22, 0.1), (32, 0.8)])

        results[-1].reorder("move-closest-first")
        self.assertListEquals(results[0], [(2, 0.8), (1, 0.1), (3, 0.6)])
        self.assertListEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertListEquals(results[2], [(32, 0.8), (22, 0.1), (12, 0.6)])

class TestArenaMoveClosestFirst(TestCase, CreateSearchResults, TestMoveClosestFirst):
    pass
    
class TestFPSMoveClosestFirst(TestCase, CreateFPSSearchResults, TestMoveClosestFirst):
    pass
        

class TestReverse(object):
    def test_empty(self):
        results = self._create(2, [])
        results.reorder_all("reverse")
        self.assertEquals(len(results), 2)
        self.assertEquals(len(results[0]), 0)
        self.assertEquals(len(results[1]), 0)
        
    def test_one(self):
        results = self._create(2, [
            (0, 9, 0.1),
            (1, 8, 0.8)])
        results.reorder_all("reverse")
        self.assertListEquals(results[0], [(9, 0.1)])
        self.assertListEquals(results[1], [(8, 0.8)])

    def test_two(self):
        results = self._create(2, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (1, 2, 0.8),
            (1, 3, 0.6)])
        results.reorder_all("reverse")
        self.assertListEquals(results[0], [(2, 0.8), (1, 0.1)])
        self.assertListEquals(results[1], [(3, 0.6), (2, 0.8)])

    def test_three(self):
        results = self._create(3, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),

            (1, 12, 0.8),
            (1, 22, 0.1),
            (1, 32, 0.6),

            (2, 12, 0.6),
            (2, 32, 0.1),
            (2, 22, 0.8),])
        
        results.reorder_all("reverse")
        
        self.assertListEquals(results[0], [(3, 0.6), (2, 0.8), (1, 0.1)])
        self.assertListEquals(results[1], [(32, 0.6), (22, 0.1), (12, 0.8)])
        self.assertListEquals(results[2], [(22, 0.8), (32, 0.1), (12, 0.6)])

    def test_row(self):
        results = self._create(3, [
            (0, 1, 0.1),
            (0, 2, 0.8),
            (0, 3, 0.6),

            (1, 12, 0.8),
            (1, 22, 0.1),
            (1, 32, 0.6),
        
            (2, 12, 0.6),
            (2, 32, 0.1),
            (2, 22, 0.8)])
        
        results[0].reorder("reverse")
        self.assertListEquals(results[0], [(3, 0.6), (2, 0.8), (1, 0.1)])
        self.assertListEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertListEquals(results[2], [(12, 0.6), (32, 0.1), (22, 0.8)])

        results[-1].reorder("reverse")
        self.assertListEquals(results[0], [(3, 0.6), (2, 0.8), (1, 0.1)])
        self.assertListEquals(results[1], [(12, 0.8), (22, 0.1), (32, 0.6)])
        self.assertListEquals(results[2], [(22, 0.8), (32, 0.1), (12, 0.6)])

class TestArenaReverse(TestCase, CreateSearchResults, TestReverse):
    pass
    
class TestFPSReverse(TestCase, CreateFPSSearchResults, TestReverse):
    pass

if __name__ == "__main__":
    unittest2.main()
