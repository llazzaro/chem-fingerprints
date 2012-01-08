# Test the symmetric code

import unittest2
from cStringIO import StringIO
import array

import chemfp
from chemfp import search, bitops

from support import fullpath, PUBCHEM_SDF, PUBCHEM_SDF_GZ

fps = chemfp.load_fingerprints(fullpath("queries.fps"))

zeros = chemfp.load_fingerprints(StringIO("""\
0000\tA
0000\tB
0001\tC
0002\tD
FFFE\tE
FFFF\tF
"""))

def slow_counts(counts, fps, threshold,
                query_start, query_end,
                target_start, target_end):
    N = len(fps)
    query_end = min(N, query_end)
    target_end = min(N, target_end)

    for row in range(query_start, query_end):
        row_fp = fps[row][1]
        for col in range(max(row+1, target_start), target_end):
            col_fp = fps[col][1]
            if bitops.byte_tanimoto(row_fp, col_fp) >= threshold:
                if row != col:
                    #print row, col, "X"
                    counts[row] += 1
                    counts[col] += 1

def slow_threshold_search(results, fps, threshold,
                          query_start, query_end,
                          target_start, target_end):
    N = len(fps)
    query_end = min(N, query_end)
    target_end = min(N, target_end)
    
    for row in range(query_start, query_end):
        row_fp = fps[row][1]
        for col in range(max(row+1, target_start), target_end):
            col_fp = fps[col][1]
            score = bitops.byte_tanimoto(row_fp, col_fp)
            if score >= threshold:
                if row != col:
                    #print row, col, "X"
                    results[row].append((col, score))
                    #results[col].append((row, score)) # only do the upper triangle


class TestCounts(unittest2.TestCase):
    #maxDiff = 0
    def test_symmetric(self):
        # query[i] always matches target[i] so x[i] will be at least one
        x = search.count_tanimoto_hits(fps, fps, 0.6)

        # This only processes the upper-triangle, and not the diagonal
        y = search.count_tanimoto_hits_symmetric(fps, 0.6)

        self.assertEquals(len(x), len(y))
        for i in range(len(x)):
            self.assertEquals(x[i]-1, y[i])


    def test_zeros(self):
        y = search.count_tanimoto_hits_symmetric(zeros, 0.9)
        self.assertEquals(list(y), [0, 0, 0, 0, 1, 1])
        y = search.count_tanimoto_hits_symmetric(zeros, 0.001)
        self.assertEquals(list(y), [0, 0, 1, 2, 2, 3])


    def test_partial_counts_with_zero_threshold(self):
        threshold = 0.0
        N = len(fps)
        counts = array.array("i", [0]*N)
        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold, 1, 5, 0, N)

        expected = [0] * N
        slow_counts(expected, fps, threshold, 1, 5, 0, N)
        self.assertSequenceEqual(counts, expected)

        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold,
                                                     3, 33, 0, N)
        slow_counts(expected, fps, threshold, 3, 33, 0, N)
        self.assertSequenceEqual(counts, expected)

    def test_partial_counts_in_columns(self):
        threshold = 0.6
        N = len(fps)
        counts = array.array("i", [0]*N)
        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold, 1, 5, 0, N)

        expected = [0] * N
        slow_counts(expected, fps, threshold, 1, 5, 0, N)
        self.assertSequenceEqual(counts, expected)

        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold,
                                                     3, 33, 0, N)
        slow_counts(expected, fps, threshold, 3, 33, 0, N)
        self.assertSequenceEqual(counts, expected)
        

    def test_partial_counts_in_rows(self):
        threshold = 0.2
        N = len(fps)
        counts = array.array("i", [0]*N)
        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold, 1, 5, 9, 30)

        expected = [0] * N
        slow_counts(expected, fps, threshold, 1, 5, 9, 30)
        self.assertSequenceEqual(counts, expected)

        search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold, 20, 26, 30, 50)
        slow_counts(expected, fps, threshold, 20, 26, 30, 50)
        self.assertSequenceEqual(counts, expected)

    def test_partial_counts(self):
        threshold = 0.2
        N = len(fps)
        counts = array.array("i", [0]*N)
        expected = [0] * N

        for i in range(0, N, 10):
            for j in range(i, N, 8):
                slow_counts(expected, fps, threshold, i, i+10, j, j+8)
                search.partial_count_tanimoto_hits_symmetric(counts, fps, threshold,
                                                             i, i+10, j, j+8)
                self.assertSequenceEqual(counts, expected)

        normal = search.count_tanimoto_hits_symmetric(fps, threshold)
        self.assertSequenceEqual(normal, expected)

def _compare_search_results(self, result, expected):
    q=list(result.iter_indices_and_scores())
    for i in range(len(result)):
        self.assertEquals(len(q[i]), len(expected[i]), "length of %d" % i)
        self.assertEquals(sorted(q[i]), sorted(expected[i]), "contents of %d" % i)


class TestThreshold(unittest2.TestCase):
    def test_upper_only(self):
        # query[i] always matches target[i] so x[i] will always contain i
        x = search.threshold_tanimoto_search(fps, fps, 0.9)
        x = list(x.iter_indices_and_scores())

        # This only processes the upper-triangle, and not the diagonal
        y = search.threshold_tanimoto_search_symmetric(fps, 0.9, include_lower_triangle=False)


        rows = list(row.get_indices_and_scores() for row in y)
        row_sizes = map(len, rows)
        # Move elements to the lower triangle
        for rowno, (row, row_size) in enumerate(zip(rows, row_sizes)):
            for (colno, score) in row[:row_size]:
                assert colno > rowno, (rowno, colno)
                rows[colno].append( (rowno, score) )

            # Fill in the diagonal
            row.append((rowno, 1.0))

            # Put into a consistent order
            row.sort()
            
            # Match with the NxM algorithm
            expected_row = x[rowno]
            expected_row.reorder("increasing-index")

            self.assertEquals(row, list(expected_row), rowno)

    def test_upper_and_lower(self):
        # query[i] always matches target[i] so x[i] will always contain i
        x = search.threshold_tanimoto_search(fps, fps, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.threshold_tanimoto_search_symmetric(fps, 0.9)

        for i, (x_row, y_row) in enumerate(zip(x, y)):
            x_row = x_row.get_indices_and_scores()
            y_row = y_row.get_indices_and_scores()
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()

            self.assertEquals(x_row, y_row)

    def test_partial_search_with_zero_threshold(self):
        threshold = 0.0
        N = len(fps)
        result = search.SearchResults(N, fps.arena_ids)
        search.partial_threshold_tanimoto_search_symmetric(result, fps, threshold,
                                                           1, 5, 0, N)

        expected = [[] for i in range(N)]
        slow_threshold_search(expected, fps, threshold, 1, 5, 0, N)
        _compare_search_results(self, result, expected)

    def test_partial_search_in_rows(self):
        threshold = 0.2
        N = len(fps)
        result = search.SearchResults(N, fps.arena_ids)
        search.partial_threshold_tanimoto_search_symmetric(result, fps, threshold,
                                                           1, 9, 0, N)

        expected = [[] for i in range(N)]
        slow_threshold_search(expected, fps, threshold, 1, 9, 0, N)
        _compare_search_results(self, result, expected)


    def test_partial_search_in_cols(self):
        threshold = 0.1
        N = len(fps)
        result = search.SearchResults(N, fps.arena_ids)
        search.partial_threshold_tanimoto_search_symmetric(result, fps, threshold,
                                                           0, N, 7, 69)

        expected = [[] for i in range(N)]
        slow_threshold_search(expected, fps, threshold, 0, N, 7, 69)
        _compare_search_results(self, result, expected)
        _compare_search_results(self, result, expected)

    
    def test_partial_threshold_search(self):
        threshold = 0.1
        N = len(fps)
        result = search.SearchResults(N, fps.arena_ids)
        expected = [[] for i in range(N)]

        for i in range(0, N, 13):
            for j in range(i, N, 8):
                search.partial_threshold_tanimoto_search_symmetric(result, fps, threshold,
                                                                   i, i+13, j, j+8)
                slow_threshold_search(expected, fps, threshold, i, i+13, j, j+8)

        _compare_search_results(self, result, expected)

        counts_before = map(len, result)
        search.fill_lower_triangle(result)
        counts_after = map(len, result)
        self.assertNotEqual(counts_before, counts_after)
        self.assertSequenceEqual(counts_after,
                                 search.count_tanimoto_hits_symmetric(fps, threshold))

        normal = search.threshold_tanimoto_search_symmetric(fps, threshold)
        _compare_search_results(self, result, list(normal.iter_indices_and_scores()))
        
        

class TestKNearest(unittest2.TestCase):
    def test_symmetric(self):
        # query[i] always matches target[i] so x[i] will always contain element[i]
        x = search.knearest_tanimoto_search(fps, fps, 31, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.knearest_tanimoto_search_symmetric(fps, 30, 0.9)
        
        for i, (x_row, y_row) in enumerate(zip(x, y)):
            x_row = x_row.get_indices_and_scores()
            y_row = y_row.get_indices_and_scores()
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()
            self.assertEquals(x_row, y_row, "Problem in %d" % i)

    def test_symmetric2(self):
        # query[i] always matches target[i] so x[i] will always contain element[i]
        x = search.knearest_tanimoto_search(fps, fps, 81, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.knearest_tanimoto_search_symmetric(fps, 80, 0.9)

        for i, (x_row, y_row) in enumerate(zip(x, y)):
            x_row = x_row.get_indices_and_scores()
            y_row = y_row.get_indices_and_scores()
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()
            self.assertEquals(x_row, y_row, "Problem in %d" % i)

if __name__ == "__main__":
    unittest2.main()
