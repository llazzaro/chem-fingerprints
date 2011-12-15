# Test the symmetric code

import unittest2
from cStringIO import StringIO

import chemfp
from chemfp import search

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

class TestCounts(unittest2.TestCase):
    def test_symmetric(self):
        # query[i] always matches target[i] so x[i] will be at least one
        x = search.count_tanimoto_hits(fps, fps, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.count_tanimoto_hits_symmetric(fps, 0.9)

        self.assertEquals(len(x), len(y))
        for i in range(len(x)):
            self.assertEquals(x[i]-1, y[i])


    def test_zeros(self):
        y = search.count_tanimoto_hits_symmetric(zeros, 0.9)
        self.assertEquals(list(y), [0, 0, 0, 0, 1, 1])
        y = search.count_tanimoto_hits_symmetric(zeros, 0.001)
        self.assertEquals(list(y), [0, 0, 1, 2, 2, 3])
        

class TestThreshold(unittest2.TestCase):
    def test_symmetric(self):
        # query[i] always matches target[i] so x[i] will always contain i
        x = search.threshold_tanimoto_search(fps, fps, 0.9)
        x = list(x.iter_hits())

        # This only processes the upper-triangle, and not the diagonal
        y = search.threshold_tanimoto_search_symmetric(fps, 0.9, include_lower_triangle=False)


        rows = list(y.iter_hits())
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
            expected_row.sort()

            self.assertEquals(row, expected_row, rowno)

    def test_lower_triangle(self):
        # query[i] always matches target[i] so x[i] will always contain i
        x = search.threshold_tanimoto_search(fps, fps, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.threshold_tanimoto_search_symmetric(fps, 0.9)

        for i, (x_row, y_row) in enumerate(zip(x.iter_hits(), y.iter_hits())):
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()

            self.assertEquals(x_row, y_row)

        

class TestKNearest(unittest2.TestCase):
    def test_symmetric(self):
        # query[i] always matches target[i] so x[i] will always contain element[i]
        x = search.knearest_tanimoto_search(fps, fps, 31, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.knearest_tanimoto_search_symmetric(fps, 30, 0.9)

        for i, (x_row, y_row) in enumerate(zip(x.iter_hits(), y.iter_hits())):
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()
            self.assertEquals(x_row, y_row, "Problem in %d" % i)

    def test_symmetric2(self):
        # query[i] always matches target[i] so x[i] will always contain element[i]
        x = search.knearest_tanimoto_search(fps, fps, 81, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = search.knearest_tanimoto_search_symmetric(fps, 80, 0.9)

        for i, (x_row, y_row) in enumerate(zip(x.iter_hits(), y.iter_hits())):
            y_row.append((i, 1.0))
            x_row.sort()
            y_row.sort()
            self.assertEquals(x_row, y_row, "Problem in %d" % i)

if __name__ == "__main__":
    unittest2.main()
