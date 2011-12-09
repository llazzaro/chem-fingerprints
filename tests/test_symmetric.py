# Test the symmetric code

import unittest2

import chemfp
from chemfp import arena

from support import fullpath, PUBCHEM_SDF, PUBCHEM_SDF_GZ

fps = chemfp.load_fingerprints(fullpath("queries.fps"))

class TestCounts(unittest2.TestCase):
    def test_upper_triangle(self):
        # query[i] always matches target[i] so x[i] will be at least one
        x = arena.count_tanimoto_hits_arena(fps, fps, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = arena.count_tanimoto_hits_arena_symmetric(fps, 0.9)

        self.assertEquals(len(x), len(y))
        for i in range(len(x)):
            self.assertEquals(x[i]-1, y[i])


class TestThreshold(unittest2.TestCase):
    def test_upper_triangle(self):
        # query[i] always matches target[i] so x[i] will always contain chei
        x = arena.threshold_tanimoto_search_arena(fps, fps, 0.9)

        # This only processes the upper-triangle, and not the diagonal
        y = arena.threshold_tanimoto_search_arena_symmetric(fps, 0.9)


        rows = list(y)
        # Move elements to the lower triangle
        for rowno, row in enumerate(rows):
            for (colno, score) in row:
                assert colno > rowno, (rowno, colno)
                rows[colno].append( (rowno, score) )

            # Fill in the diagonal
            row.append((fps[rowno][0], 1.0))

            # Put into a consistent order
            row.sort()
            print rowno, row
            
            # Match with the NxM algorithm
            expected_row = x[rowno]
            print expected_row
            expected_row.sort()

            self.assertEquals(row, expected_row, rowno)

if __name__ == "__main__":
    unittest2.main()
