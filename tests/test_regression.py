# tests/test_regression.py
# Exact second moments SM(Q) = int_{VC_Q} x^T Q x dx (Euclidean measure, cell volume 1),
# checked against an independent computation (bisectors over {-3..3}^d, Sage's own triangulation).

import contextlib
import io
import unittest
from sage.all import QQ
from VoronoiCellToolBox.voronoi_cell import second_moment, VCell
from VoronoiCellToolBox.orchestration import normalizedChamberSecondMomentPolynomial

KNOWN_SECOND_MOMENTS = [
    ([[2, -1], [-1, 2]], QQ(5) / 18),                              # hexagonal lattice
    ([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]], QQ(19) / 32),        # A3*
    ([[2, -1, -1], [-1, 3, -1], [-1, -1, 3]], QQ(23) / 48),
    ([[QQ(3) / 2, -1, -1], [-1, 3, -1], [-1, -1, 3]], QQ(149) / 384),
    ([[3, -1, -2], [-1, 3, -1], [-2, -1, 3]], QQ(9) / 20),
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], QQ(1) / 4),
]


class TestExactSecondMoment(unittest.TestCase):

    def test_second_moment_is_exact(self):
        for Q, expected in KNOWN_SECOND_MOMENTS:
            with self.subTest(Q=Q):
                self.assertEqual(second_moment(Q, range=2), expected)

    def test_cell_volume_is_one(self):
        for Q, _ in KNOWN_SECOND_MOMENTS:
            with self.subTest(Q=Q):
                self.assertEqual(VCell(Q, range=2).volume(), 1)


class TestVerboseModes(unittest.TestCase):
    """The verbose flags must not change the returned polynomial."""

    def test_verbose_flags_return_same_polynomial(self):
        Q = [[2, -1], [-1, 2]]
        plain = normalizedChamberSecondMomentPolynomial(Q)
        for flag in ["verboseComputationProgress", "verboseSecondMoments", "verboseVertex"]:
            with self.subTest(flag=flag), contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(normalizedChamberSecondMomentPolynomial(Q, **{flag: True}), plain)
                self.assertTrue(out.getvalue())


if __name__ == "__main__":
    unittest.main()
