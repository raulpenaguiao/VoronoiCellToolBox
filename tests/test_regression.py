# tests/test_regression.py
# Exact second moments SM(Q) = int_{VC_Q} x^T Q x dx (Euclidean measure, cell volume 1),
# checked against an independent computation (bisectors over {-3..3}^d, Sage's own triangulation).

import contextlib
import io
import unittest
from sage.all import QQ, Matrix, PolynomialRing, sage_eval
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


def evaluate_polynomial(polynomial, Q):
    """Evaluate the rational function string returned by Macaulay2 at the metric Q."""
    d = len(Q)
    R = PolynomialRing(QQ, ['q_%d' % k for k in range(d*(d+1)//2)]).fraction_field()
    f = sage_eval(polynomial, locals={str(g): g for g in R.gens()})
    return f(*[QQ(Q[i][j]) for i in range(d) for j in range(i, d)])


class TestAnyChamber(unittest.TestCase):
    """Skewed metrics Q = U^T Q0 U have relevant vectors outside range=2; the polynomial must still be right."""

    def test_skewed_metrics(self):
        U2 = Matrix(QQ, [[7, 2], [3, 1]])
        U3 = Matrix(QQ, [[1, 3, 0], [0, 1, 0], [2, 6, 1]])
        for Q0, U, expected in [([[2, -1], [-1, 2]], U2, QQ(5) / 18),
                                ([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]], U3, QQ(19) / 32)]:
            Q = [list(row) for row in (U.transpose() * Matrix(QQ, Q0) * U).rows()]
            with self.subTest(Q=Q):
                with self.assertRaises(Exception):
                    VCell(Q, range=2)
                self.assertEqual(evaluate_polynomial(normalizedChamberSecondMomentPolynomial(Q), Q), expected)


class TestVerboseModes(unittest.TestCase):
    """The old verbose flags are still accepted, with a deprecation warning, and change nothing;
    debug_m2 prints the raw Macaulay2 output."""

    def test_verbose_flags_return_same_polynomial(self):
        Q = [[2, -1], [-1, 2]]
        plain = normalizedChamberSecondMomentPolynomial(Q)
        for flag in ["verboseComputationProgress", "verboseSecondMoments", "verboseVertex"]:
            with self.subTest(flag=flag), self.assertWarns(DeprecationWarning):
                self.assertEqual(normalizedChamberSecondMomentPolynomial(Q, **{flag: True}), plain)

    def test_debug_m2_prints_raw_output(self):
        Q = [[2, -1], [-1, 2]]
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(normalizedChamberSecondMomentPolynomial(Q, debug_m2=True),
                             normalizedChamberSecondMomentPolynomial(Q))
        self.assertIn("Raw Macaulay2 output", out.getvalue())

    def test_unknown_keyword_raises(self):
        with self.assertRaises(TypeError):
            normalizedChamberSecondMomentPolynomial([[2, -1], [-1, 2]], verbos=True)


if __name__ == "__main__":
    unittest.main()
