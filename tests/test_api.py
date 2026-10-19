# tests/test_api.py
# Tests of the public API against independent ground truth: an independent second-moment checker,
# random metrics inside secondary cones, the d=4 chambers of Zimmermann's thesis (Section 6.4.3),
# and the known values G(A4*) and G(D4) (Conway-Sloane, Soto).
import random
import unittest
from sage.all import QQ, ZZ, matrix, identity_matrix
from VoronoiCellToolBox import (
    voronoi_cell, relevant_vectors, is_generic, secondary_cone, rayify, matrixify,
    second_moment, normalized_second_moment, chamber_second_moment, metric_coordinates,
)
from tests.independent_checker import sm_independent

HEXAGONAL = [[2, -1], [-1, 2]]
A3STAR = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
A3STAR_NEIGHBOUR = [[QQ(3)/2, -1, -1], [-1, 3, -1], [-1, -1, 3]]


def random_interior_point(Q, rng, scale=4):
    """Q plus a random nonnegative combination of the rays of its secondary cone."""
    d = len(Q)
    M = matrix(QQ, Q)
    for ray in rayify(secondary_cone(Q)):
        M += QQ(rng.randint(0, scale)) * matrix(QQ, matrixify(list(ray), d))
    return [list(row) for row in M.rows()]


def random_generic_metric(d, rng, entries=2):
    while True:
        A = matrix(ZZ, d, d, [rng.randint(-entries, entries) for _ in range(d*d)])
        if A.det() != 0:
            Q = [list(row) for row in A.transpose()*A]
            if is_generic(Q):
                return Q


# The three Delone triangulations in d=4 (Zimmermann (6.7)-(6.9), after Vallentin 4.4.1)
def _E0(i):
    M = matrix(QQ, 4, 4)
    M[i, i] = 1
    return M


def _E(i, j):
    M = _E0(i) + _E0(j)
    M[i, j] = M[j, i] = -1
    return M


EA = matrix(QQ, [[2, 1, -1, -1], [1, 2, -1, -1], [-1, -1, 2, 0], [-1, -1, 0, 2]])   # D4
EB = matrix(QQ, [[1, 1, -1, -1], [1, 1, -1, -1], [-1, -1, 1, 1], [-1, -1, 1, 1]])
_E0S = [_E0(i) for i in range(4)]
_PAIRS = [(i, j) for i in range(4) for j in range(i+1, 4)]
D4_CHAMBERS = {
    "Delta1": ([[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]],
               _E0S + [_E(i, j) for i, j in _PAIRS]),
    "Delta2": ([[5, 1, -2, -2], [1, 5, -2, -2], [-2, -2, 6, -1], [-2, -2, -1, 6]],
               _E0S + [_E(i, j) for i, j in _PAIRS if (i, j) != (0, 1)] + [EA]),
    "Delta3": ([[6, 2, -3, -3], [2, 6, -3, -3], [-3, -3, 6, 1], [-3, -3, 1, 6]],
               _E0S + [_E(0, 2), _E(0, 3), _E(1, 2), _E(1, 3), EA, EB]),
}


class TestCell(unittest.TestCase):

    def test_relevant_vector_counts(self):
        self.assertEqual(len(relevant_vectors(HEXAGONAL)), 6)
        self.assertEqual(len(relevant_vectors(A3STAR)), 14)
        self.assertEqual(len(relevant_vectors([[1, 0], [0, 1]])), 4)

    def test_genericity(self):
        self.assertTrue(is_generic(A3STAR))
        self.assertFalse(is_generic([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        self.assertFalse(is_generic([[5, 3], [3, 2]]))     # a parallelogram cell, equivalent to I

    def test_skewed_metric_needs_no_range(self):
        U = matrix(ZZ, [[7, 2], [3, 1]])
        Q = [list(row) for row in (U.transpose() * matrix(QQ, HEXAGONAL) * U).rows()]
        V = voronoi_cell(Q)
        self.assertEqual((V.volume(), V.n_facets()), (1, 6))
        self.assertEqual(second_moment(Q), QQ(5)/18)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            voronoi_cell([[1, 2], [2, 1]])                  # not positive definite
        with self.assertRaises(ValueError):
            voronoi_cell([[2, 1], [0, 2]])                  # not symmetric


class TestSecondMoment(unittest.TestCase):

    def test_against_independent_checker(self):
        rng = random.Random(1)
        metrics = [HEXAGONAL, A3STAR, A3STAR_NEIGHBOUR] + [random_generic_metric(3, rng) for _ in range(3)]
        for Q in metrics:
            self.assertEqual(second_moment(Q), sm_independent(Q), Q)

    def test_normalized_second_moment(self):
        self.assertEqual(normalized_second_moment(HEXAGONAL), 5/(36*QQ(3).sqrt()))
        self.assertAlmostEqual(float(normalized_second_moment(A3STAR)), 19/(192*2**(1/3)), places=14)
        self.assertEqual(normalized_second_moment([[1, 0], [0, 1]]), QQ(1)/12)

    def test_scale_and_basis_invariance(self):
        U = matrix(ZZ, [[1, 3, 0], [0, 1, 0], [2, 6, 1]])
        Q = [list(row) for row in (U.transpose() * matrix(QQ, A3STAR) * U).rows()]
        self.assertEqual(normalized_second_moment(Q), normalized_second_moment(A3STAR))
        self.assertEqual(normalized_second_moment([[5*x for x in row] for row in A3STAR]),
                         normalized_second_moment(A3STAR))


class TestChamberSecondMoment(unittest.TestCase):

    def test_random_points_in_chamber(self):
        rng = random.Random(2)
        metrics = [HEXAGONAL, A3STAR, A3STAR_NEIGHBOUR] + [random_generic_metric(3, rng) for _ in range(2)]
        for Q in metrics:
            f = chamber_second_moment(Q)
            for _ in range(3):
                Qp = random_interior_point(Q, rng)
                self.assertEqual(f(*metric_coordinates(Qp)), sm_independent(Qp), (Q, Qp))

    def test_continuity_across_a_wall(self):
        # A3* and A3STAR_NEIGHBOUR lie in adjacent chambers; on their common wall, both chamber
        # functions must give the second moment of the (non-generic) wall metric.
        f, g = chamber_second_moment(A3STAR), chamber_second_moment(A3STAR_NEIGHBOUR)
        wall = [[2, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        self.assertFalse(is_generic(wall))
        self.assertEqual(f(*metric_coordinates(wall)), QQ(23)/48)
        self.assertEqual(g(*metric_coordinates(wall)), QQ(23)/48)
        self.assertEqual(sm_independent(wall), QQ(23)/48)

    def test_nongeneric_raises(self):
        with self.assertRaises(ValueError):
            chamber_second_moment([[1, 0], [0, 1]])
        with self.assertRaises(ValueError):
            secondary_cone([[2, -1, -1], [-1, 3, -1], [-1, -1, 3]])


class TestDimension4(unittest.TestCase):

    def test_secondary_cones_match_zimmermann(self):
        from sage.all import Polyhedron
        for name, (Q, rays) in D4_CHAMBERS.items():
            expected = Polyhedron(rays=[[r[i, j] for i in range(4) for j in range(i, 4)] for r in rays], base_ring=QQ)
            self.assertEqual(secondary_cone(Q), expected, name)

    def test_A4star_and_D4(self):
        # G(A4*) from the principal chamber; G(D4) from the boundary ray E_a of chamber Delta2
        Q1 = D4_CHAMBERS["Delta1"][0]
        f1 = chamber_second_moment(Q1)
        self.assertEqual(f1(*metric_coordinates(Q1)), QQ(389)/375)
        self.assertAlmostEqual(float(normalized_second_moment(Q1)), 0.0775588, places=7)
        f2 = chamber_second_moment(D4_CHAMBERS["Delta2"][0])
        D4 = [list(row) for row in EA.rows()]
        self.assertEqual(f2(*metric_coordinates(D4)), QQ(13)/30)
        self.assertEqual(second_moment(D4), QQ(13)/30)
        self.assertAlmostEqual(float(normalized_second_moment(D4)), 0.0766032346, places=9)


if __name__ == '__main__':
    unittest.main()
