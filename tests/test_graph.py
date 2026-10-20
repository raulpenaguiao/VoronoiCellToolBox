# tests/test_graph.py
# Graph lattices as an independent source of truth: the relevant vectors are the simple cycles
# (known from the graph alone), and the second moment is the tropical moment of de Jong and
# Shokrieh, whose Theorem A computes it from effective resistances without any Voronoi cell.
import itertools
import random
import unittest
from sage.all import QQ, Graph, graphs
from VoronoiCellToolBox import (
    graph_edges, cycle_basis_matrix, riemann_matrix, simple_cycles, graph_relevant_vectors, graph_face_counts,
    tropical_moment, graph_second_moment_function,
    voronoi_cell, relevant_vectors, second_moment,
)


def banana(n):
    return Graph([(0, 1)] * n, multiedges=True)


GRAPHS = {
    "theta": banana(3),
    "banana4": banana(4),
    "K4": graphs.CompleteGraph(4),
    "prism": graphs.CycleGraph(3).cartesian_product(graphs.PathGraph(2)),
    "K33": graphs.CompleteBipartiteGraph(3, 3),
    "loops and parallel edges": Graph([(0, 0), (0, 1), (0, 1), (1, 2), (2, 0), (2, 2)], loops=True, multiedges=True),
}


def random_lengths(G, rng):
    return [QQ(rng.randint(1, 9)) / rng.randint(1, 4) for _ in graph_edges(G)]


class TestCycleBasis(unittest.TestCase):

    def test_basis_of_the_cycle_lattice(self):
        # g rows, each a cycle, unimodular: the maximal minors have gcd 1
        for name, G in GRAPHS.items():
            B = cycle_basis_matrix(G)
            g = len(graph_edges(G)) - G.num_verts() + 1
            self.assertEqual(B.dimensions(), (g, len(graph_edges(G))), name)
            self.assertEqual(B.elementary_divisors(), [1] * g, name)
            for c in simple_cycles(G):                # integer coordinates in the basis
                self.assertTrue(all(x in QQ and x.denominator() == 1 for x in B.solve_left(c)), name)

    def test_simple_cycle_counts(self):
        self.assertEqual(len(simple_cycles(banana(5))), 10)
        self.assertEqual(len(simple_cycles(graphs.CompleteGraph(4))), 7)
        self.assertEqual(len(simple_cycles(graphs.CompleteBipartiteGraph(3, 3))), 15)
        # against Sage's enumeration on the oriented graph (each cycle twice, 2-cycles dropped)
        for G in [graphs.CompleteGraph(5), graphs.PetersenGraph(), GRAPHS["prism"]]:
            directed = [c for c in G.to_directed().all_simple_cycles() if len(c) > 3]
            self.assertEqual(2 * len(simple_cycles(G)), len(directed))

    def test_disconnected_graph(self):
        with self.assertRaises(ValueError):
            cycle_basis_matrix(Graph([(0, 1), (2, 3)]))


class TestRelevantVectorsAreSimpleCycles(unittest.TestCase):

    def test_unit_and_random_lengths(self):
        rng = random.Random(7)
        for name, G in GRAPHS.items():
            predicted = sorted(graph_relevant_vectors(G))
            for L in [None, random_lengths(G, rng), random_lengths(G, rng)]:
                self.assertEqual(predicted, sorted(relevant_vectors(riemann_matrix(G, L))), (name, L))

    def test_sums_of_disjoint_cycles_touch_but_are_not_relevant(self):
        # c1 + c2 with disjoint supports: (c1 + c2)/2 lies on the boundary of the cell (on the
        # face where the facets of c1 and c2 meet), but its bisector is not a facet
        rng = random.Random(8)
        for name in ["banana4", "prism", "loops and parallel edges"]:
            G = GRAPHS[name]
            B = cycle_basis_matrix(G)
            Q = riemann_matrix(G, random_lengths(G, rng))
            V = voronoi_cell(Q)
            relevant = set(relevant_vectors(Q))
            sums = [c1 + c2 for c1, c2 in itertools.combinations(simple_cycles(G), 2)
                    if all(a * b == 0 for a, b in zip(c1, c2))]
            self.assertTrue(sums, name)
            for x in sums:
                a = B.solve_left(x)
                self.assertNotIn(tuple(a), relevant, name)
                self.assertTrue(V.contains(a / 2) and not V.interior_contains(a / 2), name)


class TestFacesAreStronglyConnectedOrientations(unittest.TestCase):
    # Amini: the face poset of the cell is that of the strongly connected orientations of subgraphs

    def test_f_vector_from_the_graph(self):
        rng = random.Random(11)
        for name, G in GRAPHS.items():
            predicted = graph_face_counts(G)
            for L in [None, random_lengths(G, rng)]:
                self.assertEqual(voronoi_cell(riemann_matrix(G, L)).f_vector(), predicted, (name, L))


class TestTropicalMoment(unittest.TestCase):

    def test_circle(self):
        # a circle of length L has I = L/12
        self.assertEqual(tropical_moment(graphs.CycleGraph(4), [1, 2, 3, 4]), QQ(10) / 12)
        self.assertEqual(tropical_moment(Graph([(0, 0)], loops=True), [5]), QQ(5) / 12)

    def test_theorem_A_equals_second_moment(self):
        rng = random.Random(9)
        for name, G in GRAPHS.items():
            for L in [None, random_lengths(G, rng)]:
                Q = riemann_matrix(G, L)
                self.assertEqual(tropical_moment(G, L), second_moment(Q), (name, L))

    def test_known_lattices(self):
        # theta: hexagonal, SM 5/18; K4: A3*, SM 19/32
        self.assertEqual(tropical_moment(banana(3)), QQ(5) / 18)
        self.assertEqual(tropical_moment(graphs.CompleteGraph(4)), QQ(19) / 32)

    def test_cell_combinatorics_do_not_depend_on_lengths(self):
        rng = random.Random(10)
        for name, G in GRAPHS.items():
            f = voronoi_cell(riemann_matrix(G)).f_vector()
            for _ in range(2):
                self.assertEqual(voronoi_cell(riemann_matrix(G, random_lengths(G, rng))).f_vector(), f, name)


class TestChamberFunctionOfLengths(unittest.TestCase):
    # The Macaulay2 chamber pipeline, pulled back to the edge lengths, against Theorem A as
    # rational functions: two unrelated computations of the same function

    def test_against_theorem_A(self):
        for name in ["theta", "banana4", "K4"]:
            G = GRAPHS[name]
            f = graph_second_moment_function(G)
            self.assertEqual(f, tropical_moment(G, f.parent().gens()), name)


if __name__ == "__main__":
    unittest.main()
