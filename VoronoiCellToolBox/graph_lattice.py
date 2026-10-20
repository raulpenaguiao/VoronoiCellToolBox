r"""
Lattices of graphs (tropical Jacobians)
=======================================

A connected finite graph `G` (multiple edges and loops allowed) with positive edge lengths
`\ell` defines the lattice `H_1(G, \ZZ)` of integral cycles, with the inner product
`[x, y] = \sum_e \ell(e) x_e y_e` on `\RR^E`. In the basis given by the rows of a cycle basis
matrix `B` (:func:`cycle_basis_matrix`), its Gram matrix is the Riemann matrix
`Q = B \operatorname{diag}(\ell) B^T` (:func:`riemann_matrix`), and the torus
`H_1(G, \RR) / H_1(G, \ZZ)` is the tropical Jacobian of the metric graph.

The relevant vectors of `Q` are exactly the simple cycles of `G` (Bacher, de la Harpe and
Nagnibeda), whatever the edge lengths: :func:`graph_relevant_vectors` lists them from the
graph alone. A sum `c_1 + c_2` of two edge-disjoint cycles is not relevant: `c_2 - c_1` has the
same length, so the bisector only touches the cell at `(c_1 + c_2)/2`. More precisely, by
Amini's theorem the faces of the cell are the strongly connected orientations of subgraphs of
`G` (:func:`graph_face_counts`).

The second moment of the Voronoi cell of `Q` is the tropical moment `I(\mathrm{Jac}(G))` of
de Jong and Shokrieh (arXiv:1810.02639), which their Theorem A expresses through effective
resistances on `G` (:func:`tropical_moment`). As the combinatorics of the cell depends only on
`G`, the chamber second moment is one rational function of the edge lengths
(:func:`graph_second_moment_function`).

EXAMPLES::

    sage: from VoronoiCellToolBox import riemann_matrix, graph_relevant_vectors, relevant_vectors
    sage: from VoronoiCellToolBox import tropical_moment, second_moment
    sage: K4 = graphs.CompleteGraph(4)
    sage: Q = riemann_matrix(K4); Q
    [ 3  1 -1]
    [ 1  3  1]
    [-1  1  3]
    sage: sorted(graph_relevant_vectors(K4)) == sorted(relevant_vectors(Q))   # 7 cycles
    True
    sage: tropical_moment(K4), second_moment(Q)
    (19/32, 19/32)
"""

from sage.all import QQ, ZZ, Graph, Sequence, diagonal_matrix, matrix, ones_matrix, vector


def graph_edges(G):
    r"""
    Return the edges of `G` in the order used throughout this module: ``G.edges(sort=True)``,
    each edge `(u, v, \text{label})` oriented from `u` to `v`. Edge lengths are given in this
    order.

    EXAMPLES::

        sage: from VoronoiCellToolBox import graph_edges
        sage: graph_edges(graphs.CycleGraph(3))
        [(0, 1, None), (0, 2, None), (1, 2, None)]
    """
    return list(G.edges(sort=True, labels=True))


def _incidence_matrix(G, edges):
    """The signed vertex-edge incidence matrix; a loop gives a zero column."""
    index = {v: i for i, v in enumerate(G.vertices(sort=True))}
    M = matrix(ZZ, len(index), len(edges))
    for j, (u, v, _) in enumerate(edges):
        if u != v:
            M[index[u], j] -= 1
            M[index[v], j] += 1
    return M


def _check_graph(G):
    if not isinstance(G, Graph):
        raise TypeError("G must be an undirected Sage Graph")
    if not G.is_connected():
        raise ValueError("G must be connected")


def cycle_basis_matrix(G):
    r"""
    Return a fundamental cycle basis of `G`: a `g \times |E|` integer matrix whose rows are the
    signed edge vectors of the fundamental cycles of a spanning tree, so that they form a
    `\ZZ`-basis of `H_1(G, \ZZ)`. Here `g = |E| - |V| + 1` is the genus (first Betti number).

    EXAMPLES::

        sage: from VoronoiCellToolBox import cycle_basis_matrix
        sage: cycle_basis_matrix(graphs.CompleteGraph(4))
        [ 1 -1  0  1  0  0]
        [ 1  0 -1  0  1  0]
        [ 0  1 -1  0  0  1]
        sage: G = Graph([(0, 0), (0, 1), (0, 1)], loops=True, multiedges=True)
        sage: cycle_basis_matrix(G)        # a loop and a pair of parallel edges
        [ 1  0  0]
        [ 0 -1  1]
    """
    _check_graph(G)
    edges = graph_edges(G)
    inc = _incidence_matrix(G, edges)
    # a spanning tree by union-find, scanning the edges in order
    parent = {v: v for v in G.vertices(sort=True)}

    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x

    tree, others = [], []
    for j, (u, v, _) in enumerate(edges):
        ru, rv = root(u), root(v)
        if ru != rv:
            parent[ru] = rv
            tree.append(j)
        else:
            others.append(j)
    # the fundamental cycle of e: e itself plus the unique flow on the tree that closes it
    inc_tree = inc.matrix_from_columns(tree)
    rows = []
    for j in others:
        row = [0] * len(edges)
        row[j] = 1
        if edges[j][0] != edges[j][1]:
            for t, x in zip(tree, inc_tree.solve_right(-inc.column(j))):
                row[t] = x
        rows.append(row)
    return matrix(ZZ, len(rows), len(edges), rows)


def riemann_matrix(G, edge_lengths=None, cycle_basis=None):
    r"""
    Return the Riemann matrix `Q = B \operatorname{diag}(\ell) B^T` of `G`: the Gram matrix of
    `H_1(G, \ZZ)` in the basis given by the rows of `B`.

    INPUT:

    - ``G`` -- a connected Sage ``Graph`` (multiple edges and loops allowed)
    - ``edge_lengths`` -- (default: all 1) the lengths, in the order of :func:`graph_edges`;
      they can be numbers or elements of a polynomial ring
    - ``cycle_basis`` -- (default: :func:`cycle_basis_matrix`) the matrix `B`

    EXAMPLES::

        sage: from VoronoiCellToolBox import riemann_matrix
        sage: riemann_matrix(graphs.CycleGraph(3), [1, 2, 3])        # one cycle of length 6
        [6]
        sage: theta = Graph([(0, 1), (0, 1), (0, 1)], multiedges=True)
        sage: riemann_matrix(theta)                                 # the hexagonal lattice
        [2 1]
        [1 2]
        sage: R.<a,b,c> = QQ[]
        sage: riemann_matrix(theta, [a, b, c])
        [a + b     a]
        [    a a + c]
    """
    B = cycle_basis_matrix(G) if cycle_basis is None else matrix(ZZ, cycle_basis)
    lengths = [1] * B.ncols() if edge_lengths is None else list(edge_lengths)
    if len(lengths) != B.ncols():
        raise ValueError("edge_lengths must have one entry per edge")
    K = Sequence(lengths + [QQ(1)]).universe()
    return B.change_ring(K) * diagonal_matrix(K, lengths) * B.transpose()


def simple_cycles(G):
    r"""
    Return the simple cycles of `G` (loops, pairs of parallel edges, and cycles through
    distinct vertices) as signed edge vectors in `\ZZ^E`, one orientation of each.

    EXAMPLES::

        sage: from VoronoiCellToolBox import simple_cycles
        sage: len(simple_cycles(graphs.CompleteGraph(4))), len(simple_cycles(graphs.CompleteGraph(5)))
        (7, 37)
        sage: simple_cycles(Graph([(0, 1), (0, 1), (0, 1)], multiedges=True))
        [(1, -1, 0), (1, 0, -1), (0, 1, -1)]
    """
    edges = graph_edges(G)
    order = {v: i for i, v in enumerate(G.vertices(sort=True))}
    incident = {v: [] for v in order}
    for j, (u, v, _) in enumerate(edges):
        if u != v:
            incident[u].append((j, v))
            incident[v].append((j, u))

    cycles, seen = [], set()

    def record(walk):
        support = frozenset(j for j, _, _ in walk)
        if support not in seen:
            seen.add(support)
            x = [0] * len(edges)
            for j, a, b in walk:
                x[j] += 1 if (a, b) == (edges[j][0], edges[j][1]) else -1
            cycles.append(vector(ZZ, x))

    for j, (u, v, _) in enumerate(edges):
        if u == v:
            record([(j, u, v)])

    # each cycle is found from its smallest vertex, once in each direction
    def extend(start, here, used, visited, walk):
        for j, w in incident[here]:
            if j in used:
                continue
            if w == start:
                record(walk + [(j, here, w)])
            elif w not in visited and order[w] > order[start]:
                extend(start, w, used | {j}, visited | {w}, walk + [(j, here, w)])

    for s in sorted(order, key=order.get):
        extend(s, s, frozenset(), {s}, [])
    return cycles


def graph_relevant_vectors(G, cycle_basis=None):
    r"""
    Return the relevant vectors of the Riemann matrix of `G` predicted by its simple cycles,
    in the coordinates of the cycle basis: the `a \in \ZZ^g` such that `a B` is `\pm` a simple
    cycle. They do not depend on the edge lengths.

    EXAMPLES::

        sage: from VoronoiCellToolBox import graph_relevant_vectors, riemann_matrix, relevant_vectors
        sage: G = graphs.CompleteBipartiteGraph(3, 3)
        sage: P = graph_relevant_vectors(G); len(P)          # 9 four-cycles and 6 six-cycles, two signs
        30
        sage: Q = riemann_matrix(G, [1, 2, 3, 4, 5, 6, 7, 8, 9])
        sage: sorted(P) == sorted(relevant_vectors(Q))
        True
    """
    B = cycle_basis_matrix(G) if cycle_basis is None else matrix(ZZ, cycle_basis)
    found = []
    for c in simple_cycles(G):
        a = B.solve_left(c)
        found += [tuple(a), tuple(-a)]
    return found


def graph_face_counts(G):
    r"""
    Return the f-vector of the Voronoi cell of `H_1(G, \ZZ)` predicted by the graph alone, in
    the format of ``Polyhedron.f_vector()``: `(1, f_0, \dots, f_{g-1}, 1)`.

    By Amini's theorem (arXiv:1007.2456, conjectured by Caporaso and Viviani), the faces of the
    cell correspond to the strongly connected orientations of subgraphs of `G` (each component
    strongly connected, isolated vertices ignored), a subgraph of genus `k` giving a face of
    codimension `k`. Whatever the edge lengths, the cell has this f-vector. All `3^{|E|}`
    choices (edge absent, or one of its two orientations) are enumerated, so this is for small
    graphs.

    EXAMPLES::

        sage: from VoronoiCellToolBox import graph_face_counts, voronoi_cell, riemann_matrix
        sage: K4 = graphs.CompleteGraph(4)
        sage: graph_face_counts(K4)                                 # the truncated octahedron
        (1, 24, 36, 14, 1)
        sage: graph_face_counts(K4) == voronoi_cell(riemann_matrix(K4, [1, 2, 3, 4, 5, 6])).f_vector()
        True
    """
    from itertools import product
    from sage.all import DiGraph
    _check_graph(G)
    edges = graph_edges(G)
    g = len(edges) - G.num_verts() + 1
    codim = [0] * (g + 1)
    for choice in product((0, 1, -1), repeat=len(edges)):
        arcs = [(u, v) if s == 1 else (v, u) for (u, v, _), s in zip(edges, choice) if s]
        if not arcs:
            codim[0] += 1
            continue
        D = DiGraph(arcs, multiedges=True, loops=True)
        component = {v: i for i, C in enumerate(D.strongly_connected_components()) for v in C}
        if all(component[a] == component[b] for a, b in arcs):
            codim[len(arcs) - D.num_verts() + D.connected_components_number()] += 1
    return vector(ZZ, [1] + [codim[g - k] for k in range(g)] + [1])


def tropical_moment(G, edge_lengths=None):
    r"""
    Return the tropical moment `I(\mathrm{Jac}(G)) = \int_{\mathrm{Vor}(0)} [z, z] \, d\mu`
    by the formula of de Jong and Shokrieh, Theorem A of arXiv:1810.02639: for any vertex `q`,

    .. MATH::

        I = \frac{1}{12} \sum_e F(e)^2 \ell(e)
            + \frac{1}{4} \sum_{e = \{u, v\}} \Big( r(u, v)
              - \frac{j_u(v, q)^2 + j_v(u, q)^2}{\ell(e)} \Big),

    with the effective resistance `r`, the potentials `j` and the Foster coefficients
    `F(e) = 1 - r(u, v)/\ell(e)` of the electrical network with resistances `\ell`, all read
    off the pseudoinverse of the weighted Laplacian. No Voronoi cell is computed; the result
    equals ``second_moment(riemann_matrix(G, edge_lengths))``, which makes this an
    independent check of the package.

    INPUT: as for :func:`riemann_matrix`. With symbolic edge lengths the result is a rational
    function of them.

    EXAMPLES::

        sage: from VoronoiCellToolBox import tropical_moment, riemann_matrix, second_moment
        sage: tropical_moment(graphs.CycleGraph(5), [1, 1, 1, 1, 2])   # a circle of length L: L/12
        1/2
        sage: G = graphs.CompleteBipartiteGraph(3, 3); L = [1, 2, 3, 1, 2, 3, 3, 2, 1]
        sage: tropical_moment(G, L) == second_moment(riemann_matrix(G, L))
        True
        sage: R.<a,b,c> = QQ[]
        sage: tropical_moment(Graph([(0, 1), (0, 1), (0, 1)], multiedges=True), [a, b, c])
        (1/12*a^2*b + 1/12*a*b^2 + 1/12*a^2*c + 1/3*a*b*c + 1/12*b^2*c + 1/12*a*c^2 + 1/12*b*c^2)/(a*b + a*c + b*c)
    """
    _check_graph(G)
    edges = graph_edges(G)
    lengths = [1] * len(edges) if edge_lengths is None else list(edge_lengths)
    if len(lengths) != len(edges):
        raise ValueError("edge_lengths must have one entry per edge")
    K = Sequence(lengths + [QQ(1)]).universe().fraction_field()
    lengths = [K(x) for x in lengths]
    index = {v: i for i, v in enumerate(G.vertices(sort=True))}
    n = len(index)
    lap = matrix(K, n, n)
    for (u, v, _), l in zip(edges, lengths):
        if u != v:
            i, j = index[u], index[v]
            lap[i, i] += 1/l
            lap[j, j] += 1/l
            lap[i, j] -= 1/l
            lap[j, i] -= 1/l
    J = ones_matrix(K, n, n) / n
    Gp = (lap + J).inverse() - J            # the pseudoinverse of the Laplacian

    def j(z, x, y):                          # potential at x, unit current in at y, out at z
        return Gp[x, y] - Gp[x, z] - Gp[z, y] + Gp[z, z]

    q = 0
    total = K(0)
    for (u, v, _), l in zip(edges, lengths):
        u, v = index[u], index[v]
        r = j(v, u, u)
        F = 1 - r/l
        total += F**2 * l / 12 + (r - (j(u, v, q)**2 + j(v, u, q)**2) / l) / 4
    return _monic_denominator(total)


def _monic_denominator(f):
    """Write a rational function with a monic denominator (Sage keeps stray constants)."""
    if not hasattr(f, "denominator") or f.parent() is QQ:
        return f
    c = f.denominator().lc()
    return f.parent()(f.numerator() / c) / f.parent()(f.denominator() / c)


def graph_second_moment_function(G, cycle_basis=None):
    r"""
    Return the second moment of the Voronoi cell of `H_1(G, \ZZ)` as a rational function of
    the edge lengths, computed by the chamber pipeline (:func:`chamber_second_moment`) and
    pulled back along `\ell \mapsto B \operatorname{diag}(\ell) B^T`.

    The Riemann matrices of `G` are usually not generic: they lie on a face of a secondary cone
    that does not depend on the lengths. The face formula is computed at unit lengths and is
    valid for all positive lengths. It agrees with de Jong and Shokrieh's formula
    (:func:`tropical_moment`) as a rational function.

    OUTPUT: an element of the fraction field of `\QQ[\ell_0, \dots, \ell_{|E|-1}]`

    EXAMPLES::

        sage: from VoronoiCellToolBox import graph_second_moment_function, tropical_moment
        sage: G = graphs.CompleteGraph(4)
        sage: f = graph_second_moment_function(G)
        sage: f == tropical_moment(G, f.parent().gens())
        True
    """
    from sage.all import PolynomialRing
    from VoronoiCellToolBox.orchestration import chamber_second_moment
    B = cycle_basis_matrix(G) if cycle_basis is None else matrix(ZZ, cycle_basis)
    ring = PolynomialRing(QQ, ['l_%d' % k for k in range(B.ncols())])
    Q1 = riemann_matrix(G, cycle_basis=B)
    f = chamber_second_moment(Q1, allow_nongeneric=True)
    Ql = riemann_matrix(G, ring.gens(), cycle_basis=B)
    g = Ql.nrows()
    point = [Ql[i, k] for i in range(g) for k in range(i, g)]   # as in metric_coordinates
    return _monic_denominator(ring.fraction_field()(f.numerator()(*point)) / f.denominator()(*point))
