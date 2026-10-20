r"""
Voronoi cells of the lattice `\ZZ^d` with a metric `Q`
======================================================

A positive definite symmetric matrix `Q` (a Gram matrix) defines the inner product
`\langle x, y \rangle_Q = x^T Q y` on `\RR^d`. This module computes, exactly over `\QQ`:

- the Voronoi cell `V_Q = \{x : x^T Q x \le (x-v)^T Q (x-v) \text{ for all } v \in \ZZ^d\}`
  (:func:`voronoi_cell`), its relevant vectors (:func:`relevant_vectors`), and whether it is
  primitive, that is, whether `Q` lies in the interior of a secondary cone (:func:`is_generic`);
- the secondary cone of `Q`: the metrics whose Voronoi cell has the same combinatorics
  (:func:`secondary_cone`);
- a pulling triangulation (:func:`pulling_triangulation`) and the second moment
  `SM(Q) = \int_{V_Q} x^T Q x \, dx` (:func:`second_moment`), with Euclidean measure, so the
  cell has volume 1, and the dimensionless `G(Q) = SM(Q) / (d \det(Q)^{1/d})`
  (:func:`normalized_second_moment`).

Every cell is certified: :func:`VCell` checks that the computed polytope has volume 1, which
fails if a relevant vector was missed.

EXAMPLES::

    sage: from VoronoiCellToolBox import voronoi_cell, second_moment, normalized_second_moment
    sage: Q = [[2, -1], [-1, 2]]                     # the hexagonal lattice
    sage: voronoi_cell(Q).n_vertices()
    6
    sage: second_moment(Q)
    5/18
    sage: normalized_second_moment(Q)
    5/108*sqrt(3)
"""

import math
import numpy
import itertools
from sage.all import Polyhedron, vector, dimension, Matrix, Rational, QQ, ZZ, SR, gcd as sage_gcd
from sage.arith.functions import LCM_list


def _rows(Q):
    """The metric as a list of rows, so that a Sage matrix and a list of lists are accepted alike."""
    return [list(row) for row in Q.rows()] if hasattr(Q, "rows") else Q


def _candidates(d, rangeorlist, default_range=None):
    r"""
    The list of candidate relevant vectors given by ``range=r`` (all of `\{-r..r\}^d \setminus 0`)
    or ``list=[...]``. ``default_range`` is used when neither is given.
    """
    if "range" in rangeorlist or ("list" not in rangeorlist and default_range is not None):
        r = rangeorlist.get("range", default_range)
        if (r <= 0):
            raise ValueError("The range must be a positive integer.")
        p = [list(vec) for vec in itertools.product(range(-r, r+1), repeat=d)]
        p.remove([0] * d)
        return p
    if "list" in rangeorlist:
        return rangeorlist["list"]
    raise ValueError("A range or list of potential relevant vectors needs to be given.")

def VCell(Q, **rangeorlist):
    r"""
    Computes the Voronoi cell of the Z^d lattice with metric `Q`

    This condition is crucial for validating that the matrix `Q` represents
    a valid transformation in the context of lattice quantizers. Specifically,
    the determinant of `Q` must be greater than zero to ensure that the
    transformation preserves orientation and does not collapse the lattice.

    Parameters:
    - `Q`: A 2D array-like object (e.g., list of lists, NumPy array, or Sage matrix)
      representing the matrix. The matrix should be square (i.e., have the same number of rows and columns).
    - input2 (range = int or list = list): 
        If range, an integer, it must be a positive value representing a range. 
        If a list, it should contain vertices for processing.
    
    Usage:
    - `VCell([[2, -1], [-1, 2]], range=2)` or `VCell([[2, -1], [-1, 2]], list=[...])`

    Behavior:
    - If `Q` is provided as a list of lists, it will be converted to a NumPy array
      for determinant computation.
    - If `Q` is already a NumPy array, it will be used directly.
    - If the determinant of `Q` is less than or equal to zero, this indicates
      an invalid or degenerate transformation, and appropriate handling should
      be implemented in the surrounding code.

    Note:
    - Ensure that `Q` is a valid square matrix before calling this function.
    - Symmetry and positive definiteness are checked exactly; a ``ValueError`` is raised otherwise.
    - Raises an exception with message "relevant vector list not complete" when the polytope
      cut out by the candidates does not have volume 1. :func:`voronoi_cell` chooses the
      candidates automatically.

    EXAMPLES::

        sage: from VoronoiCellToolBox import VCell
        sage: sorted(VCell([[2, -1], [-1, 2]], range=2).vertices_list())
        [[-2/3, -1/3], [-1/3, -2/3], [-1/3, 1/3], [1/3, -1/3], [1/3, 2/3], [2/3, 1/3]]
        sage: VCell([[61, 11], [11, 2]], range=2)
        Traceback (most recent call last):
        ...
        Exception: relevant vector list not complete
    """
    Q = _rows(Q)
    M = Matrix(QQ, Q)
    if not M.is_symmetric():
        raise ValueError(f"Matrix {Q} is not symmetric.")
    if not M.is_positive_definite():
        raise ValueError(f"Matrix {Q} is not positive-definite.")
    d = len(Q)
    p = _candidates(d, rangeorlist)
    ineqs = []
    for vert in p:
        #Collect an inequality for each potential relevant vector (not all may be useful)
        #equation is -2 v^T.Q.x + v^T . Q . v >= 0
        b = 0
        for i in range(d):
            for j in range(d):
                b += vert[i]*Q[i][j]*vert[j]
        vec = [-2*sum([vert[i]*Q[i][j] for i in range(d)]) for j in range(d)]
        ineqs += [tuple([b] + vec)]
    VC = Polyhedron(ieqs = ineqs)
    if VC.volume() != 1:
        raise Exception("relevant vector list not complete")
    return VC


def reduce_metric(Q):
    r"""
    LLL-reduces the metric `Q` by a unimodular change of basis U, so that Qr = U^T Q U.

    The Voronoi cell of Qr is U^{-1} times the Voronoi cell of Q, so both have the
    same second moment, and Q' -> U^T Q' U maps the secondary cone of Q onto that of Qr.
    The relevant vectors of a reduced metric have small entries, which is what makes
    `VCell(Qr, range=2)` succeed where `VCell(Q, range=2)` misses relevant vectors.

    Parameters:
        Q (Matrix): The positive definite metric matrix.

    Returns:
        tuple: (Qr, U), with Qr a list of lists and U a unimodular Sage matrix.

    EXAMPLES::

        sage: from VoronoiCellToolBox import reduce_metric
        sage: Qr, U = reduce_metric([[5, 3], [3, 2]])
        sage: Qr
        [[1, 0], [0, 1]]
        sage: U
        [-1  1]
        [ 1 -2]
    """
    Q = _rows(Q)
    U = Matrix(QQ, Q).LLL_gram()
    Qr = U.transpose() * Matrix(QQ, Q) * U
    return [list(row) for row in Qr.rows()], U


def voronoi_cell(Q, candidate_range=None):
    r"""
    Return the Voronoi cell of `\ZZ^d` with metric `Q`, as an exact Sage polyhedron.

    The candidate relevant vectors are chosen automatically: `Q` is LLL-reduced to
    `Q_r = U^T Q U` (:func:`reduce_metric`), the cell of `Q_r` is computed from the candidates
    `\{-r..r\}^d` with `r = 2, 3, \dots` until it has volume 1, and it is mapped back by
    `V_Q = U V_{Q_r}`. The volume-1 check certifies that no relevant vector was missed.

    INPUT:

    - ``Q`` -- a positive definite symmetric matrix (list of lists or Sage matrix)
    - ``candidate_range`` -- (default: ``None``) if given, use exactly the candidates
      `\{-r..r\}^d` for `Q` itself, as ``VCell(Q, range=candidate_range)`` does

    OUTPUT: a ``Polyhedron`` over `\QQ`

    EXAMPLES::

        sage: from VoronoiCellToolBox import voronoi_cell
        sage: voronoi_cell([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]).f_vector()  # A_3^*: truncated octahedron
        (1, 24, 36, 14, 1)
        sage: V = voronoi_cell([[61, 11], [11, 2]])   # relevant vectors far outside {-2..2}^2
        sage: V.volume(), V.n_facets()
        (1, 4)
    """
    Q = _rows(Q)
    if candidate_range is not None:
        return VCell(Q, range=candidate_range)
    Qr, U = reduce_metric(Q)
    r = 2
    while True:
        try:
            Vr = VCell(Qr, range=r)
            break
        except Exception as e:
            if "relevant vector list not complete" not in str(e):
                raise
            r += 1
    if U == 1:
        return Vr
    return Polyhedron(vertices=[U * vector(QQ, v) for v in Vr.vertices()], base_ring=QQ)


def relevant_vectors(Q):
    r"""
    Return the relevant vectors of `Q`: the `v \in \ZZ^d` whose bisector
    `\{x : x^T Q v = v^T Q v / 2\}` contains a facet of the Voronoi cell.

    There are `2(2^d - 1)` of them exactly when `Q` is generic (:func:`is_generic`).

    EXAMPLES::

        sage: from VoronoiCellToolBox import relevant_vectors
        sage: sorted(relevant_vectors([[2, -1], [-1, 2]]))
        [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, 0), (1, 1)]
        sage: len(relevant_vectors([[1, 0], [0, 1]]))       # the square lattice is not generic
        4
    """
    Q = _rows(Q)
    VC = voronoi_cell(Q)
    return [tuple(relevant_vector(Q, F)) for F in VC.facets()]


def is_generic(Q):
    r"""
    Return whether `Q` lies in the interior of a secondary cone, that is, whether its Voronoi
    cell is primitive (has the maximal number `2(2^d - 1)` of facets).

    Only for generic `Q` is the chamber second moment one rational function on a full
    neighbourhood of `Q`. (Testing that every vertex lies on exactly `d` facets is not enough
    in general: the square lattice passes that test.)

    EXAMPLES::

        sage: from VoronoiCellToolBox import is_generic
        sage: is_generic([[2, -1], [-1, 2]]), is_generic([[1, 0], [0, 1]])
        (True, False)
        sage: is_generic([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]])
        True
    """
    Q = _rows(Q)
    d = len(Q)
    return voronoi_cell(Q).n_facets() == 2 * (2**d - 1)


def FromMatrixPositionToListPosition(i: int, j: int, d: int) -> int:
    r"""
    Maps the position (i, j) on the $d \times d$ matrix to a list index, so that
    each symmetric entry can be encoded in a list of size $\frac{d\times (d+1)}{2}$.

    Parameters:
        i (int): The row index of the matrix.
        j (int): The column index of the matrix.
        d (int): The dimension of the matrix.

    Usage:
        FromMatrixPositionToListPosition(0, 1, 3)  # Returns 1
        FromMatrixPositionToListPosition(1, 2, 3)  # Returns 4
        FromMatrixPositionToListPosition(2, 2, 3)  # Returns 5

    Behavior:
        The function should return a unique index for each pair (i, j) such that
        all pairs (i, j) with i < j are assigned indices before any pairs with
        i = j.

    Note:
        Also works when i > j. Undefined when $j < 0$, $i < 0$, $i \geq d$ or $j \geq d$.
    """
    if i>j: return FromMatrixPositionToListPosition(j, i, d)
    return i*d+j-(((i+1)*i)//2)



def matrixify(list_of_ute, d):
    r"""
    Converts a list of upper triangular entries into a $d \times d$ matrix.

    Parameters:
        list_of_ute (list): The list of entries of the symmetric matrix.
        d (int): The dimension of the matrix.
    
    Returns:
        list: A $d \times d$ matrix represented as a list of lists.
    
    EXAMPLES::

        sage: from VoronoiCellToolBox import matrixify
        sage: matrixify([3, 1, 5, 3, -1, 2], 3)
        [[3, 1, 5], [1, 3, -1], [5, -1, 2]]
    """
    return [[list_of_ute[int(FromMatrixPositionToListPosition(i, j, d))]  for j in range(d)] for i in range(d)]

#calculates equation of wall in secondary cone
def eq_of_wall(neighborVectors: list[vector], isolatedVector: vector, d: int, verbose: bool = False) -> list:
    r"""
    Finds the linear inequality that a matrix must satisfy in order for a
    particular potential relevant vector nonRelevantVectors to not be relevant due to the 
    conditions imposed by the list of relevant vectors neighborVectors
    Computes the equation of the wall in the secondary cone.

    Parameters:
        neighborVectors (list(vector)): The list of neighbor vectors.
        isolatedVector (vector): The isolated vector.
        d (int): The dimension of the matrix.
        verbose (bool): Whether to print verbose output.

    Returns:
        list: The coefficients of the inequality on the matrix coefficients
        defining the wall that decides if isolatedVector is relevant.

    EXAMPLES::

        sage: from VoronoiCellToolBox import eq_of_wall
        sage: eq_of_wall([(0, 1, 1), (1, 1, 1), (0, 0, 1)], [-1, 0, 0], 3)
        [2, 2, 2, 0, 0, 0]

    """
    if verbose:
        print(neighborVectors, isolatedVector, d)
    Tinv = Matrix(neighborVectors).inverse()
    inequalityCoefficients = [0]*(d*(d+1)//2)
    for i in range(0, d):
        for j in range(0, d):
            inequalityCoefficients[FromMatrixPositionToListPosition(i, j, d)] += isolatedVector[i]*isolatedVector[j]
            if(verbose):
                print(isolatedVector[i]*isolatedVector[j], " - > (", i, ", ", j, ")")
    if(verbose): 
        print(Tinv)
    for i in range(0, d):
        for j in range(0, d):
            for k in range(0, d):
                for l in range(0, d):
                    inequalityCoefficients[FromMatrixPositionToListPosition(k, l, d)] -= isolatedVector[i]*Tinv[i][j]*neighborVectors[j][k]*neighborVectors[j][l]
                    if ( verbose): 
                        print( isolatedVector[i]*Tinv[i][j]*neighborVectors[j][k]*neighborVectors[j][l], " -> (", k,", ", l,")")
                        print(" (", i, " , ", j, ", ", k ,", ", l, ") - ", isolatedVector[i], ", ", Tinv[i][j], ", ", neighborVectors[j][k], ", ", neighborVectors[j][l], " -> (", k, ", ", l, ") ")
    return inequalityCoefficients #do not include b becuase it's homogeneous and b is always 0
#\sum a_{i, j} * q_{i, j} <= 0 for newv to NOT be relevant   

def lcmReduce(v):
    r"""
    Reduces a vector by finding the least common multiple of its denominators
    and scaling the vector accordingly.

    Parameters:
        v (list): The vector to be reduced.

    Returns:
        vector: The reduced vector.

    Usage:
        lcmReduce([1/2, 2/3, 3/4]) #Returns [6, 8, 9]
    """
    try:
        lcm =  LCM_list([vi.denominator() for vi in v])
    except:
        return lcmReduce([Rational(int(round(float(vi)*10**9)),int(10**9)) for vi in v])
    return vector(gcdReduce([vi*lcm for vi in v]))

def gcdReduce(v):
    r"""
    Reduces a vector of integers by finding the greatest common divisor of its components
    and scaling the vector accordingly.

    Parameters:
        v (list): The vector to be reduced.

    Returns:
        vector: The reduced vector.

    Usage:
        gcdReduce([12, 24, 34]) #Returns [6, 12, 17]
    """
    gcd = GCD_list(v)
    return [vi/gcd for vi in v]


def gcd(a, b):
    r"""
    Computes the greatest common divisor of two integers.

    Parameters:
        a (int): The first integer.
        b (int): The second integer.

    Returns:
        int: The greatest common divisor of a and b.

    Usage:
        gcd(48, 18) #Returns 6
    """
    if (a < 0): return gcd(-a, b)
    if (a > b): return gcd(b, a)
    if (a < 1 and a > -1): return b
    return gcd(b%a, a)

def GCD_list(v):
    r"""
    Computes the greatest common divisor of a list of integers.

    Parameters:
        v (list): The list of integers.

    Returns:
        int: The greatest common divisor of the integers in the list.

    Usage:
        GCD_list([102, 304, 54, 80]) #Returns 2
    """
    if (len(v) == 0): return 1
    if (len(v) == 1): return v[0]
    if (len(v) == 2): return gcd(v[0], v[1])
    return GCD_list([gcd(v[0], v[1])] + v[2:] )

def relevant_vector(Q, F):
    r"""
    For a facet F of the voronoi cell V, this function computes 
    the corresponding relevant vector.

    Parameters:
        Q (Matrix): The metric matrix defining the Voronoi cell.
        F (Facet): The facet for which to compute the relevant vector.

    Returns:
        Vector: The corresponding relevant vector.

    Usage:
        #The following code prints the relevant vectors of each facet for d = 2
        Q = [[2, -1], [-1, 2]]
        VC = VCell(Q, range = 2)
        list_facets = VC.facets()
        for F in list_facets:
            print(F.vertices(), " - ", relevant_vector(Q, F))

    Notes:
        We assume that the facet F is a face of the Voronoi cell V, 
        defined by the metric matrix Q.
        The behaviour of this function is not defined when F is not a facet of V.
    """
    # The facet is b + a.x >= 0 with (a, b) = c (-2 Q v, v^T Q v) for some unknown c > 0.
    # So w = -Q^{-1} a = 2 c v, and b = w^T Q w / (4 c) gives c, hence v = 2 b w / (w^T Q w).
    Qm = Matrix(QQ, Q)
    ieq = F.ambient_Hrepresentation()[0]
    b, a = QQ(ieq[0]), vector(QQ, ieq[1:])
    w = -Qm.inverse() * a
    return vector(ZZ, 2 * b * w / (w * Qm * w))

def Qform(v, Q):
    r"""
    For a matrix $Q$ defining a metric on $Z^d$, and a vector v, this function 
    computes the quadratic form $v^T Q v$.
    Parameters:
        v (Vector): The vector to be evaluated.
        Q (Matrix): The metric matrix.
    Returns:
        Scalar: The value of the quadratic form.
    EXAMPLES::

        sage: from VoronoiCellToolBox import Qform
        sage: Qform(vector([1, 2]), Matrix([[2, 0], [0, 2]]))
        10
    """
    return (Matrix([[vi] for vi in v]).transpose()*Matrix(Q)*Matrix([v]).transpose())[0][0]


def equidistantPoint(vlist, Q):
    r"""
    Given a list of non-zero lattice points, 
    this function computes the equidistant point to all points and zero
    Parameters:
        vlist (list): A list of non-zero lattice points.
        Q (Matrix): The metric matrix.

    Returns:
        Vector: The equidistant point to all points and zero.
    Usage:
        equidistantPoint([vector([1, 0]), vector([0, 1])], Matrix([[1, 0], [0, 1]])) #Returns [1/2, 1/2]
    
    Errors:
        vlist is expected to be a square matrix
    """
    T = Matrix([list(v) for v in vlist])
    return (2*T*Matrix(Q)).inverse()*vector([Qform(v, Q) for v in vlist])
    
def relevantVectorDictionary(P, Q):
    r"""
    Computes a dictionary, exposing the relevant vectors associated with each vertex of P.
    For each vertex v, the relevant vector relative to each facet incident to v is computed.

    Parameters:
        P (Polyhedron): The polyhedron to be analyzed.
        Q (Matrix): The metric matrix.

    Returns:
        dict: A dictionary mapping each vertex of P to its incident relevant vectors.

    Usage:
        Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        P = VCell(Q, range = 2)
        len(relevantVectorDictionary(P, Q))# Returns 24
    """
    vdict = {}
    for facet in P.facets():
        r_vec = relevant_vector(Q, facet)
        for v in facet.vertices():
            if v in vdict:
                vdict[v] += [r_vec]
            else:
                vdict[v] = [r_vec]
    return vdict


def secondary_cone(Q, verbose = False, **rangeorlist):
    r"""
    Return the secondary cone of a generic metric `Q`: the closed cone of metrics `Q'` whose
    Voronoi cell has the same relevant vectors and the same vertex structure as that of `Q`.

    The cone lives in the space of symmetric matrices with coordinates the upper-triangular
    entries in row-major order, `(q_{11}, q_{12}, \dots, q_{1d}, q_{22}, \dots, q_{dd})`, the
    variables `q_0, q_1, \dots` of :func:`chamber_second_moment`; :func:`matrixify` converts a
    ray back to a matrix. For each vertex of the cell (determined by its `d` relevant vectors)
    and each candidate vector, the wall equation of :func:`eq_of_wall` must stay nonnegative.

    INPUT:

    - ``Q`` -- a generic positive definite metric
    - ``range=r`` or ``list=[...]`` -- (optional) the candidate relevant vectors for `Q` itself.
      If neither is given, `Q` is LLL-reduced, the cone of the reduced metric is computed with
      ``range=2``, and it is mapped back (`Q' \mapsto U^{-T} Q' U^{-1}`).

    OUTPUT: a ``Polyhedron`` (a pointed cone over `\QQ`)

    A ``ValueError`` is raised when `Q` is not generic, since its secondary cone is then not
    full-dimensional.

    EXAMPLES::

        sage: from VoronoiCellToolBox import secondary_cone, rayify, matrixify
        sage: C = secondary_cone([[2, -1], [-1, 2]])
        sage: sorted(matrixify(list(r), 2) for r in rayify(C))
        [[[0, 0], [0, 1]], [[1, -1], [-1, 1]], [[1, 0], [0, 0]]]
        sage: secondary_cone([[1, 0], [0, 1]])
        Traceback (most recent call last):
        ...
        ValueError: Q is not generic: its Voronoi cell has 4 facets, not 2(2^d - 1) = 6
    """
    Q = _rows(Q)
    d = len(Q)
    if not rangeorlist:
        Qr, U = reduce_metric(Q)
        Cr = secondary_cone(Qr, verbose, range=2)
        if U == 1:
            return Cr
        Uinv = U.inverse()
        rays = []
        for ray in rayify(Cr):
            X = Uinv.transpose() * Matrix(QQ, matrixify(list(ray), d)) * Uinv
            rays.append([X[i, j] for i in range(d) for j in range(i, d)])
        return Polyhedron(rays=rays, base_ring=QQ)

    prv = _candidates(d, rangeorlist)
    VC = VCell(Q, list=prv)
    if VC.n_facets() != 2 * (2**d - 1):
        raise ValueError(f"Q is not generic: its Voronoi cell has {VC.n_facets()} facets, not 2(2^d - 1) = {2 * (2**d - 1)}")
    ifs = relevantVectorDictionary(VC, Q)
    ineqs = set()
    for vert, opposite in _regulator_points(ifs, d).items():
        neighbors = ifs[vert]
        B = Matrix(QQ, neighbors)
        Binv = B.inverse()
        outer = [B.row(j).column() * B.row(j).row() for j in range(d)]
        for vec in opposite:
            # Same functional as eq_of_wall(neighbors, vec, d), computed with one inverse per vertex
            x = vector(QQ, vec)
            c = x * Binv
            M = x.column() * x.row() - sum(c[j] * outer[j] for j in range(d))
            coeffs = vector(QQ, [M[i, j] if i == j else 2 * M[i, j] for i in range(d) for j in range(i, d)])
            if coeffs == 0:
                continue
            coeffs = coeffs * coeffs.denominator()
            coeffs = coeffs / sage_gcd(list(coeffs))
            ineqs.add(tuple(coeffs))
    if verbose:
        print(len(ifs), "vertices,", len(ineqs), "distinct wall inequalities")
    return Polyhedron(ieqs=[[0] + list(ieq) for ieq in ineqs], base_ring=QQ)


def _regulator_points(ifs, d):
    r"""
    For each vertex `v` of a primitive Voronoi cell, with relevant vectors `N(v)`, return the
    lattice points opposite to the `d+1` facets of the Delaunay simplex `S_v = \{0\} \cup N(v)`,
    that is, the points `p'` such that `S'` = (facet of `S_v`) `\cup \{p'\}` is the adjacent
    Delaunay simplex. By Voronoi's theory the secondary cone is cut out by the "regulator"
    inequalities "`p'` is not inside the empty sphere of `S_v`" (:func:`eq_of_wall`) over these
    pairs alone.

    - Facets containing `0` are shared with `S_w` for the vertices `w` adjacent to `v` along an
      edge of the cell: `p'` is the vector of `N(w)` not in `N(v)`.
    - The facet `N(v)` opposite `0`: translating by `-n` for some `n \in N(v)`, the adjacent
      simplex is `S_w` for the vertex `w \ne v - n` with `N(w) \supset \{m - n : m \in N(v), m \ne n\}`,
      and `p' = n + e` for the remaining vector `e` of `N(w)`.
    """
    sets = {vert: {tuple(x) for x in N} for vert, N in ifs.items()}
    result = {}
    for vert, N in sets.items():
        opposite = set()
        n = next(iter(N))
        shifted = {tuple(m[i] - n[i] for i in range(d)) for m in N if m != n}
        minus_n = tuple(-c for c in n)
        for w, Nw in sets.items():
            if len(N & Nw) == d - 1:
                opposite |= Nw - N
            if shifted <= Nw:
                rest = Nw - shifted
                if rest != {minus_n}:
                    e = next(iter(rest))
                    opposite.add(tuple(n[i] + e[i] for i in range(d)))
        result[vert] = opposite
    return result


def rayify(cone):
    r"""
    This function is a much better way of describing the secondary cone.
    Parameters:
        cone (Polyhedron): The cone to be rayified.
    Returns:
        list: A list of rays representing the cone.

    EXAMPLES::

        sage: from VoronoiCellToolBox import secondary_cone, rayify
        sage: sorted(rayify(secondary_cone([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]])))
        [(0, 0, 0, 0, 0, 1), (0, 0, 0, 1, -1, 1), (0, 0, 0, 1, 0, 0), (1, -1, 0, 1, 0, 0), (1, 0, -1, 0, 0, 1), (1, 0, 0, 0, 0, 0)]
    """
    return [ray[:] for ray in cone.rays()]


def pulling_triangulation(P):
    r"""
    Return a pulling triangulation of the polytope `P`: pull the first vertex `v` and cone it
    over a pulling triangulation of every facet not containing `v`, recursively. No new
    vertices are added, so for a Voronoi cell every simplex vertex is a cell vertex, whose
    coordinates are rational functions of the metric on the whole secondary cone.

    INPUT: ``P`` -- a polytope

    OUTPUT: a list of simplices, each a list of `\dim P + 1` vertices

    EXAMPLES::

        sage: from VoronoiCellToolBox import voronoi_cell, pulling_triangulation
        sage: T = pulling_triangulation(voronoi_cell([[2, -1], [-1, 2]])); len(T)
        4
        sage: T[0]
        [A vertex at (-2/3, -1/3), A vertex at (1/3, -1/3), A vertex at (2/3, 1/3)]
        sage: len(pulling_triangulation(voronoi_cell([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]])))
        34
    """
    #print(P.vertices())
    if(dimension(P) < 2):
        return [list(P.vertices())]
    # Pull the first vertex, so that the triangulation (and everything built on it) is reproducible
    v = P.vertices()[0]
    ans = []
    for F in P.facets():
        if not(v in F.vertices()):
            for triangle in pulling_triangulation(F.as_polyhedron()):
                ans += [[v] + triangle]
    return ans


#function to output Delaunay sets associated to each vertex, outputs as a dictionary with vertices as keys
def Delsets(VC, Q, **rangeorlist):
    r"""
    Computes a dictionary whereby each vertex in the Voronoi cell VC is associated with its Delaunay set.
    That is the collection of lattice points that are closest to it.

    Parameters:
        VC (VCell): The Voronoi cell.
        Q (Matrix): The quadratic form matrix.
        rangeorlist (list or range): A list of potential relevant vectors, or a range to search for them.

    Returns:
        dict: A dictionary mapping each vertex to its Delaunay set.

    EXAMPLES::

        sage: from VoronoiCellToolBox import voronoi_cell, Delsets
        sage: Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        sage: D = Delsets(voronoi_cell(Q), Q, range=2)
        sage: len(D), {P.n_vertices() for P in D.values()}     # 24 vertices, each a Delaunay simplex
        (24, {4})
    """
    Q = _rows(Q)
    d = len(Q)
    prv = None
    if "range" in rangeorlist:
        r = rangeorlist["range"]
        if (r <= 0):
            raise Exception("The range must be a positive integer.")
        prv = [list(vec) for vec in itertools.product(range(-r, r+1), repeat=d)]
        prv.remove([0] * d)
    elif "list" in rangeorlist:
        prv = rangeorlist["list"]
    else:
        raise Exception("A range or list of potential relevant vectors needs to be given.")
    

    vertices = VC.vertices_list()
    #loop through points p in pot_rel_vectors to check if distance from v to 0 is equal to distance to p
    # Ensure Q is a SageMath matrix
    Q = Matrix(QQ, Q)
    
    # Initialize the dictionary
    output_dict = {}
    
    # Loop over each vertex
    for v in vertices:
        v = vector(QQ, v)  # Ensure v is a vector
        distance_origin = (v * Q * v) # squared distance of vertex from the origin
        
        # Initialize an empty set for this vertex
        delset = [tuple([0]*d)]
        
        # Loop over each point
        for p in prv:
            p = vector(QQ, p)  # Ensure p is a vector
            distance_point = ((v-p) * Q * (v-p))  # squared distance of vertex to the point
            
            # Check if the distances match
            if distance_origin == distance_point:
                delset.append(tuple(p))  # Add the point to the set (convert to tuple for immutability)
                
        D = Polyhedron(vertices = delset)
        # Add the set to the dictionary
        output_dict[tuple(v)] = D
    
    return output_dict

#function to output Delaunay set associated to any point in the VC (do not use this if running Delsets anyway, just use the dictionary)
def Delset(Q, v, **rangeorlist):
    r"""
    For each vertex in the Voronoi cell VC computes the corresponding Delaunay set.

    Parameters:
        Q (Matrix): The quadratic form matrix.
        v (Vector): The vertex for which to compute the Delaunay set.
        rangeorlist (list or range): A list of potential relevant vectors, or a range to search for them.

    Returns:
        list: A list of vectors representing the Delaunay set.

    Usage:
        import random
        Q = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
        VC = VCell(Q, range=1)
        vertices = VC.vertices()
        v = vertices[random.randint(0, len(vertices))]
        Delset(Q, v, range = 2)#[(0, 0, 0), (-1, 0, 1), (-1, 1, 0), (-1, 1, 1), (0, 0, 1), (0, 1, 0)]
    """
    Q = _rows(Q)
    d = len(Q)
    prv = None
    if "range" in rangeorlist:
        r = rangeorlist["range"]
        if (r <= 0):
            raise Exception("The range must be a positive integer.")
        prv = [list(vec) for vec in itertools.product(range(-r, r+1), repeat=d)]
        prv.remove([0] * d)
    elif "list" in rangeorlist:
        prv = rangeorlist["list"]
    else:
        raise Exception("A range or list of potential relevant vectors needs to be given.")
    
    Q = Matrix(QQ, Q)
    v = vector(QQ, v)  # Ensure v is a vector
    distance_origin = (v * Q * v) # squared distance of vertex from the origin
    # Initialize a del set for this vertex
    delset = [tuple([0]*d)]
        
    # Loop over each point
    for p in prv:
        p = vector(QQ, p)  # Ensure p is a vector
        distance_point = ((v-p) * Q * (v-p))  # squared distance of vertex to the point
            
        # Check if the distances match
        if distance_origin == distance_point:
            delset.append(tuple(p))  # Add the point to the set (convert to tuple for immutability)
                
    return delset


def VertexFromRelevantVectors(B, Q, verbose = False):
    r"""
    Computes the vertex coordinates from a matrix of relevant vectors using the formula:
    $\vec{v} = \frac{1}{2} Q^{-1}(B^T)^{-1}b_Q$

    where b_Q is a vector whose i-th component is the Q-norm of the i-th column of B.

    Parameters:
        B (Matrix): A d x d matrix where columns are relevant vectors
        Q (Matrix): The positive definite metric matrix

    Returns:
        Vector: The computed vertex coordinates

    Usage:
        Q = [[2, -1], [-1, 2]]
        B = [[1, 0], [0, 1]]
        VertexFromRelevantVectors(B, Q)
    """
    d = len(Q)
    Q_matrix = Matrix(Q)
    B_matrix = Matrix(B)

    # Compute b_Q: Q-norms of each column of B
    b_Q = Matrix([[Qform(B_matrix.column(i), Q_matrix)] for i in range(d)])
    if verbose:
        print(f"b_Q = \n{b_Q}")

    # Compute (B^T)^{-1} * b_Q
    B_T_inv = B_matrix.transpose().inverse()
    intermediate = B_T_inv * b_Q
    if verbose:
        print(f" B transpose inverse = \n{B_T_inv}")
        print(f"intermediate = bt inv * bq \n{intermediate}")

    # Compute Q^{-1} * intermediate
    Q_inv = Q_matrix.inverse()
    result = Q_inv * intermediate
    if verbose:
        print(f" Q inverse = {Q_inv}")
        print(f"result before scaling = {result}")
        print(f"result after scaling = {[result[i,0]/2 for i in range(d)]}")

    return vector( [result[i,0]/2 for i in range(d)] )


def smMatrixSimplex(Q, verbose = False, **rangeorlist):
    r"""
    Computes the second moment of the Voronoi cell defined by the quadratic form Q,
    in the inherited metric given by Q.

    The computation uses the pulling triangulation and the formula:
    SM_Q(T) = (det(T) / (sqrt(det Q) * (d+2)!)) * (||(d+1)*s_hat||_Q^2 + sum_i ||s_i||_Q^2)

    where T is the matrix of edge vectors, s_hat is the barycenter, and s_i are the vertices.

    Parameters:
        Q (Matrix): The positive definite metric matrix.
        rangeorlist (list or range): A list of potential relevant vectors, or a range to search for them.

    Returns:
        float: The second moment of the Voronoi cell.

    Usage:
        Q = [[2, -1], [-1, 2]]
        smMatrixSimplex(Q, range=2)
    """
    Q = _rows(Q)
    VC = None
    if "range" in rangeorlist:
        r = rangeorlist["range"]
        VC = VCell(Q, range=r)
    elif "list" in rangeorlist:
        prv = rangeorlist["list"]
        VC = VCell(Q, list=prv)
    else:
        raise Exception("A range or list of potential relevant vectors needs to be given.")

    pt = pulling_triangulation(VC)
    if verbose:
        for index, triangle in enumerate(pt):
            print(f"triangles[{index}] = {triangle}")
    total = 0
    d = len(Q)
    detQ = numpy.linalg.det(numpy.array(Q))

    dp2fact = math.factorial(d+2)

    for triangle in pt:
        # Compute the second moment of the simplex defined by the vertices in triangle
        # Compute the barycenter
        s_hat = sum([vector(v) for v in triangle])/len(triangle)

        # Compute the matrix T of edge vectors (determinant)
        T = [list(vector(triangle[i]) - vector(triangle[0])) for i in range(1, len(triangle))]
        detT = abs(Matrix(QQ, T).det())

        # Compute the second moment contribution: det(T) / (sqrt(det Q) * (d+2)!) * (||(d+1)*s||_Q^2 + sum ||v_i||_Q^2)
        partial_total = Qform((d+1)*s_hat, Q) + sum([Qform(v, Q) for v in triangle])
        total += detT * partial_total
        if verbose:
            print(f"partial_total = {partial_total}")

    return total / (math.sqrt(detQ) * dp2fact)


def second_moment(Q, **rangeorlist):
    r"""
    Return the exact second moment `SM(Q) = \int_{V_Q} x^T Q x \, dx` of the Voronoi cell of `Q`,
    with Euclidean measure (the cell has volume 1).

    The cell is triangulated (:func:`pulling_triangulation`) and each simplex with vertices
    `s_0, \dots, s_d` contributes `|\det T| (\|\sum_i s_i\|_Q^2 + \sum_i \|s_i\|_Q^2) / (d+2)!`,
    where `T` has rows `s_i - s_0`.

    INPUT:

    - ``Q`` -- a positive definite symmetric matrix
    - ``range=r`` or ``list=[...]`` -- (optional) candidate relevant vectors, as in :func:`VCell`;
      by default they are chosen automatically (:func:`voronoi_cell`)

    OUTPUT: a rational number

    EXAMPLES::

        sage: from VoronoiCellToolBox import second_moment
        sage: second_moment([[2, -1], [-1, 2]])
        5/18
        sage: second_moment([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]])
        19/32
        sage: second_moment([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        1/4
    """
    Q = _rows(Q)
    VC = VCell(Q, **rangeorlist) if rangeorlist else voronoi_cell(Q)
    pt = pulling_triangulation(VC)
    total = 0
    d = len(Q)

    import math
    dp2fact = math.factorial(d+2)
    for triangle in pt:
        r"""
        det T / (sqrt(det Q) (d+2)! ) ( || (d+1) s ||_Q^2 + sum ||v_i||_Q^2 )
        """
        # Compute the second moment of the simplex defined by the vertices in triangle
        #  Compute the baricentre and the determinant of the matrix T
        s = sum([vector(v) for v in triangle])/len(triangle)
        T = [list(vector(triangle[i]) - vector(triangle[0]))  for i in range(1, len(triangle))]
        detT = abs(Matrix(QQ, T).det())
        partial_total = Qform((d+1)*s, Q) + sum([Qform(v, Q) for v in triangle])
        total += detT * partial_total
    return total / ( dp2fact)


def normalized_second_moment(Q, **rangeorlist):
    r"""
    Return the dimensionless (normalized) second moment `G(Q) = SM(Q) / (d \det(Q)^{1/d})`,
    exactly, as a symbolic expression. It is invariant under scaling `Q` and under changes of
    basis `Q \mapsto U^T Q U`, `U \in GL_d(\ZZ)`; it is the quantity minimized in the lattice
    quantizer problem.

    INPUT: as for :func:`second_moment`

    EXAMPLES::

        sage: from VoronoiCellToolBox import normalized_second_moment
        sage: G = normalized_second_moment([[2, -1], [-1, 2]]); G
        5/108*sqrt(3)
        sage: G.n(digits=8)
        0.080187537
        sage: normalized_second_moment([[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]).n(digits=6)  # A_4^*
        0.0775588
    """
    Q = _rows(Q)
    d = len(Q)
    return second_moment(Q, **rangeorlist) / (d * SR(Matrix(QQ, Q).det())**(QQ(1)/d)) 