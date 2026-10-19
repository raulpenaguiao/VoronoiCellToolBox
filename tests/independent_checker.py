# Independent checker for SM(Q) = int_{V_Q} x^T Q x dx (Euclidean measure), used as ground truth.
# It shares no code with the package: bisectors over {-r..r}^d, Sage's own triangulation, exact arithmetic.
import itertools
from sage.all import QQ, matrix, vector, Polyhedron, factorial


def sm_independent(Q, r=3):
    """Exact SM(Q). Asserts that the cell has volume 1 (no relevant vector missed)."""
    Q = matrix(QQ, Q)
    d = Q.nrows()
    ieqs = []
    for v in itertools.product(range(-r, r+1), repeat=d):
        v = vector(QQ, v)
        if v == 0:
            continue
        ieqs.append([v*Q*v/2] + list(-(Q*v)))    # x^T Q v <= v^T Q v / 2
    P = Polyhedron(ieqs=ieqs, base_ring=QQ)
    assert P.volume() == 1, "candidate range too small"
    total = 0
    for T in P.triangulate():
        S = [vector(QQ, P.Vrepresentation(i)) for i in T]
        vol = abs(matrix([s - S[0] for s in S[1:]]).det())/factorial(d)
        c = sum(S)
        total += vol/((d+1)*(d+2)) * (sum(s*Q*s for s in S) + c*Q*c)
    return total
