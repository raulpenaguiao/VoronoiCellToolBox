from VoronoiCellToolBox.voronoi_cell import VCell, pulling_triangulation, relevantVectorDictionary, _rows
from itertools import chain, combinations
from sage.all import Matrix, QQ

def macaulifyMatrix(lot):
    return "matrix" + str(lot).replace(')', '}').replace('(', '{').replace(']', '}').replace('[', '{')

def AddBreaks(string, nChar):
    ans = ""
    counter = 0
    for c in string:
        ans += c
        counter += 1
        if counter == nChar:
            if c in ["-", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "m", "a", "t", "r", "i", "x"]:
                counter -= 1
            else:
                ans += "\n "
                counter = 0
    return ans

def Transpose(mat):
    return [[mat[i][j] for i in range(len(mat))] for j in range(len(mat[0]))]

def selectDIndependentVectors(vecs, d):
    """
    Select exactly d linearly independent vectors from vecs.

    On the boundary of a secondary cone, a Voronoi vertex can be incident to more
    than d facets, giving more than d relevant vectors. fromRelevantVectorsToVertex
    requires a square d x d matrix B, so we must take exactly d independent ones.
    Any such subset yields the same vertex coordinates because the vertex is
    equidistant to all its relevant lattice points.
    """
    selected = []
    for v in vecs:
        candidate = selected + [list(v)]
        if Matrix(QQ, candidate).rank() == len(candidate):
            selected.append(list(v))
        if len(selected) == d:
            break
    return selected

def FormatPullingTrigMatrix(Q):
    Q = _rows(Q)
    d = len(Q)
    VC = VCell(Q, range = 2)
    irv = relevantVectorDictionary(VC, Q)
    pt = pulling_triangulation(VC)
    string = "{"
    for triangle in pt:
        string += "{"
        for v in triangle:
            rvecs = selectDIndependentVectors(irv[v], d)
            string += macaulifyMatrix(Transpose(rvecs))
            string += ","
        string = string[:-1]
        string += "},"
    string = string[:-1]
    string += "}"
    return AddBreaks(string, 100)

