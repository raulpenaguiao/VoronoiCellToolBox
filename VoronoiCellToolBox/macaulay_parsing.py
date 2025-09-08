from VoronoiCellToolBox.voronoi_cell import VCell, pulling_triangulation, relevantVectorDictionary
from itertools import chain, combinations

def macaulifyMatrix(lot):
    return "matrix" + str(lot).replace(')', '}').replace('(', '{').replace(']', '}').replace('[', '{')

def PrintWithBreaks(string, nChar):
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

def FormatPullingTrigMatrix(Q):
    VC = VCell(Q, range = 2)
    irv = relevantVectorDictionary(VC, Q)
    pt = pulling_triangulation(VC)
    string = "{"
    for triangle in pt:
        string += "{"
        for v in triangle:
            string += macaulifyMatrix(Transpose(irv[v]))
            string += ","
        string = string[:-1]
        string += "},"
    string = string[:-1]
    string += "}"
    return PrintWithBreaks(string, 100)