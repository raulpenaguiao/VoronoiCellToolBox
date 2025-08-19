r"""
Computes the Voronoi cell of a lattice given by a matrix
===============================================

TODO description

EXAMPLES::
TODO examples
"""

# ****************************************************************************
#                           Copyright (C)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numpy
import itertools
import random
from sage.all import Polyhedron, vector, dimension, Matrix, Rational
from sage.arith.functions import LCM_list

def VCell(Q, **rangeorlist):
    """
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
    - The determinant computation uses NumPy's `linalg.det` method, which may
      introduce floating-point inaccuracies for very large or very small values.
    """
    if numpy.linalg.det(numpy.array(Q)) <=0: 
        raise Exception("your matrix is not positive-definite.")
        return "your matrix is not positive-definite."
    if (numpy.array(Q).transpose() != numpy.array(Q)).any():
        raise Exception("your matrix is not symmetric.")
        return "your matrix is not symmetric."
    d = len(Q)
    if "range" in rangeorlist:
        r = rangeorlist["range"]
        if (r <= 0):
            raise Exception("The range must be a positive integer.")
        p = [list(vec) for vec in itertools.product(range(-r, r+1), repeat=d)]
        p.remove([0] * d)
    elif "list" in rangeorlist:
        p = rangeorlist["list"]
    else:
        raise Exception("A range or list of potential relevant vectors needs to be given.")
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


def pos(i: int, j: int, d: int) -> int:
    """
    Maps the position (i, j) on the $d \times d$ matrix to a list index, so that 
    each symmetric entry can be encoded in a list of size $\frac{d\times (d+1)}{2}$.

    Parameters:
        i (int): The row index of the matrix.
        j (int): The column index of the matrix.
        d (int): The dimension of the matrix.

    Usage:
        pos(0, 1, 3)  # Returns 1
        pos(1, 2, 3)  # Returns 4
        pos(2, 2, 3)  # Returns 5

    Behavior:
        The function should return a unique index for each pair (i, j) such that
        all pairs (i, j) with i < j are assigned indices before any pairs with
        i = j.

    Note:
        Also works when i > j. Undefined when $j < 0$, $i < 0$, $i \geq d$ or $j \geq d$.
    """
    if i>j: return pos(j, i, d)
    return i*d+j-(((i+1)*i)//2)



def matrixify(list_of_ute, d):
    """
    Converts a list of upper triangular entries into a $d \times d$ matrix.

    Parameters:
        list_of_ute (list): The list of entries of the symmetric matrix.
        d (int): The dimension of the matrix.

    Returns:
        list: A $d \times d$ matrix represented as a list of lists.
    """
    return [[list_of_ute[int(pos(i, j, d))]  for j in range(d)] for i in range(d)]

#calculates equation of wall in secondary cone
def eq_of_wall(rvecs, nrvec, d, verbose = False):
    if verbose:
        print(rvecs, nrvec, d)
    Tinv = Matrix(rvecs).inverse()
    ineq = [0]*(d*(d+1)//2)
    for i in range(0, d):
        for j in range(0, d):
            ineq[pos(i, j, d)] += nrvec[i]*nrvec[j]
            if(verbose):
                print(nrvec[i]*nrvec[j], " - > (", i, ", ", j, ")")
    if(verbose): 
        print(Tinv)
    for i in range(0, d):
        for j in range(0, d):
            for k in range(0, d):
                for l in range(0, d):
                    ineq[pos(k, l, d)] -= nrvec[i]*Tinv[i][j]*rvecs[j][k]*rvecs[j][l]
                    if ( verbose): 
                        print( nrvec[i]*Tinv[i][j]*rvecs[j][k]*rvecs[j][l], " -> (", k,", ", l,")")
                        print(" (", i, " , ", j, ", ", k ,", ", l, ") - ", nrvec[i], ", ", Tinv[i][j], ", ", rvecs[j][k], ", ", rvecs[j][l], " -> (", k, ", ", l, ") ")
    return ineq #do not include b becuase it's homogeneous and b is always 0
#\sum a_{i, j} * q_{i, j} <= 0 for newv to NOT be relevant   

def reduce(v):
    #lcm =  LCM_list([vi.denominator() for vi in v])
    #return vector(gcdReduce([vi*lcm for vi in v]))
    try:
        lcm =  LCM_list([vi.denominator() for vi in v])
    except:
        return reduce([Rational(int(round(vi*10**9)),int(10**9)) for vi in v])
        #return vector(gcdReduce([vi for vi in v]))
    return vector(gcdReduce([vi*lcm for vi in v]))

def gcdReduce(v):
    gcd = GCD_list(v)
    return [vi/gcd for vi in v]


def gcd(a, b):
    if (a < 0): return gcd(-a, b)
    if (a > b): return gcd(b, a)
    if (a < 1 and a > -1): return b
    return gcd(b%a, a)

def GCD_list(v):
    if (len(v) == 0): return 1
    if (len(v) == 1): return v[0]
    if (len(v) == 2): return gcd(v[0], v[1])
    return GCD_list([gcd(v[0], v[1])] + v[2:] )

def relevant_vector(Q, F):
    v = (1/2)*vector(F.ambient_Hrepresentation()[0][1:])*(Matrix(Q).inverse())
    return -reduce(v) #We have no idea why - seems to work every time but
    #it is either 1 or -1
    #1/2 * t * Q-1
    #This gives a scalar multiple of the relevant vector, but which one?

def Qnorm(v, Q):
    return (Matrix([[vi] for vi in v]).transpose()*Matrix(Q)*Matrix([v]).transpose())[0][0]


def eps(vlist, Q):
    #Given a list ov non-zero lattice points, eps = the equidistant tuple to all points and zero
    T = Matrix([list(v) for v in vlist])
    return (2*T*Matrix(Q)).inverse()*vector([Qnorm(v, Q) for v in vlist])
    
def inc_fac(P, Q):
    vdict = {}
    for facet in P.facets():
        r_vec = relevant_vector(Q, facet)
        for v in facet.vertices():
            if v in vdict:
                vdict[v] += [r_vec]
            else:
                vdict[v] = [r_vec]
    return vdict


def secondary_cone(Q, prv, verbose = False):
    """
    Computes the polyhedral cone that describes all PSD matrices with the same voronoi triangulation as Q
    
    Parameters
    ----------
    Q : Square matrix
    prv : potential relevant vectors

    Returns
    -------

    See Also
    --------

    Examples
    --------
    """
    d = len(Q)
    VC = VCell(prv, Q)
    ifs = inc_fac(VC, Q)
    nrvs = []
    rvs = [] #nrvs = prv - relevant vectors 
    for v in prv:
        isRelevant = False
        for F in VC.facets():
            if vector(v) == vector(relevant_vector(Q, F)):
                isRelevant = True
        if isRelevant == True:
            rvs += [v]
        else:
            nrvs += [v]
    ineqs = []
    if verbose:
        print("v4 - len of ifs = ", len(ifs) , " and len of n relevant vecs = ", len(nrvs))
    for vert in ifs:
        for vec in prv:
            if not(vec in ifs[vert]):
                #compute wall crossing between ifs[vert] and vec
                eq_wall = eq_of_wall(ifs[vert], vec, d, verbose)
                ineqs += [ [0] + eq_wall]
                Qxeq_wall = 0
                for i in range(d):
                    for j in range(i+1):
                        Qxeq_wall += Q[i][j]*eq_wall[pos(i, j, d)]
                if verbose and Qxeq_wall < 0:
                    print(vec, " - ", vert, " - ", ifs[vert], " - ", eq_of_wall(ifs[vert], vec, d))
    return Polyhedron(ieqs = ineqs)

def rayify(cone):
    return [ray[:] for ray in cone.Vrepresentation()[1:]]


def pulling_triangulation(P):
    #print(P.vertices())
    if(dimension(P) < 2):
        return [list(P.vertices())]
    vert = P.vertices()
    index = random.randint(0, len(vert)-1)
    v = vert[index]
    ans = []
    for F in P.facets():
        if not(v in F.vertices()):
            for triangle in pulling_triangulation(F.as_polyhedron()):
                ans += [[v] + triangle]
    return ans