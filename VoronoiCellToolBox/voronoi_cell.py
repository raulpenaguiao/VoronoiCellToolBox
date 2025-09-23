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

import math
import numpy
import itertools
import random
from sage.all import Polyhedron, vector, dimension, Matrix, Rational, QQ
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


def FromMatrixPositionToListPosition(i: int, j: int, d: int) -> int:
    """
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
    """
    Converts a list of upper triangular entries into a $d \times d$ matrix.

    Parameters:
        list_of_ute (list): The list of entries of the symmetric matrix.
        d (int): The dimension of the matrix.
    
    Returns:
        list: A $d \times d$ matrix represented as a list of lists.
    
    Usage:
        matrixify([3, 1, 5, 3, -1, 2], 3) #Returns [[3, 1, 5], [1, 3, -1], [5, -1, 2]]
    """
    return [[list_of_ute[int(FromMatrixPositionToListPosition(i, j, d))]  for j in range(d)] for i in range(d)]

#calculates equation of wall in secondary cone
def eq_of_wall(neighborVectors: list[vector], isolatedVector: vector, d: int, verbose: bool = False) -> list:
    """
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

    Usage:
        rvecs = [(0, 1, 1), (1, 1, 1), (0, 0, 1)]
        nrvec = [-1, 0, 0]
        eq_of_wall(rvecs, nrvec, 3, True) #Returns [2, 2, 2, 0, 0, 0]

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
    """
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
    """
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
    """
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
    """
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
    """
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
    v = (1/2)*vector(F.ambient_Hrepresentation()[0][1:])*(Matrix(Q).inverse())
    return -lcmReduce(v) 
    #We have no idea why - seems to work every time but
    #it is either 1 or -1
    #1/2 * t * Q-1
    #This gives a scalar multiple of the relevant vector, but which one?

def Qform(v, Q):
    """
    For a matrix $Q$ defining a metric on $Z^d$, and a vector v, this function 
    computes the quadratic form $v^T Q v$.
    Parameters:
        v (Vector): The vector to be evaluated.
        Q (Matrix): The metric matrix.
    Returns:
        Scalar: The value of the quadratic form.
    Usage:
        Qform(vector([1, 2]), Matrix([[2, 0], [0, 2]])) #Returns 8
    """
    return (Matrix([[vi] for vi in v]).transpose()*Matrix(Q)*Matrix([v]).transpose())[0][0]


def equidistantPoint(vlist, Q):
    """
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
    """
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
    """
    Computes a Polyhedron type object that describes all matrices Q' 
    that have the same relevant vector structure as Q 
    
    Parameters:
        Q (Matrix): The metric matrix defining the Voronoi cell.
        range = r (int): The range of potential relevant vectors.
        OR
        list = prv (list): A list of potential relevant vectors.

    Returns:
        Polyhedron: The secondary cone of Q.

    Usage:
        Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        prv = [list(vec) for vec in itertools.product(range(-2, 3), repeat=3)]
        prv.remove([0] * 3)
        secondary_cone(Q, prv)
    """

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
    
    VC = VCell(Q, list=prv)
    ifs = relevantVectorDictionary(VC, Q)
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
                        Qxeq_wall += Q[i][j]*eq_wall[FromMatrixPositionToListPosition(i, j, d)]
                if verbose and Qxeq_wall < 0:
                    print(vec, " - ", vert, " - ", ifs[vert], " - ", eq_of_wall(ifs[vert], vec, d))
    return Polyhedron(ieqs = ineqs)

def rayify(cone):
    """
    This function is a much better way of describing the secondary cone.
    Parameters:
        cone (Polyhedron): The cone to be rayified.
    Returns:
        list: A list of rays representing the cone.

    Usage:
        Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        prv = [list(vec) for vec in itertools.product(range(-2, 3), repeat=3)]
        prv.remove([0] * 3)
        SC = secondary_cone(Q, prv)
        rayify(SC)
    """
    return [ray[:] for ray in cone.Vrepresentation()[1:]]


def pulling_triangulation(P):
    """
    Computes a pulling triangulation of a polyhedron P.
    Parameters:
        P (Polyhedron): The polyhedron to be triangulated.

    Returns:
        list: A list of triangles representing the triangulation.

    Usage:
        Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        P = VCell(Q, range=2)
        len(pulling_triangulation(P)) #Returns 34
    """
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


#function to output Delaunay sets associated to each vertex, outputs as a dictionary with vertices as keys
def Delsets(VC, Q, **rangeorlist):
    """
    Computes a dictionary whereby each vertex in the Voronoi cell VC is associated with its Delaunay set.
    That is the collection of lattice points that are closest to it.

    Parameters:
        VC (VCell): The Voronoi cell.
        Q (Matrix): The quadratic form matrix.
        rangeorlist (list or range): A list of potential relevant vectors, or a range to search for them.

    Returns:
        dict: A dictionary mapping each vertex to its Delaunay set.

    Usage:
        Q = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
        VC = VCell(Q, range=2)
        len(Delsets(VC, Q, range=2))
    """
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
    """
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


def second_moment(Q, **rangeorlist):
    """
    Computes the second moment of the Voronoi cell defined by the quadratic form Q.

    Parameters:
        Q (Matrix): The quadratic form matrix.
        rangeorlist (list or range): A list of potential relevant vectors, or a range to search for them.

    Returns:
        float: The second moment of the Voronoi cell.

    Usage:
        Q = [[2, -1], [-1, 2]]
        second_moment(Q, range=2) #Returns 
    """
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
    total = 0
    d = len(Q)
    detQ = numpy.linalg.det(numpy.array(Q))

    import math
    dp2fact = math.factorial(d+2)
    print("d = ", d, " detQ = ", detQ, " (d+2)! = ", dp2fact)
    for triangle in pt:
        """
        det T / (sqrt(det Q) (d+2)! ) ( || (d+1) s ||_Q^2 + sum ||v_i||_Q^2 )
        """
        print("triangle = ", triangle)
        # Compute the second moment of the simplex defined by the vertices in triangle
        #  Compute the baricentre and the determinant of the matrix T
        s = sum([vector(v) for v in triangle])/len(triangle)
        T = [list(vector(triangle[i]) - vector(triangle[0]))  for i in range(1, len(triangle))]
        print("T = ", T, " s = ", s)
        detT = abs(numpy.linalg.det(numpy.array(T)))
        print("detT = ", detT)
        print(" ||v_i||_Q^2 = ", [Qform(v, Q) for v in triangle])
        print(" (d+1)*s = ", (d+1)*s)
        partial_total = Qform((d+1)*s, Q) + sum([Qform(v, Q) for v in triangle])
        print("partial_total = ",   partial_total)
        total += detT * partial_total
        print("total = ", total)
    print("final total = ", total)
    return total / ( dp2fact) 