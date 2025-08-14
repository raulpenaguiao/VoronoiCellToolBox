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
