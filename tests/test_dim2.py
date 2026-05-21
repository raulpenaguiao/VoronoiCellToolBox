# tests/test_dim2.py
# This file tests all voronoi cell functions on dimension 2
# We will be using the matrices Q = [[2, -1], [-1, 2]],  [[1, 0], [0, 1]],  [[2, 1], [1, 2]]
# We will compute the secondary cone, the voronoi cell, its second moment and the sm polynomial

import random
import unittest
from sage.all import Matrix, QQ, vector, Rational, Polyhedron
from VoronoiCellToolBox.voronoi_cell import (
    second_moment,
    VCell
)

Q21 = [[2, -1], [-1, 2]]
Q22 = [[1, 0], [0, 1]]
Q23 = [[2, 1], [1, 2]]

Q31 = [[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]
Q32 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Q33 = [[9, 2, 5], [2, 4, 4], [5, 4, 8]]

class TestSecondMoment(unittest.TestCase):
    """Test suite for second moment function"""
    def test_secondmoment_is_expected(self):
        result21 = second_moment(Q21, range = 2)
        self.assertTrue( result21 == 1/12 )


class TestVoronoiCells(unittest.TestCase):
    """Test suite for VCell function"""

    def test_vcell_Q21_is_expected(self):
        result21 = VCell(Q21, range = 2)
        verts21 =[ 
            [-1/3, 1/3],
            [1/3, -1/3],
            [2/3, 1/3],
            [1/3, 2/3],
            [-2/3, -1/3],
            [-1/3, -2/3]
        ]
        expected21 = Polyhedron( vertices = verts21 )
        self.assertTrue( result21 == expected21 )
    
    def test_vcell_Q22_is_expected(self):
        result22 = VCell(Q22, range = 2)
        verts22 = [ 
            [1/2, 1/2],
            [1/2, -1/2],
            [-1/2, -1/2],
            [-1/2, 1/2]
        ]
        expected22 = Polyhedron( vertices = verts22 )
        self.assertTrue( result22 == expected22 )
    
    def test_vcell_Q23_is_expected(self):
        result23 = VCell(Q23, range = 2)
        verts23 = [ 
            [1/3, 1/3],
            [-1/3, -1/3],
            [-2/3, 1/3],
            [-1/3, 2/3],
            [2/3, -1/3],
            [1/3, -2/3]
        ]
        expected23 = Polyhedron( vertices = verts23 )
        self.assertTrue( result23 == expected23 )
    
    def test_vcell_Q31_is_expected(self):
        result31 = VCell(Q31, range = 2)
        verts31 =[ 
            [1/2, 1/2, 1/2],
            [1/2, 1/2, -1/2],
            [1/2, -1/2, 1/2],
            [1/2, -1/2, -1/2],
            [-1/2, 1/2, 1/2],
            [-1/2, 1/2, -1/2],
            [-1/2, -1/2, 1/2],
            [-1/2, -1/2, -1/2],
        ]
        expected31 = Polyhedron( vertices = verts31 )
        self.assertTrue( result31 == expected31 )
    
    def test_vcell_Q32_is_expected(self):
        result32 = VCell(Q32, range = 2)
        verts32 =[ 
            [1/2, 1/2, 1/2],
            [1/2, 1/2, -1/2],
            [1/2, -1/2, 1/2],
            [1/2, -1/2, -1/2],
            [-1/2, 1/2, 1/2],
            [-1/2, 1/2, -1/2],
            [-1/2, -1/2, 1/2],
            [-1/2, -1/2, -1/2],
        ]
        expected32 = Polyhedron( vertices = verts32 )
        self.assertTrue( result32 == expected32 )

