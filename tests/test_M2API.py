#/tests/test_M2API.py
#This files tests all M2 functions called from the orchestration module

#Functions to test
### barycentre
### favouriteMatrix
### fromRelevantVectorsToVertex
### inverseCofactorMatrix
### isSame
### Listify
### listOf1Minors
### makepos
### qnorm
### qNormMatrixFormat
### secondMoment
### SmPoly
### symMatricesRing
### VectorizedVertex

import random
import unittest
from sage.all import Matrix, QQ
from VoronoiCellToolBox.orchestration import barycentre, VectorizedVertex_m2


class TestM2API(unittest.TestCase):

    def test_barycentre(self):
        V = [[1, 2, 3], [4 ,5, 6]]
        result = barycentre(V)
        
        # Check if the result is a non-empty string (basic sanity check)
        self.assertIsInstance(result, list)
        self.assertTrue(len(result) == 2)
        self.assertTrue(result == [[2], [5]])
    
    def test_barycentre_randomized_inputs(self):
        for _ in range(10):  # Run the test 10 times with different random inputs
            dim = random.randint(1, 5)  # Random dimension between 1 and 5
            num_points = random.randint(2, 10)  # Random number of points between 2 and 10
            
            # Generate random points
            V = [[random.randint(-10, 10) for _ in range(num_points)] for _ in range(dim)]
            
            # Calculate expected barycentre manually
            expected = [[sum(V[j][i] for i in range(num_points)) / num_points] for j in range(dim)]
            
            # Get result from barycentre function
            result = barycentre(V)
            
            #print("Tested barycentre with V =", V)
            #print("Expected =", expected)
            #print("Got =", result)
            # Assert that the result matches the expected value
            for j in range(dim):
                self.assertAlmostEqual(result[j][0], expected[j][0], delta=1e-6)


class TestVectorizedVertex(unittest.TestCase):

    def test_single_matrix(self):
        Q = [[1, 0], [0, 1]]
        listMatrices = [[[1, 0], [0, 1]]]
        result = VectorizedVertex_m2(listMatrices, Q)
        # One input matrix -> one vertex column: "matrix {{1/2}, {1/2}}"
        self.assertIsInstance(result, str)
        self.assertIn("matrix", result)
        self.assertIn("1/2", result)

    def test_two_matrices(self):
        Q = [[1, 0], [0, 1]]
        listMatrices = [[[1, 0], [0, 1]], [[0, 1], [1, 0]]]
        result = VectorizedVertex_m2(listMatrices, Q)
        # Two input matrices -> two vertex columns horizontally concatenated
        self.assertIsInstance(result, str)
        self.assertIn("matrix", result)

    def test_empty_list(self):
        Q = [[1, 0], [0, 1]]
        result = VectorizedVertex_m2([], Q)
        # Empty input -> d x 0 empty matrix
        self.assertIsInstance(result, str)
        self.assertIn("matrix", result)
