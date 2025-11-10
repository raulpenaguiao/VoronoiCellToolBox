-- Takes a dimension d and matVertices structure
-- matVertices encodes a triangulation of the voronoi cell
-- matvertices is a matrix, where each row corresponds to a triangle
-- in each row an entry corresponds to a fromRelevantVectorsToVertex of the triangle
-- a fromRelevantVectorsToVertex is encoded by a matrix
-- output is the rational function giving the second moment for any Voronoi cell in 
-- the same chamber as our favorite matrix in terms of the entries of the matrix Q
-- SmPoly(d, matVertices)
-- Input:
--   d: dimension
--   matVertices: list of matrices encoding triangulation
--   A: metric symmetric matrix
-- Output:
--   Rational function for the second moment polynomial of the Voronoi cell
-- Description:
--   Computes the second moment polynomial for a Voronoi cell using its triangulation.
-- Example:
--   SmPoly(2, { {matrix{{1,0},{0,1}}, matrix{{0,1},{1,0}}} })
SmPoly = (d, matVertices, A, verbose) -> (
    G = d*(d+1)//2;
    R = QQ[q_0..q_(G-1)];
    Q = genericSymmetricMatrix(R, q_0, d);
    Zpoly = 0;
    lvalues = Listify(A);
    l = # matVertices;
    for i from 0 to l-1 do(
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - new loop")) else null;
        concatVertices = VectorizedVertex(matVertices_i, Q, d);
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - loop compute secondMoment")) else null;
        ply = secondMoment(concatVertices, d, Q);
        if(verbose) then print(concatenate("l values = ", toString lvalues)) else null;
        newPolyTriangle = makepos(ply, lvalues, d, R);
        if(verbose) then print(concatenate("newPolyTriangle = ", toString newPolyTriangle)) else null;
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - loop add ")) else null;
        Zpoly = Zpoly + newPolyTriangle;
    );
    return Zpoly
);
