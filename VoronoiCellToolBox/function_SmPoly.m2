-- Takes a dimension d and matVertices structure
-- matVertices encodes a triangulation of the voronoi cell
-- matvertices is a matrix, where each row corresponds to a triangle
-- in each row an entry corresponds to a fromRelevantVectorsToVertex of the triangle
-- a fromRelevantVectorsToVertex is encoded by a matrix
-- output is the rational function giving the second moment for any Voronoi cell in
-- the same chamber as our favorite matrix in terms of the entries of the matrix Q
-- SmPoly(d, matVertices, A, verbose)
-- Input:
--   d: dimension
--   matVertices: list of matrices encoding triangulation
--   A: metric symmetric matrix
--   verbose: boolean flag for verbose output
-- Output:
--   Rational function for the second moment polynomial of the Voronoi cell
-- Description:
--   Computes the second moment polynomial for a Voronoi cell using its triangulation.
-- Example:
--   SmPoly(2, { {matrix{{1,0},{0,1}}, matrix{{0,1},{1,0}}} }, matrix{{2,0},{0,3}}, false)
SmPoly = (d, matVertices, A, verboseComputationProgress, verboseSecondMoments, verboseVertex) -> (
    if(verboseComputationProgress) then print("Starting SmPoly computation for dimension ", toString(d)) else null;
    if(verboseComputationProgress) then print("Number of triangles: ", toString(#matVertices)) else null;
    if(verboseComputationProgress) then print("Verbose mode enabled") else null;
    if(verboseComputationProgress) then print("Vertex matrices", toString matVertices) else null;
    if(verboseComputationProgress) then print("Metric matrix", toString A) else null;
    G = d*(d+1)//2;
    R = QQ[q_0..q_(G-1)];
    Q = genericSymmetricMatrix(R, q_0, d);
    Zpoly = 0;
    lvalues = Listify(A, d);
    l = # matVertices;
    for i from 0 to l-1 do(
        if(verboseComputationProgress) then print(concatenate(toString i, " of total ", toString l, " - new loop")) else null;
        concatVertices = VectorizedVertex(matVertices_i, Q, d, false);
        if(verboseVertex) then print(concatenate("concatVertices = \n", toString concatVertices)) else null;
        if(verboseComputationProgress) then print(concatenate(toString i, " of total ", toString l, " - loop compute secondMoment")) else null;
        ply = secondMoment(concatVertices, d, Q);
        if(verboseComputationProgress) then print(concatenate("l values = ", toString lvalues)) else null;
        newPolyTriangle = makepos(ply, lvalues, d, R);
        if(verboseSecondMoments) then print(concatenate("newPolyTriangle = \n", toString newPolyTriangle)) else null;
        if(verboseComputationProgress) then print(concatenate(toString i, " of total ", toString l, " - loop add ")) else null;
        Zpoly = Zpoly + newPolyTriangle;
    );
    return Zpoly
);
