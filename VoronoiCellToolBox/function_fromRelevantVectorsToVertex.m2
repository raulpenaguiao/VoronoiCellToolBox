-- fromRelevantVectorsToVertex(B, Q, d)
-- Input:
--   B: matrix of relevant vectors, as column vectors
--   Q: symmetric matrix
--   d: dimension
-- Output:
--   Vertex coordinates as a column vector and as a polynomial in entries of Q
-- Description:
--   Computes the vertex of a Voronoi cell from relevant vectors and Q.
-- Example:
--   toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{2,0},{0,3}}, 2)
-- -- Output: matrix {{1/2}, {1/2}}
--   toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{q11, q12},{q12,q22}}, 2)
-- -- Output: matrix {{(-q11*q22+q12*q22)/(2*q12^2-2*q11*q22)}, {(q11*q12-q11*q22)/(2*q12^2-2*q11*q22)}}
fromRelevantVectorsToVertex = (B, Q, d, verbose)->(
    if(verbose) then print("Computing vertex from relevant vectors with B = ", toString B, " and Q = ", toString Q) else null;
    if(verbose) then print("Dimension d = ", toString d) else null;
    if(verbose) then print("Computing Q-norm matrix format of B") else null;
    if(verbose) then print("Q-norm matrix format: ", toString qNormMatrixFormat(B, Q, d, verbose)) else null;
    if(verbose) then print("Computing inverse of transpose of B") else null;
    if(verbose) then print("Inverse of transpose of B: ", toString inverse(transpose(B))) else null;
    if(verbose) then print("Computing inverse of cofactor matrix of Q") else null;
    if(verbose) then print("Inverse of cofactor matrix of Q: ", toString inverseCofactorMatrix(Q, d)) else null;
    if(verbose) then print("Combination 1 ", toString(inverse(transpose(B))*qNormMatrixFormat(B, Q, d, verbose))) else null;
    if(verbose) then print("Combination 2 ", toString(inverseCofactorMatrix(Q, d)*inverse(transpose(B))*qNormMatrixFormat(B, Q, d, verbose))) else null;
    if(verbose) then print("Combination 3 ", toString(1/2*inverseCofactorMatrix(Q, d)*inverse(transpose(B))*qNormMatrixFormat(B, Q, d, verbose))) else null;
    if(verbose) then print("Combination 4 ", toString((1/2)*inverseCofactorMatrix(Q, d)*inverse(transpose(B))*qNormMatrixFormat(B, Q, d, verbose))) else null;
    return 1/2*inverseCofactorMatrix(Q, d)*inverse(transpose(B))*qNormMatrixFormat(B, Q, d, verbose)
);
