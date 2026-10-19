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
--   Every vertex is v = adj(Q) w / (2 det Q) with w = B^(-T) (b_i^T Q b_i)_i linear in q, so
--   the second moment of a simplex with vertices v_0..v_d is
--     det(w_1-w_0 | ... | w_d-w_0) * (s^T adj(Q) s + sum_i w_i^T adj(Q) w_i) / (2^(d+2) (d+2)! det(Q)^2)
--   with s = sum_i w_i. All simplices share the denominator, so the numerators are summed as
--   polynomials and there is a single division at the end. Each simplex term gets the sign
--   that makes it positive at A (as makepos does).
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
    L := toList(0..d-1);
    adjQ := matrix table(d, d, (i, j) -> (-1)^(i+j) * det submatrix(Q, drop(L, {j, j}), drop(L, {i, i})));
    detQ := det Q;
    lvalues = Listify(A, d);
    atA := apply(G, k -> R_k => lvalues_k);
    wCache := new MutableHashTable;
    wOf := B -> (
        key := entries B;
        if not wCache#?key then wCache#key = sub(inverse(transpose(sub(B, QQ))), R) * qNormMatrixFormat(B, Q, d, false);
        wCache#key
    );
    denom := 2^(d+2) * (d+2)! * detQ^2;
    Znum := 0_R;
    l = # matVertices;
    for i from 0 to l-1 do(
        if(verboseComputationProgress) then print(concatenate(toString i, " of total ", toString l, " - new loop")) else null;
        W := apply(matVertices_i, wOf);
        if(verboseVertex) then print(concatenate("concatVertices = \n", toString((1/2*inverseCofactorMatrix(Q, d)) * fold((a, b) -> a | b, W)))) else null;
        s := sum W;
        volumeFactor := det fold((a, b) -> a | b, apply(drop(W, 1), w -> w - W_0));
        if sub(volumeFactor, atA) < 0 then volumeFactor = -volumeFactor;
        term := volumeFactor * ((transpose(s)*adjQ*s)_(0,0) + sum(W, w -> (transpose(w)*adjQ*w)_(0,0)));
        if(verboseSecondMoments) then print(concatenate("newPolyTriangle = \n", toString(term/denom))) else null;
        Znum = Znum + term;
    );
    Zpoly = Znum / denom;
    return Zpoly
);
