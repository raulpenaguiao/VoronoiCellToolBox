-- The second moment on a whole secondary cone, from a triangulation of the Voronoi cell, as a
-- rational function of the entries q_0, q_1, ... of the generic symmetric matrix.
-- SmPoly(d, matVertices, A)
-- Input:
--   d: dimension
--   matVertices: list of simplices; each is a list of d x d matrices B, one per vertex, whose
--     columns are d relevant vectors through that vertex
--   A: the metric whose secondary cone is wanted (fixes the orientation of each simplex)
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
--   SmPoly(2, {{matrix{{1,0},{0,1}}, matrix{{1,1},{0,1}}, matrix{{0,1},{1,1}}}}, matrix{{2,-1},{-1,2}})
SmPoly = (d, matVertices, A) -> (
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
    Znum := 0_R;
    for simplex in matVertices do (
        W := apply(simplex, wOf);
        s := sum W;
        volumeFactor := det fold((a, b) -> a | b, apply(drop(W, 1), w -> w - W_0));
        if sub(volumeFactor, atA) < 0 then volumeFactor = -volumeFactor;
        Znum = Znum + volumeFactor * ((transpose(s)*adjQ*s)_(0,0) + sum(W, w -> (transpose(w)*adjQ*w)_(0,0)));
    );
    Zpoly = Znum / (2^(d+2) * (d+2)! * detQ^2);
    return Zpoly
);
