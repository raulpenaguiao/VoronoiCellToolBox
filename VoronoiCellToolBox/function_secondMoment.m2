-- compute second moment of simplex given by d+1 vertices as columns of L
-- actual second moment is also divided by dth root(det(Q))
-- make sure L has a fraction so we're in the right field
-- secondMoment(L, d, Q)
-- Input:
--   L: matrix of d+1 vertices (columns)
--   d: dimension
--   Q: symmetric matrix
-- Output:
--   Rational function for the second moment of the simplex
-- Description:
--   Computes the second moment of a simplex given by L, normalized by det(Q).
-- Example:
--   secondMoment(matrix{{1,2,3},{4,5,6}}, 2, matrix{{2,0},{0,3}})
--   -- Output: rational number
secondMoment = (L, d, Q)->(
    T = submatrix'(L-matrix(for i from 0 to d list (L)_0),,{0});
    lis =  for i from 0 to d list qnorm(matrix(L_i), Q);
    return (det(T)*(qnorm((d+1)*barycentre(L, d), Q)+sum(lis))*(1/(d+2)!))_0_0
);
