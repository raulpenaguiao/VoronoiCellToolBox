-- compute all 1-minors of a d x d matrix A
-- listOf1Minors(A, d)
-- Input:
--   A: a d x d matrix
--   d: dimension
-- Output:
--   A list of lists containing all 1-minors of A
-- Description:
--   Computes the matrix of cofactors (1-minors) of a square matrix A.
--   The (i,j)th entry is the determinant of the (d-1)x(d-1) submatrix corresponding to removing row i and column j.
-- Example:
-- R = QQ[a,b,c,d,e,f,g,h,i]
-- toString listOf1Minors(matrix{{a, b, c}, {d, e, f}, {g, h, i}}, 3);
--   -- Output: {{-f*h+e*i, -f*g+d*i, -e*g+d*h}, {-c*h+b*i, -c*g+a*i, -b*g+a*h}, {-c*e+b*f, -c*d+a*f, -b*d+a*e}}
listOf1Minors = (A, d) -> (
    L = toList(0..d-1);
    return for i from 0 to d-1 list (for j from 0 to d-1 list det submatrix(A, drop(L, {i, i}), drop(L, {j, j})));
);
