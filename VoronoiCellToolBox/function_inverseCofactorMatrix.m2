-- inverseCofactorMatrix(A, d)
-- Input:
--   A: a d x d matrix
--   d: dimension
-- Output:
--   The inverse of matrix A using the cofactor method
-- Description:
--   Computes the inverse of a square matrix using cofactors and determinant.
-- Example:
--  toString inverseCofactorMatrix(matrix{{2,1},{1,2}}, 2)
-- -- Output: matrix {{2/3, -1/3}, {-1/3, 2/3}}
--  toString inverseCofactorMatrix(matrix{{a, b}, {c, d}}, 2)
-- -- Output: matrix {{-d/(b*c-a*d), c/(b*c-a*d)}, {b/(b*c-a*d), -a/(b*c-a*d)}}
inverseCofactorMatrix = (A, d) ->(
    de = det A;
    mins = flatten (listOf1Minors(A, d));
    answer = for i from 0 to d-1 list( for j from 0 to d-1 list mins_(j+d*i)*(-1)^(i+j)/de);
    return matrix(answer)
);
