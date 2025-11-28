-- qnorm(v, Q)
-- Input:
--   v: a column vector (matrix)
--   Q: a symmetric matrix
-- Output:
--   The quadratic form vᵗ Q v (a scalar)
-- Description:
--   Computes the Q-norm (quadratic form) of vector v with respect to matrix Q.
-- Example:
--   qnorm(matrix{{1},{2}}, matrix{{2,0},{0,3}})
-- Output: 1*2*1 + 2*3*2 = 2 + 12 = 14
qnorm = (v, Q)-> (
    return transpose(v)*Q*v
);
