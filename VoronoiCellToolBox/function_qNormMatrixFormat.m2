-- qNormMatrixFormat(B, Q, d)
-- Input:
--   B: a matrix whose columns are the desired vectors
--   Q: a symmetric matrix
--   d: dimension (number of columns in B)
-- Output:
--   Row vector of Q-norms for each column of B
-- Description:
--   Computes the Q-norm for each column vector in B and returns as a row matrix.
-- Example:
--   qNormMatrixFormat(matrix{{1,1},{0,1}}, matrix{{2,0},{0,3}}, 2)
--   -- Output: matrix {{2}, {5}}
qNormMatrixFormat = (B, Q, d) ->(
    return transpose( matrix{for i from 0 to d-1 list qnorm(matrix B_i, Q)})
);
-- print(toString qNormMatrixFormat(matrix{{1,1},{0,1}}, matrix{{2,0},{0,3}}, 2));