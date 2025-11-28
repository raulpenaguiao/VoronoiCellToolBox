-- symMatricesRing(d)
-- Input:
--   d: dimension
-- Output:
--   Polynomial ring with variables for symmetric d x d matrix
-- Description:
--   Creates a polynomial ring with variables for entries of a symmetric matrix, enforcing symmetry by identifying x_(i,j) with x_(j,i).
-- Example:
--   symMatricesRing(2)
--   -- Output: QQ[x_(1,1), x_(1,2), x_(2,2)]/(x_(1,2)-x_(2,1))
symMatricesRing = (d) -> (
    R = QQ[x_(1, 1)..x_(d, d)];
    I = ideal();
    for i from 1 to d do(
        for j from i+1 to d do(
            I = I+ideal(x_(i, j)-x_(j, i))
        );
    );
    return R/I
);
