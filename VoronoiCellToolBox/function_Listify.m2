-- listifying a matrix
-- Listify(matVertices, d)
-- Input:
--   matVertices: square matrix
--   d: dimension
-- Output:
--   List of values for symmetric matrix entries listed in an n*(n+1)/2 sized list,
--      values from the diagonal to the right
-- Description:
--   Produces a list encoding the matrix for use in substitutions.
-- Example:
--   Listify(matrix{{2, 1}, {1, 3}}, 2) -- Output: {2, 1, 3}
Listify = (matVertices, d) -> (
    lvalues = {};
    for i from 0 to d-1 do(
        for j from i to d-1 do(
            lvalues = append(lvalues, matVertices_(i,j))
        );
    );
    return lvalues
);
