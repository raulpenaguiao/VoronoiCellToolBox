-- listifying favorite matrix
-- Listify(d)
-- Input:
--   matVertices: matrix with 
-- Output:
--   List of values for symmetric matrix entries: d for diagonal, -1 for off-diagonal
-- Description:
--   Produces a list encoding the favorite matrix for use in substitutions.
-- Example:
--   Listify(2) -- Output: {2, -1, 2}
Listify = (matVertices) -> (
    lvalues = {};
    for i from 0 to d-1 do(
        for j from i to d-1 do(
            if i === j then lvalues = append(lvalues, d)
            else lvalues = append(lvalues, -1);
        );
    );
    return lvalues
);
-- print(toString Listify(2)); -- {2, -1, 2}