-- voronoi cell is a generic d-dim permutohedron
-- favouriteMatrix(d)
-- Input:
--   d: dimension
-- Output:
--   d x d matrix with d on diagonal, -1 elsewhere
-- Description:
--   Constructs a matrix with d on the diagonal and -1 off-diagonal, used for Voronoi cell computations.
-- Example:
--   favouriteMatrix(3)
--   -- Output: matrix {{3,-1,-1},{-1,3,-1},{-1,-1,3}}
favouriteMatrix = (d) -> (
    ans = mutableMatrix map(QQ^d,QQ^d,(i,j)->isSame(i, j)*(d+1) - 1);
    for i from 0 to d-1 do(
        for j from 0 to d-1 do(
            ans_(i, j) = -1
        );
    );
    for i from 0 to d-1 do(
        ans_(i, i) = d
    );
    return matrix ans
);
