-- barycenter of vectors which are columns of L
-- barycentre(L, d)
-- Input:
--   L: a matrix whose columns are vectors
--   d: dimension (number of columns in L minus 1)
-- Output:
--   The barycenter (average) of the vectors as a column vector
-- Description:
--   Computes the barycenter of d+1 vectors (columns of L) by averaging each coordinate.
-- Example:
--   barycentre(matrix{{1,2,3},{4,5,6}}, 2)
--   -- Output: matrix {{2}, {5}}
barycentre = (L, d)->(
    1/(d+1)*transpose(matrix{for i from 0 to d-1 list sum((entries(transpose(L))_i))})
);
-- print(toString barycentre(matrix{{1,2,3},{4,5,6}}, 2)) -- matrix {{2}, {5}}
