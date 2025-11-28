-- helper function
-- isSame(i, j)
-- Input:
--   i, j: integers
-- Output:
--   1 if i == j, else 0
-- Description:
--   Returns 1 if the two indices are the same, 0 otherwise.
-- Example:
--   isSame(1,1) -- Output: 1
--   isSame(1,2) -- Output: 0
isSame = (i, j) ->(
    if i === j then return 1
    else return 0
);
