-- output is a list of vertices of a triangle suitable for input to secondMoment function
-- listMatrices is a list of matrices, each of which is a B matrix for a fromRelevantVectorsToVertex
-- VectorizedVertex(listMatrices, Q, d)
-- Input:
--   listMatrices: list of matrices (each for a vertex), may be empty
--   Q: symmetric matrix
--   d: dimension
-- Output:
--   Matrix of concatenated vertices; a d x 0 empty matrix if listMatrices is empty
-- Description:
--   Computes and concatenates vertices from a list of relevant vector matrices.
-- Example:
--   VectorizedVertex({matrix{{1,0},{0,1}}, matrix{{0,1},{1,0}}}, matrix{{2,0},{0,3}}, 2)
VectorizedVertex = (listMatrices, Q, d, verbose) -> (
    if #listMatrices == 0 then return map(QQ^d, QQ^0, {});
    concatVertices = fromRelevantVectorsToVertex(listMatrices_0, Q, d, verbose);
    for j from 1 to #listMatrices - 1 do(
        concatVertices = concatVertices | fromRelevantVectorsToVertex(listMatrices_j, Q, d, verbose);
    );
    return concatVertices
);
