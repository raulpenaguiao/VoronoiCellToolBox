-- output is a list of vertices of a triangle suitable for input to secondMoment function
-- listMatrices is a list of matrices, each of which is a B matrix for a fromRelevantVectorsToVertex
-- VectorizedVertex(listMatrices, Q, d)
-- Input:
--   listMatrices: list of matrices (each for a vertex)
--   Q: symmetric matrix
--   d: dimension
-- Output:
--   Matrix of concatenated vertices
-- Description:
--   Computes and concatenates vertices from a list of relevant vector matrices.
-- Example:
--   VectorizedVertex({matrix{{1,0},{0,1}}, matrix{{0,1},{1,0}}}, matrix{{2,0},{0,3}}, 2)
VectorizedVertex = (listMatrices, Q, d) -> ( 
    concatVertices = fromRelevantVectorsToVertex(listMatrices_0, Q, d);
    for j from 1 to d do(
        concatVertices = concatVertices | fromRelevantVectorsToVertex(listMatrices_j, Q, d);
    );
    return concatVertices
);