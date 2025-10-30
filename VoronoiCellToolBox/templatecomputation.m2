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
-- compute Q norms
qnorm = (v, Q)-> (
    return transpose(v)*Q*v
);
-- print(qnorm(matrix{{1},{2}}, matrix{{2,0},{0,3}}));


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


-- compute all 1-minors of a d x d matrix A
-- listOf1Minors(A, d)
-- Input:
--   A: a d x d matrix
--   d: dimension
-- Output:
--   A list of lists containing all 1-minors of A
-- Description:
--   Computes the matrix of cofactors (1-minors) of a square matrix A.
--   The (i,j)th entry is the determinant of the (d-1)x(d-1) submatrix corresponding to removing row i and column j.
-- Example:
-- R = QQ[a,b,c,d,e,f,g,h,i]
-- toString listOf1Minors(matrix{{a, b, c}, {d, e, f}, {g, h, i}}, 3);
--   -- Output: {{-f*h+e*i, -f*g+d*i, -e*g+d*h}, {-c*h+b*i, -c*g+a*i, -b*g+a*h}, {-c*e+b*f, -c*d+a*f, -b*d+a*e}}
listOf1Minors = (A, d) -> (
    L = toList(0..d-1);
    return for i from 0 to d-1 list (for j from 0 to d-1 list det submatrix(A, drop(L, {i, i}), drop(L, {j, j})));
);
-- R = QQ[a, b, c, d, e, f, g, h, i];
-- print(toString listOf1Minors(matrix{{a, b, c}, {d, e, f}, {g, h, i}}, 3));


-- inverseCofactorMatrix(A, d)
-- Input:
--   A: a d x d matrix
--   d: dimension
-- Output:
--   The inverse of matrix A using the cofactor method
-- Description:
--   Computes the inverse of a square matrix using cofactors and determinant.
-- Example:
--  toString inverseCofactorMatrix(matrix{{2,1},{1,2}}, 2) 
-- -- Output: matrix {{2/3, -1/3}, {-1/3, 2/3}}
--  toString inverseCofactorMatrix(matrix{{a, b}, {c, d}}, 2)
-- -- Output: matrix {{-d/(b*c-a*d), c/(b*c-a*d)}, {b/(b*c-a*d), -a/(b*c-a*d)}}
inverseCofactorMatrix = (A, d) ->(
    de = det A;
    mins = flatten (listOf1Minors(A, d)); -- minors transpose
    answer = for i from 0 to d-1 list( for j from 0 to d-1 list mins_(j+d*i)*(-1)^(i+j)/de);
    return matrix(answer)
);
-- R = QQ[a, b, c, d];
-- print(toString inverseCofactorMatrix(matrix{{2,1},{1,2}}, 2));
-- print(toString inverseCofactorMatrix(matrix{{a, b}, {c, d}}, 2));


-- fromRelevantVectorsToVertex(B, Q, d)
-- Input:
--   B: matrix of relevant vectors, as column vectors
--   Q: symmetric matrix
--   d: dimension
-- Output:
--   Vertex coordinates as a column vector and as a polynomial in entries of Q
-- Description:
--   Computes the vertex of a Voronoi cell from relevant vectors and Q.
-- Example:
--   toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{2,0},{0,3}}, 2)
-- -- Output: matrix {{1/2}, {1/2}}
--   toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{q11, q12},{q12,q22}}, 2)
-- -- Output: matrix {{(-q11*q22+q12*q22)/(2*q12^2-2*q11*q22)}, {(q11*q12-q11*q22)/(2*q12^2-2*q11*q22)}}
fromRelevantVectorsToVertex = (B, Q, d)->(
    return 1/2*inverseCofactorMatrix(Q, d)*inverse(transpose(B))*qNormMatrixFormat(B, Q, d)
);
-- print(toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{2,0},{0,3}}, 2));
-- R = QQ[q11, q12, q22];
-- print(toString fromRelevantVectorsToVertex(matrix{{1,0},{0,1}}, matrix{{q11, q12},{q12,q22}}, 2));

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
-- print(symMatricesRing(4));

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
-- print(toString isSame(1,1)); -- 1
-- print(toString isSame(1,2)); -- 0

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
-- print(toString favouriteMatrix(3)); -- matrix {{3,-1,-1},{-1,3,-1},{-1,-1,3}}


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

-- it takes a polynomial, a list of values, and flips the sign of the polynomial if 
-- in the evaluation in these values it is negative
-- makepos(polynom, lvalues, d, RingR)
-- Input:
--   polynom: a polynomial (rational function)
--   lvalues: list of values for substitution
--   d: dimension
--   RingR: polynomial ring
-- Output:
--   The polynomial, possibly with sign flipped to ensure positivity at lvalues
-- Description:
--   Substitutes lvalues into polynom and flips sign if result is negative.
-- Example:
--   makepos(q_0+q_1, {2,-1}, 2, QQ[q_0,q_1])
--   -- Output: q_0+q_1 or -(q_0+q_1) depending on sign
makepos = (polynom, lvalues, d, RingR) -> (
    G = d*(d+1)//2;
    slst = {};
    genList = gens RingR;
    for i from 0 to G-1 do (
        slst = append(slst, genList_i => lvalues_i); 
    );
    subvalue1 = sub(numerator polynom, slst);
    subvalue2 = sub(denominator polynom, slst);
    subvalue = subvalue1/subvalue2;
    if subvalue < 0 then return -polynom
    else return polynom
);


-- compute second moment of simplex given by d+1 vertices as columns of L
-- actual second moment is also divided by dth root(det(Q))
-- make sure L has a fraction so we're in the right field
-- secondMoment(L, d, Q)
-- Input:
--   L: matrix of d+1 vertices (columns)
--   d: dimension
--   Q: symmetric matrix
-- Output:
--   Rational function for the second moment of the simplex
-- Description:
--   Computes the second moment of a simplex given by L, normalized by det(Q).
-- Example:
--   secondMoment(matrix{{1,2,3},{4,5,6}}, 2, matrix{{2,0},{0,3}})
--   -- Output: rational number
secondMoment = (L, d, Q)->(
    T = submatrix'(L-matrix(for i from 0 to d list (L)_0),,{0});
    lis =  for i from 0 to d list qnorm(matrix(L_i), Q);
    return (det(T)*(qnorm((d+1)*barycentre(L, d), Q)+sum(lis))*(1/(d+2)!))_0_0
);


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


-- Takes a dimension d and matVertices structure
-- matVertices encodes a triangulation of the voronoi cell
-- matvertices is a matrix, where each row corresponds to a triangle
-- in each row an entry corresponds to a fromRelevantVectorsToVertex of the triangle
-- a fromRelevantVectorsToVertex is encoded by a matrix
-- output is the rational function giving the second moment for any Voronoi cell in 
-- the same chamber as our favorite matrix in terms of the entries of the matrix Q
-- SmPoly(d, matVertices)
-- Input:
--   d: dimension
--   matVertices: list of matrices encoding triangulation
--   A: metric symmetric matrix
-- Output:
--   Rational function for the second moment polynomial of the Voronoi cell
-- Description:
--   Computes the second moment polynomial for a Voronoi cell using its triangulation.
-- Example:
--   SmPoly(2, { {matrix{{1,0},{0,1}}, matrix{{0,1},{1,0}}} })
SmPoly = (d, matVertices, A, verbose) -> (
    G = d*(d+1)//2;
    R = QQ[q_0..q_(G-1)];
    Q = genericSymmetricMatrix(R, q_0, d);
    Zpoly = 0;
    lvalues = Listify(A);
    l = # matVertices;
    for i from 0 to l-1 do(
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - new loop")) else null;
        concatVertices = VectorizedVertex(matVertices_i, Q, d);
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - loop compute secondMoment")) else null;
        ply = secondMoment(concatVertices, d, Q);
        if(verbose) then print(concatenate("l values = ", toString lvalues)) else null;
        newPolyTriangle = makepos(ply, lvalues, d, R);
        if(verbose) then print(concatenate("newPolyTriangle = ", toString newPolyTriangle)) else null;
        if(verbose) then print(concatenate(toString i, " of total ", toString l, " - loop add ")) else null;
        Zpoly = Zpoly + newPolyTriangle;
    );
    return Zpoly
);
