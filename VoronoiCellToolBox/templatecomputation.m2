-- compute Q norms
qnorm = (v, Q)-> (
    return transpose(v)*Q*v
);

qnormmat = (B, Q, d) ->(
    return transpose( matrix{for i from 0 to d-1 list qnorm(matrix B_i, Q)})
);


-- compute inverse of dxd matrix
invCofactor = (A, d) ->(
    de = det A;
    mins = reverse (entries gens minors(d-1,A))_0; -- minors transpose
    answer = for i from 0 to d-1 list( for j from 0 to d-1 list mins_(j+d*i)*(-1)^(i+j)/de);
    return matrix(answer)
);

-- compute vertices
vertex = (B, Q, d)->(
    return 1/2*invCofactor(Q, d)*inverse(transpose(B))*qnormmat(B, Q, d)
);

-- barycenter of vectors which are columns of L
bary = (L, d)->(
    1/(d+1)*transpose(matrix{for i from 0 to d-1 list sum((entries(transpose(L))_i))})
);

--c reate ring for symmetric matrix variables
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

-- helper function
isSame = (i, j) ->(
    if i === j then return 1
    else return 0
);

-- voronoi cell is a generic d-dim permutohedron
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

-- listifying favorite matrix
Listify = (d) -> (
    lvalues = {};
    for i from 0 to d-1 do(
        for j from i to d-1 do(
            if i === j then lvalues = append(lvalues, d)
            else lvalues = append(lvalues, -1);
        );
    );
    return lvalues
);



-- it takes a polynomial, a list of values, and flips the sign of the polynomial if 
-- in the evaluation in these values it is negative
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
sm = (L, d, Q)->(
    T = submatrix'(L-matrix(for i from 0 to d list (L)_0),,{0});
    lis =  for i from 0 to d list qnorm(matrix(L_i), Q);
    return (det(T)*(qnorm((d+1)*bary(L, d), Q)+sum(lis))*(1/(d+2)!))_0_0
);

-- output is a list of vertices of a triangle suitable for input to sm function
-- listMatrices is a list of matrices, each of which is a B matrix for a vertex
VectorizedVertex = (listMatrices, Q, d) -> ( 
    concatVertices = vertex(listMatrices_0, Q, d);
    for j from 1 to d do(
        concatVertices = concatVertices | vertex(listMatrices_j, Q, d);
    );
    return concatVertices
);

-- Takes a dimension d and matVertices structure
-- matVertices encodes a triangulation of the voronoi cell
-- matvertices is a matrix, where each row corresponds to a triangle
-- in each row an entry corresponds to a vertex of the triangle
-- a vertex is encoded by a matrix
-- output is the rational function giving the second moment for any Voronoi cell in 
-- the same chamber as our favorite matrix in terms of the entries of the matrix Q
SmPoly = (d, matVertices) -> (
    -- A = favouriteMatrix(d);
    G = d*(d+1)//2;
    R = QQ[q_0..q_(G-1)];
    Q = genericSymmetricMatrix(R, q_0, d);
    Zpoly = 0;
    l = # matVertices;
    for i from 0 to l-1 do(
        if({{VERBOSE}}) then print(concatenate(toString i, " of total ", toString l, " - new loop")) else null;
        concatVertices = VectorizedVertex(matVertices_i, Q, d);
        if({{VERBOSE}}) then print(concatenate(toString i, " of total ", toString l, " - loop compute sm")) else null;
        ply = sm(concatVertices, d, Q);
        lvalues = Listify(d);
        if({{VERBOSE}}) then print(concatenate("l values = ", toString lvalues)) else null;
        newPolyTriangle = makepos(ply, lvalues, d, R);
        if({{VERBOSE}}) then print(concatenate("newPolyTriangle = ", toString newPolyTriangle)) else null;
        if({{VERBOSE}}) then print(concatenate(toString i, " of total ", toString l, " - loop add ")) else null;
        Zpoly = Zpoly + newPolyTriangle;
    );
    return Zpoly
);

mat = {{SAGESTRING}};
d = rank source  (mat_0)_0;
Zpoly = SmPoly(d, mat);
toString Zpoly