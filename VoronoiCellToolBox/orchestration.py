import subprocess
from importlib.resources import files
from sage.interfaces.macaulay2 import macaulay2
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix, macaulifyMatrix

def load_m2_template(template_name="templatecomputation.m2"):
    """
    Load the templatecomputation.m2 file from your package
    """
    try:
        template_content = files("VoronoiCellToolBox") \
            .joinpath(template_name) \
            .read_text()
        
        return template_content
    except Exception as e:
        print(f"Error loading resource: {e}")
        raise FileNotFoundError(f"Could not find {template_name} in package")


def barycentre(V, verbose=False):
    """
    Uses macaulay2 to compute the barycentre of a collection of vectors.
    
    Parameters:
    -----------
    V : a list of vectors
        
    Returns:
    --------
    A vector representing the barycentre of the list of vectors V.
    """
    #print("\nEntering barycentre function")
    # Step 1: Parse input vectors into Macaulay2 format
    Vmatrix = [[V[i][j] for j in range(len(V[0]))] for i in range(len(V))]#We don't assume V is a matrix
    string_vectors = "matrix " + str(Vmatrix).replace("]", "}").replace("[", "{")

    # Step 2: Prepare Macaulay2 input file
    m2_input_string = load_m2_template("function_barycentre.m2")

    #Step 3: Add the code to compute barycentre
    barycentreCode = """
toString barycentre({{VECTORS}}, {{DIMENSION}})
    """.replace("{{VECTORS}}", string_vectors) \
       .replace("{{DIMENSION}}", str(len(V))) \
       .replace("{{VERBOSE}}", "true" if verbose else "false")
    
    #print("Debug 3: m2_input_string = " + m2_input_string)
    #print("Debug 4: barycentreCode = " + barycentreCode)

    #Step 4: Run Macaulay2
    result = macaulay2.eval(m2_input_string + barycentreCode)

    #Step 5: Parse output back into Sage vector
    result = result.splitlines()[-1]  # get only the last line
    result = result.replace("matrix {{", "").replace("}}", "")
    result = [[eval(str(x))] for x in result.split("}, {")]

    return result

def normalizedChamberSecondMomentPolynomial(Q, verbose=False):
    """
    Execute a Sage computation and pass results to Macaulay2.

    Parameters:
    -----------
    Q : matrix

    Returns:
    --------
    the string displaying the polynomial of the second moment of the Voronoi cell
    for any matrix P in the same secondary cone as Q.
    """
    # Step 1: Run Sage computation
    sage_string = FormatPullingTrigMatrix(Q)
    matrix_m2 = macaulifyMatrix(Q)
    if(verbose):
        print("Debug 1: matrix_m2 = " + matrix_m2)
    
    # Step 2: Prepare Macaulay2 input file
    m2_input_string = load_m2_template()
    #print("Debug 2: m2_input_string = " + m2_input_string)

    secondMomentCode = """
mat = {{SAGESTRING}};
metric_matrix = {{SAGESTRING2}};
d = numrows (mat_0)_0;
Zpoly = SmPoly(d, mat, metric_matrix, {{VERBOSE}});
toString Zpoly
"""

    #parse m2 for this function
    m2_input_string += secondMomentCode
    m2_input_string = m2_input_string.replace("{{SAGESTRING}};", sage_string.replace("\n", "") + ";")
    m2_input_string = m2_input_string.replace("{{SAGESTRING2}}", matrix_m2 )
    if( verbose ):
        m2_input_string = m2_input_string.replace("{{VERBOSE}}", "true")
    else:
        m2_input_string = m2_input_string.replace("{{VERBOSE}}", "false")
    print("Debug 3: m2_input_string = " + m2_input_string)
    # Step 3: Run Macaulay2
    result = macaulay2.eval(m2_input_string)
    if( verbose ):
        print("Debug 4: result = " + str(result))
    return str(result.splitlines()[-1])  # return only the last line which contains the polynomial


def qnorm_m2(v, Q):
    """
    Compute the Q-norm (quadratic form) of vector v with respect to matrix Q using Macaulay2.

    Parameters:
    -----------
    v : list or vector
        A column vector
    Q : matrix
        A symmetric matrix

    Returns:
    --------
    str : The quadratic form v^T Q v as a polynomial expression
    """
    v_m2 = macaulifyMatrix([[vi] for vi in v])
    Q_m2 = macaulifyMatrix(Q)

    m2_input_string = load_m2_template()

    code = f"""
v = {v_m2};
Q = {Q_m2};
result = qnorm(v, Q);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def qNormMatrixFormat_m2(B, Q):
    """
    Compute Q-norm for each column vector in B using Macaulay2.

    Parameters:
    -----------
    B : matrix
        A matrix whose columns are the desired vectors
    Q : matrix
        A symmetric matrix

    Returns:
    --------
    str : Row vector of Q-norms for each column of B
    """
    d = len(B[0])  # number of columns in B
    B_m2 = macaulifyMatrix(B)
    Q_m2 = macaulifyMatrix(Q)

    m2_input_string = load_m2_template()

    code = f"""
B = {B_m2};
Q = {Q_m2};
d = {d};
result = qNormMatrixFormat(B, Q, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def listOf1Minors_m2(A):
    """
    Compute all 1-minors (cofactors) of a d x d matrix A using Macaulay2.

    Parameters:
    -----------
    A : matrix
        A d x d square matrix

    Returns:
    --------
    str : A list of lists containing all 1-minors of A
    """
    d = len(A)
    A_m2 = macaulifyMatrix(A)

    m2_input_string = load_m2_template()

    code = f"""
A = {A_m2};
d = {d};
result = listOf1Minors(A, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def inverseCofactorMatrix_m2(A):
    """
    Compute the inverse of matrix A using the cofactor method in Macaulay2.

    Parameters:
    -----------
    A : matrix
        A square d x d matrix

    Returns:
    --------
    str : The inverse of matrix A as a string representation
    """
    d = len(A)
    A_m2 = macaulifyMatrix(A)

    m2_input_string = load_m2_template()

    code = f"""
A = {A_m2};
d = {d};
result = inverseCofactorMatrix(A, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def fromRelevantVectorsToVertex_m2(B, Q):
    """
    Compute vertex coordinates from relevant vectors using Macaulay2.

    Parameters:
    -----------
    B : matrix
        Matrix of relevant vectors as column vectors
    Q : matrix
        Symmetric matrix

    Returns:
    --------
    str : Vertex coordinates as a column vector (polynomial in entries of Q)
    """
    d = len(Q)
    B_m2 = macaulifyMatrix(B)
    Q_m2 = macaulifyMatrix(Q)

    m2_input_string = load_m2_template()

    code = f"""
B = {B_m2};
Q = {Q_m2};
d = {d};
result = fromRelevantVectorsToVertex(B, Q, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def barycentre_m2(L):
    """
    Compute the barycenter of vectors (columns of L) using Macaulay2.

    Parameters:
    -----------
    L : matrix
        A matrix whose columns are vectors

    Returns:
    --------
    str : The barycenter (average) of the vectors as a column vector
    """
    d = len(L[0]) - 1  # d is number of columns minus 1
    L_m2 = macaulifyMatrix(L)

    m2_input_string = load_m2_template()

    code = f"""
L = {L_m2};
d = {d};
result = barycentre(L, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def secondMoment_m2(L, Q):
    """
    Compute the second moment of a simplex given by d+1 vertices using Macaulay2.

    Parameters:
    -----------
    L : matrix
        Matrix of d+1 vertices (columns)
    Q : matrix
        Symmetric matrix

    Returns:
    --------
    str : Rational function for the second moment of the simplex
    """
    d = len(L) - 1
    L_m2 = macaulifyMatrix(L)
    Q_m2 = macaulifyMatrix(Q)

    m2_input_string = load_m2_template()

    code = f"""
L = {L_m2};
Q = {Q_m2};
d = {d};
result = secondMoment(L, d, Q);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def VectorizedVertex_m2(listMatrices, Q):
    """
    Compute and concatenate vertices from a list of relevant vector matrices using Macaulay2.

    Parameters:
    -----------
    listMatrices : list of matrices
        List of matrices, each representing relevant vectors for a vertex
    Q : matrix
        Symmetric matrix

    Returns:
    --------
    str : Matrix of concatenated vertices
    """
    d = len(Q)
    Q_m2 = macaulifyMatrix(Q)

    # Convert list of matrices to Macaulay2 format
    list_m2 = "{" + ", ".join([macaulifyMatrix(mat) for mat in listMatrices]) + "}"

    m2_input_string = load_m2_template()

    code = f"""
listMatrices = {list_m2};
Q = {Q_m2};
d = {d};
result = VectorizedVertex(listMatrices, Q, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])

