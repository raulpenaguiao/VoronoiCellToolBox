import re
import subprocess
from importlib.resources import files
from numbers import Integral as numbersIntegral
from sage.all import PolynomialRing, QQ
from sage.misc.parser import Parser
from sage.interfaces.macaulay2 import macaulay2
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix, macaulifyMatrix
from VoronoiCellToolBox.voronoi_cell import VCell, reduce_metric, is_generic, _rows

_M2_POLYNOMIAL_RE = re.compile(r'^[-+*/^()\d\sq_]+$')
_M2_OUTPUT_LABEL_RE = re.compile(r'^oo\d+\s*=\s*')

# Define function dependencies for proper loading order
FUNCTION_DEPENDENCIES = {
    "qnorm": [],
    "qNormMatrixFormat": ["qnorm"],
    "listOf1Minors": [],
    "inverseCofactorMatrix": ["listOf1Minors"],
    "isSame": [],
    "favouriteMatrix": ["isSame"],
    "fromRelevantVectorsToVertex": ["inverseCofactorMatrix", "qNormMatrixFormat", "qnorm"],
    "barycentre": [],
    "Listify": [],
    "secondMoment": ["qnorm", "barycentre"],
    "VectorizedVertex": ["fromRelevantVectorsToVertex"],
    "makepos": [],
    "SmPoly": ["Listify", "qNormMatrixFormat", "inverseCofactorMatrix"],
    "symMatricesRing": [],
}

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


def load_m2_function_with_dependencies(function_name):
    """
    Load a specific M2 function file along with all its dependencies in the correct order.

    Parameters:
    -----------
    function_name : str
        Name of the function to load (e.g., "barycentre", "SmPoly")

    Returns:
    --------
    str : Combined M2 code with the function and all its dependencies
    """
    if function_name not in FUNCTION_DEPENDENCIES:
        raise ValueError(f"Unknown function: {function_name}")

    # Use a set to track already loaded functions (avoid duplicates)
    loaded = set()
    code_parts = []

    def load_recursive(func_name):
        if func_name in loaded:
            return
        loaded.add(func_name)

        # Load dependencies first
        for dep in FUNCTION_DEPENDENCIES.get(func_name, []):
            load_recursive(dep)

        # Then load this function
        try:
            file_name = f"function_{func_name}.m2"
            func_content = load_m2_template(file_name)
            code_parts.append(func_content)
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find {file_name} for function {func_name}")

    load_recursive(function_name)
    return "\n".join(code_parts)


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
    # Step 1: Parse input vectors into Macaulay2 format
    Vmatrix = [[V[i][j] for j in range(len(V[0]))] for i in range(len(V))]
    string_vectors = "matrix " + str(Vmatrix).replace("]", "}").replace("[", "{")

    # Step 2: Load barycentre function (no dependencies)
    m2_input_string = load_m2_function_with_dependencies("barycentre")

    # Step 3: Add the code to compute barycentre
    barycentreCode = """
toString barycentre({{VECTORS}}, {{DIMENSION}})
    """.replace("{{VECTORS}}", string_vectors) \
       .replace("{{DIMENSION}}", str(len(V)))

    # Step 4: Run Macaulay2
    result = macaulay2.eval(m2_input_string + barycentreCode)

    # Step 5: Parse output back into Sage vector
    result = result.splitlines()[-1]
    result = result.replace("matrix {{", "").replace("}}", "")
    result = [[eval(str(x))] for x in result.split("}, {")]

    return result

def normalizedChamberSecondMomentPolynomialResult(Q, verboseComputationProgress=False, verboseSecondMoments=False, verboseVertex=False, debug_m2=False):
    """
    Execute a Sage computation and pass results to Macaulay2.

    Parameters:
    -----------
    Q : matrix
    debug_m2 : bool
        Print the full raw Macaulay2 output before any processing.

    Returns:
    --------
    the polynomial string with logs and debug info
    """
    result, _ = _normalizedChamberSecondMomentPolynomialRaw(Q, verboseComputationProgress, verboseSecondMoments, verboseVertex, debug_m2)
    return result


def _normalizedChamberSecondMomentPolynomialRaw(Q, verboseComputationProgress, verboseSecondMoments, verboseVertex, debug_m2):
    """
    Run SmPoly in Macaulay2 and return (polynomial string, full raw M2 output).
    The raw output is needed by the verbose modes, which print sections of it.
    """
    # Step 1: Run Sage computation
    # The triangulation is computed with range=2 candidate relevant vectors. If that misses
    # some relevant vectors of Q, work with the LLL-reduced metric Qr = U^T Q U instead and
    # substitute back in Macaulay2: SM_Q = SM_Qr(U^T Q U) on the whole chamber of Q.
    U = None
    try:
        VCell(Q, range=2)
    except Exception as e:
        if "relevant vector list not complete" not in str(e):
            raise
        Q, U = reduce_metric(Q)
    sage_string = FormatPullingTrigMatrix(Q)
    matrix_m2 = macaulifyMatrix(Q)

    # Step 2: Load SmPoly function with all its dependencies
    m2_input_string = load_m2_function_with_dependencies("SmPoly")

    secondMomentCode = """
mat = {{SAGESTRING}};
metric_matrix = {{SAGESTRING2}};
d = numrows (mat_0)_0;
Zpoly = SmPoly(d, mat, metric_matrix, {{verboseComputationProgress}}, {{verboseSecondMoments}}, {{verboseVertecExpressions}});
{{BASISCHANGE}}toString Zpoly
"""
    basisChangeCode = ""
    if U is not None:
        basisChangeCode = """
Ured = sub({{UMATRIX}}, R);
Qback = transpose(Ured) * genericSymmetricMatrix(R, q_0, d) * Ured;
substitutions = apply(numgens R, k -> R_k => (Listify(Qback, d))_k);
Zpoly = sub(numerator Zpoly, substitutions) / sub(denominator Zpoly, substitutions);
""".replace("{{UMATRIX}}", macaulifyMatrix([list(row) for row in U.rows()]))
    secondMomentCode = secondMomentCode.replace("{{BASISCHANGE}}", basisChangeCode)

    # Parse m2 for this function
    m2_input_string += secondMomentCode
    m2_input_string = m2_input_string.replace("{{SAGESTRING}};", sage_string.replace("\n", "") + ";")
    m2_input_string = m2_input_string.replace("{{SAGESTRING2}}", matrix_m2)
    if verboseComputationProgress:
        m2_input_string = m2_input_string.replace("{{verboseComputationProgress}}", "true")
    else:
        m2_input_string = m2_input_string.replace("{{verboseComputationProgress}}", "false")

    if verboseSecondMoments:
        m2_input_string = m2_input_string.replace("{{verboseSecondMoments}}", "true")
    else:
        m2_input_string = m2_input_string.replace("{{verboseSecondMoments}}", "false")

    if verboseVertex:
        m2_input_string = m2_input_string.replace("{{verboseVertecExpressions}}", "true")
    else:
        m2_input_string = m2_input_string.replace("{{verboseVertecExpressions}}", "false")

    # Step 3: Run Macaulay2
    raw_result = macaulay2.eval(m2_input_string)

    # C: debug flag — show raw output before any processing
    if debug_m2:
        print("=== Raw Macaulay2 output ===\n" + raw_result + "\n============================")

    # A: detect known M2 fatal-error markers
    if "--back trace--" in raw_result:
        raise RuntimeError(f"Macaulay2 computation failed (back trace):\n{raw_result}")

    result = raw_result.split('\n')[-1].strip()
    # When M2 printed anything, the last line carries an output label: "oo1000000018 = (...)/(...)"
    result = _M2_OUTPUT_LABEL_RE.sub('', result)

    # B: validate output looks like a rational function in q_N variables
    if result and not _M2_POLYNOMIAL_RE.fullmatch(result):
        raise RuntimeError(
            f"Macaulay2 output does not look like a polynomial/rational function: {result!r}\n"
            f"Full output:\n{raw_result}"
        )

    return result, raw_result


def normalizedChamberSecondMomentPolynomial(Q, verboseComputationProgress=False, verboseSecondMoments=False, verboseVertex=False, debug_m2=False):
    Q = _rows(Q)
    result, raw_result = _normalizedChamberSecondMomentPolynomialRaw(Q, verboseComputationProgress, verboseSecondMoments, verboseVertex, debug_m2)
    if verboseComputationProgress:
        print("Debug computation progress: result = " + str(result))
    elif verboseSecondMoments:
        lines = raw_result.splitlines()
        secondMoment = []
        for i, line in enumerate(lines):
            if "newPolyTriangle = " in line and i + 1 < len(lines):
                secondMoment.append(f' sm{i} = "' + lines[i + 1] + '"')
                secondMoment.append(f'print(rpl(sm{i}, {Q}))')
        print("Debug second moment expressions: \n " + "\n".join(secondMoment))
    elif verboseVertex:
        #find all instances of "concatVertices = " and print the line after it
        lines = raw_result.splitlines()
        vertexCoordinates = []
        for i, line in enumerate(lines):
            if "concatVertices = " in line and i + 1 < len(lines):
                vertexCoordinates.append(f' smtriang{i} = "' + lines[i + 1] + '"')
                vertexCoordinates.append(f'print(rpl(m2_matrix_to_sage(smtriang{i}), {Q}))')
        print("Debug vertex expressions: \n " + "\n".join(vertexCoordinates))
    if '=' in result:
        result = result.split('=')[-1].strip()
    return result  # return only the last line which contains the polynomial


def chamber_ring(d):
    r"""
    Return the fraction field `\QQ(q_0, \dots, q_{N-1})`, `N = d(d+1)/2`, in which
    :func:`chamber_second_moment` returns its result. The variable `q_k` is the `k`-th
    upper-triangular entry of the metric in row-major order: for `d = 2`,
    `Q = [[q_0, q_1], [q_1, q_2]]`.

    EXAMPLES::

        sage: from VoronoiCellToolBox import chamber_ring
        sage: chamber_ring(2)
        Fraction Field of Multivariate Polynomial Ring in q_0, q_1, q_2 over Rational Field
    """
    return PolynomialRing(QQ, ['q_%d' % k for k in range(d*(d+1)//2)]).fraction_field()


def metric_coordinates(Q):
    r"""
    Return the upper-triangular entries of `Q` in row-major order, the point at which a
    chamber function in :func:`chamber_ring` is evaluated: ``f(*metric_coordinates(Q))``.

    EXAMPLES::

        sage: from VoronoiCellToolBox import metric_coordinates
        sage: metric_coordinates([[2, -1], [-1, 3]])
        [2, -1, 3]
    """
    d = len(Q)
    return [QQ(Q[i][j]) for i in range(d) for j in range(i, d)]


def chamber_second_moment(Q, allow_nongeneric=False):
    r"""
    Return the second moment on the whole secondary cone of `Q`, as a rational function of
    the entries of the metric.

    For every `Q'` in the secondary cone of `Q` (:func:`secondary_cone`), the second moment
    ``second_moment(Q')`` equals ``f(*metric_coordinates(Q'))``, where ``f`` is the output.
    The vertices of the Voronoi cell are rational functions of `Q'` on the cone, so the
    integral over a fixed triangulation is one too; the symbolic integration runs in
    Macaulay2. Metrics with relevant vectors outside `\{-2..2\}^d` are LLL-reduced first and
    the result is substituted back.

    The normalized second moment on the cone is ``f / (d * det(Q')^(1/d))``.

    INPUT:

    - ``Q`` -- a positive definite symmetric matrix
    - ``allow_nongeneric`` -- (default: ``False``) a non-generic `Q` lies on a wall between
      chambers, and its formula is only valid on that face. By default a ``ValueError`` is
      raised; with ``True`` the face formula is returned.

    OUTPUT: an element of :func:`chamber_ring`

    EXAMPLES::

        sage: from VoronoiCellToolBox import chamber_second_moment, metric_coordinates, second_moment
        sage: f = chamber_second_moment([[2, -1], [-1, 2]]); f
        (2*q_0*q_1^2 + 2*q_1^3 - q_0^2*q_2 + 2*q_1^2*q_2 - q_0*q_2^2)/(12*q_1^2 - 12*q_0*q_2)
        sage: Qp = [[3, -1], [-1, 5]]                       # another metric in the same chamber
        sage: f(*metric_coordinates(Qp)) == second_moment(Qp)
        True
        sage: chamber_second_moment([[1, 0], [0, 1]])
        Traceback (most recent call last):
        ...
        ValueError: Q is not generic (it lies on a wall between chambers); pass allow_nongeneric=True for the formula on that face
    """
    Q = _rows(Q)
    if not allow_nongeneric and not is_generic(Q):
        raise ValueError("Q is not generic (it lies on a wall between chambers); pass allow_nongeneric=True for the formula on that face")
    R = chamber_ring(len(Q))
    result = normalizedChamberSecondMomentPolynomial(Q)
    # Sage's parser, not sage_eval: Python's compiler overflows its stack on long sums
    return R(Parser(make_var={str(g): g for g in R.gens()}).parse(result))



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

    m2_input_string = load_m2_function_with_dependencies("qnorm")

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

    m2_input_string = load_m2_function_with_dependencies("qNormMatrixFormat")

    code = f"""
B = {B_m2};
Q = {Q_m2};
d = {d};
result = qNormMatrixFormat(B, Q, d, false);
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

    m2_input_string = load_m2_function_with_dependencies("listOf1Minors")

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

    m2_input_string = load_m2_function_with_dependencies("inverseCofactorMatrix")

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

    m2_input_string = load_m2_function_with_dependencies("fromRelevantVectorsToVertex")

    code = f"""
B = {B_m2};
Q = {Q_m2};
d = {d};
result = fromRelevantVectorsToVertex(B, Q, d, false);
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

    m2_input_string = load_m2_function_with_dependencies("barycentre")

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

    m2_input_string = load_m2_function_with_dependencies("secondMoment")

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

    m2_input_string = load_m2_function_with_dependencies("VectorizedVertex")

    code = f"""
listMatrices = {list_m2};
Q = {Q_m2};
d = {d};
result = VectorizedVertex(listMatrices, Q, d, false);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def isSame_m2(i, j):
    """
    Check if two indices are the same using Macaulay2.

    Parameters:
    -----------
    i : int
        First index
    j : int
        Second index

    Returns:
    --------
    int : 1 if i == j, else 0
    """
    # Input validation
    # TODO

    m2_input_string = load_m2_function_with_dependencies("isSame")

    code = f"""
result = isSame({i}, {j});
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    result_str = str(result.splitlines()[-1]).strip()
    try:
        return int(result_str)
    except ValueError:
        raise ValueError(f"Could not parse M2 output as integer: {result_str}")


def favouriteMatrix_m2(d):
    """
    Compute the favorite matrix for Voronoi cell computations using Macaulay2.

    Parameters:
    -----------
    d : int
        Dimension

    Returns:
    --------
    str : The favorite matrix with d on diagonal and -1 elsewhere
    """
    m2_input_string = load_m2_function_with_dependencies("favouriteMatrix")

    code = f"""
result = favouriteMatrix({d});
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def Listify_m2(matVertices, d):
    """
    Encode a symmetric matrix into a list of values using Macaulay2.

    Parameters:
    -----------
    matVertices : matrix
        A symmetric matrix to encode
    d : int
        Dimension

    Returns:
    --------
    str : List representation of symmetric matrix entries
    """
    matVertices_m2 = macaulifyMatrix(matVertices)

    m2_input_string = load_m2_function_with_dependencies("Listify")

    code = f"""
matVertices = {matVertices_m2};
d = {d};
result = Listify(matVertices, d);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def makepos_m2(polynom_str, lvalues, d, ring_str):
    """
    Adjust polynomial sign to ensure positivity at given values using Macaulay2.

    Parameters:
    -----------
    polynom_str : str
        Polynomial expression as a string (e.g., "q_0 + q_1")
    lvalues : list
        List of values for substitution
    d : int
        Dimension
    ring_str : str
        Polynomial ring definition as a string (e.g., "QQ[q_0, q_1]")

    Returns:
    --------
    str : The polynomial, possibly with sign flipped to ensure positivity at lvalues
    """
    # Input validation
    if not isinstance(polynom_str, str):
        raise TypeError("polynom_str must be a string")
    if not isinstance(ring_str, str):
        raise TypeError("ring_str must be a string")
    if not isinstance(lvalues, (list, tuple)):
        raise TypeError("lvalues must be a list or tuple")
    if not isinstance(d, numbersIntegral):
        raise TypeError("d must be an integer")

    m2_input_string = load_m2_function_with_dependencies("makepos")

    # Convert lvalues list to Macaulay2 format
    lvalues_str = "{" + ", ".join(str(v) for v in lvalues) + "}"

    code = f"""
R = {ring_str};
polynom = {polynom_str};
lvalues = {lvalues_str};
d = {d};
result = makepos(polynom, lvalues, d, R);
toString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])


def symMatricesRing_m2(d):
    """
    Create a polynomial ring for symmetric d x d matrices using Macaulay2.

    Parameters:
    -----------
    d : int
        Dimension

    Returns:
    --------
    str : String representation of the symmetric matrix ring with symmetry relations
    """
    m2_input_string = load_m2_function_with_dependencies("symMatricesRing")

    code = f"""
result = symMatricesRing({d});
toExternalString result
"""

    m2_input_string += code
    result = macaulay2.eval(m2_input_string)
    return str(result.splitlines()[-1])

