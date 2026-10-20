import re
import warnings
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

def load_m2_template(file_name):
    """Return the text of a Macaulay2 file shipped with the package."""
    try:
        return files("VoronoiCellToolBox").joinpath(file_name).read_text()
    except Exception as e:
        raise FileNotFoundError(f"Could not find {file_name} in package") from e


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

_VERBOSE_FLAGS = ("verboseComputationProgress", "verboseSecondMoments", "verboseVertex")


def _ignore_verbose_flags(flags):
    """The old verbose flags printed Macaulay2 expressions for helpers that are not part of the
    package. They are accepted, and ignored, so that old calls keep working."""
    unknown = set(flags) - set(_VERBOSE_FLAGS)
    if unknown:
        raise TypeError(f"unexpected keyword arguments: {', '.join(sorted(unknown))}")
    if any(flags.values()):
        warnings.warn("the verbose flags of normalizedChamberSecondMomentPolynomial have no effect "
                      "any more; use debug_m2=True to see the raw Macaulay2 output",
                      DeprecationWarning, stacklevel=3)


def normalizedChamberSecondMomentPolynomial(Q, debug_m2=False, **verbose_flags):
    r"""
    Return the second moment on the secondary cone of `Q` as a string: a rational function of
    the upper-triangular entries ``q_0, q_1, ...`` of the metric, computed in Macaulay2. This is
    the engine of :func:`chamber_second_moment`, which returns the same function as an element of
    :func:`chamber_ring` and should be preferred.

    The triangulation is computed with the candidates `\{-2..2\}^d`. If they miss relevant
    vectors of `Q`, the LLL-reduced metric `Q_r = U^T Q U` is used instead, and `U^T Q U` is
    substituted back in Macaulay2, so that the result is a function of the entries of `Q` itself.

    INPUT:

    - ``Q`` -- a positive definite symmetric matrix
    - ``debug_m2`` -- (default: ``False``) print the raw Macaulay2 output

    EXAMPLES::

        sage: from VoronoiCellToolBox import normalizedChamberSecondMomentPolynomial
        sage: normalizedChamberSecondMomentPolynomial([[2, -1], [-1, 2]])
        '(2*q_0*q_1^2+2*q_1^3-q_0^2*q_2+2*q_1^2*q_2-q_0*q_2^2)/(12*q_1^2-12*q_0*q_2)'
    """
    _ignore_verbose_flags(verbose_flags)
    Q = _rows(Q)
    U = None
    try:
        VCell(Q, range=2)
    except Exception as e:
        if "relevant vector list not complete" not in str(e):
            raise
        Q, U = reduce_metric(Q)

    code = load_m2_function_with_dependencies("SmPoly") + """
mat = %s;
metric_matrix = %s;
d = numrows (mat_0)_0;
Zpoly = SmPoly(d, mat, metric_matrix);
""" % (FormatPullingTrigMatrix(Q).replace("\n", ""), macaulifyMatrix(Q))
    if U is not None:
        code += """
Ured = sub(%s, R);
Qback = transpose(Ured) * genericSymmetricMatrix(R, q_0, d) * Ured;
substitutions = apply(numgens R, k -> R_k => (Listify(Qback, d))_k);
Zpoly = sub(numerator Zpoly, substitutions) / sub(denominator Zpoly, substitutions);
""" % macaulifyMatrix([list(row) for row in U.rows()])
    code += "toString Zpoly\n"

    raw = macaulay2.eval(code)
    if debug_m2:
        print("=== Raw Macaulay2 output ===\n" + raw + "\n============================")
    if "--back trace--" in raw:
        raise RuntimeError(f"Macaulay2 computation failed (back trace):\n{raw}")
    # the last line holds the result, possibly after an output label "oo12 = "
    result = _M2_OUTPUT_LABEL_RE.sub('', raw.split('\n')[-1].strip())
    if not _M2_POLYNOMIAL_RE.fullmatch(result):
        raise RuntimeError(f"Macaulay2 output does not look like a rational function: {result!r}\n"
                           f"Full output:\n{raw}")
    return result


def normalizedChamberSecondMomentPolynomialResult(Q, debug_m2=False, **verbose_flags):
    """Deprecated alias of :func:`normalizedChamberSecondMomentPolynomial`."""
    return normalizedChamberSecondMomentPolynomial(Q, debug_m2=debug_m2, **verbose_flags)


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

