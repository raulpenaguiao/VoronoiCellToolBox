import subprocess
from importlib.resources import files
from sage.interfaces.macaulay2 import macaulay2
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix, macaulifyMatrix

def load_m2_template():
    """
    Load the templatecomputation.m2 file from your package
    """
    try:
        template_content = files("VoronoiCellToolBox") \
            .joinpath("templatecomputation.m2") \
            .read_text()
        
        return template_content
    except Exception as e:
        print(f"Error loading resource: {e}")
        raise FileNotFoundError("Could not find templatecomputation.m2 in package")


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

    # Step 2: Prepare Macaulay2 input file
    m2_input_string = load_m2_template()
    #print("Debug 2: m2_input_string = " + m2_input_string)

    secondMomentCode = """
mat = {{SAGESTRING}};
print(1 + 1);
metric_matrix = {{SAGESTRING2}};
print(1 + 2);
d = numrows (mat_0)_0;
print(1 + 3);
Zpoly = SmPoly(d, mat, metric_matrix, {{VERBOSE}});
toString Zpoly
"""

    #parse m2 for this function
    m2_input_string += secondMomentCode
    m2_input_string = m2_input_string.replace("{{SAGESTRING}};", sage_string.replace("\n", "") + ";")
    m2_input_string = m2_input_string.replace("{{SAGESTRING2}}", matrix_m2 )
    m2_input_string = m2_input_string.replace("{{VERBOSE}}", "true")
    #print("Debug 3: m2_input_string = " + m2_input_string)

    # Step 3: Run Macaulay2
    result = macaulay2.eval(m2_input_string)
    return str(result)

