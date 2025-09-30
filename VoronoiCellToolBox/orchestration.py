import subprocess
from importlib.resources import files
from sage.interfaces.macaulay2 import macaulay2
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix

def load_m2_template(inputString, verbose=False):
    """
    Load the templatecomputation.m2 file from your package
    """
    try:
        template_content = files("VoronoiCellToolBox") \
            .joinpath("templatecomputation.m2") \
            .read_text()
        if inputString is None:
            raise ValueError("inputString is None. FormatPullingTrigMatrix(Q) may have returned None.")
        tc = template_content.replace("{{SAGESTRING}};", inputString.replace("\n", "") + ";")
        if( verbose ):
            tc = tc.replace( "{{VERBOSE}}", "true")
        else:
            tc = tc.replace( "{{VERBOSE}}", "false")
        return tc
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
    
    # Step 2: Prepare Macaulay2 input file
    m2_input_string = load_m2_template(sage_string, verbose)
    # print("Debug 2: m2_input_string = " + m2_input_string)

    # Step 3: Run Macaulay2
    result = macaulay2.eval(m2_input_string)
    return str(result)

