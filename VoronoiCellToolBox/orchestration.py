import subprocess
import pkg_resources
from sage.interfaces.macaulay2 import macaulay2
from VoronoiCellToolBox.macaulay_parsing import FormatPullingTrigMatrix

def load_m2_template(inputString):
    """
    Load the templatecomputation.m2 file from your package
    """
    print("Input String:", inputString)  # Debugging line
    try:
        template_content = pkg_resources.resource_string(
            'VoronoiCellToolBox', 'templatecomputation.m2'
        ).decode('utf-8')
        if inputString is None:
            raise ValueError("inputString is None. FormatPullingTrigMatrix(Q) may have returned None.")
        return template_content.replace("{{SAGESTRING}}", inputString)
    except Exception as e:
        print(f"Error loading resource: {e}")
        raise FileNotFoundError("Could not find templatecomputation.m2 in package")


def sage_to_macaulay2(Q):
    """
    Execute a Sage computation and pass results to Macaulay2.
    
    Parameters:
    -----------
    Q : matrix
        
    Returns:
    --------
    a string
    """
    # Step 1: Run Sage computation
    sage_string = FormatPullingTrigMatrix(Q)
    print(sage_string)
    
    # Step 2: Prepare Macaulay2 input file
    m2_input_string = load_m2_template(sage_string)
    print(m2_input_string)
    
    # Step 3: Run Macaulay2
    return macaulay2.eval(m2_input_string)

