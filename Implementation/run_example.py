#!/usr/bin/env python
"""
Example script showing how to run the inverse modeling with specific input files
"""

import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from main import main

if __name__ == "__main__":
    # Example file paths - modify these to match your setup
    input_file1 = r"D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd"
    input_file2 = r"D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\TwoMountainsFar_2Mohoflexure.grd"
    output_dir = "output/"
    
    # Set command line arguments
    sys.argv = ["main.py", input_file1, input_file2, output_dir]
    
    # Run main program
    exit_code = main()
    sys.exit(exit_code)






<<<<<<< HEAD


=======
>>>>>>> e164d71a87412ccae8297068cbd15b9f29754aad
