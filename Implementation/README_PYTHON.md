# Inverse Modeling of Elastic Thickness - Python Implementation

This is a Python implementation of the convolution-based inverse modeling method described in:

**Braitenberg, C., Ebbing, J., & Götze, H.-J. (2002).** Inverse modelling of elastic thickness by convolution method – the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(3-4), 387-404.

## Installation

### 1. Install Python Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- `numpy` (>=1.19.0) - Numerical operations
- `scipy` (>=1.5.0) - Scientific computing (convolution, special functions)

### 2. Verify Installation

```bash
python -c "import numpy, scipy; print('All dependencies installed successfully!')"
```

## Usage

### Quick Start

#### Interactive Mode (Recommended for first-time users):

```bash
python main.py
```

Then enter the paths when prompted:
```
Enter path to Input File 1 (Topography/Bouguer Gravity GRD file): path/to/input1.grd
Enter path to Input File 2 (Observed Moho Depth/Flexure GRD file): path/to/input2.grd
Enter output directory [output/]: output/
```

#### Command Line Mode:

```bash
python main.py input_file1.grd input_file2.grd output_dir/
```

#### Example with Your Files:

```bash
python main.py "D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd" "D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\TwoMountainsFar_2Mohoflexure.grd" "output/"
```

## Input Files

- **Input File 1**: Topography or Bouguer Gravity data (GMT GRD format, ASCII)
- **Input File 2**: Observed Moho depth/flexure data (GMT GRD format, ASCII)

Both files must:
- Be in GMT GRD ASCII format
- Have the same grid dimensions (nx × ny)
- Use the same coordinate system

## Output Files

The program generates three output files in the specified output directory:

1. **`Te_estimated.grd`**: Estimated effective elastic thickness map (km)
2. **`deflection_modeled.grd`**: Modeled Moho deflection (m)
3. **`misfit.grd`**: Spatial misfit map (m²)

## Program Output

The program will display:

1. **Accepted Input**: Summary of input files with dimensions, spatial extent, and data ranges
2. **Inverse Modeling Progress**: For each Te value tested:
   - Te value (km)
   - Misfit (RMS error)
   - Correlation coefficient
3. **Results Summary**:
   - Best Te value
   - Best misfit
   - Best correlation coefficient

## Example Output

```
======================================================================
                     ACCEPTED INPUT
======================================================================

[INPUT FILE 1]
  Path: input1.grd
  Type: Input Data (Topography/Bouguer Gravity)
  Status: ✓ Found
  Dimensions: 200 columns × 200 rows
  Spatial extent:
    X: -100.000000 to 100.000000 (spacing: 1.000000)
    Y: -100.000000 to 100.000000 (spacing: 1.000000)
  Data range: -5000.000 to 8000.000

[INPUT FILE 2]
  Path: input2.grd
  Type: Observed Moho Depth/Flexure
  Status: ✓ Found
  Dimensions: 200 columns × 200 rows
  Spatial extent:
    X: -100.000000 to 100.000000 (spacing: 1.000000)
    Y: -100.000000 to 100.000000 (spacing: 1.000000)
  Data range: 30000.000 to 45000.000

Starting inverse modeling...
Testing Te values from 5.0 to 50.0 km (step: 1.0 km)...
  Te =    5.0 km: Misfit =   1234.567, Correlation =  0.850
  Te =    6.0 km: Misfit =    987.654, Correlation =  0.870
  ...
  Te =   25.0 km: Misfit =    456.789, Correlation =  0.920

Best Te = 25.00 km
Best misfit (RMS) = 456.789
Best correlation = 0.920
```

## Customizing Parameters

Edit `main.py` and modify the parameters:

```python
params = InverseParams(
    Te_min=5.0,        # Minimum Te (km)
    Te_max=50.0,       # Maximum Te (km)
    Te_step=1.0,       # Step size (km)
    rho_m=3300.0,      # Mantle density (kg/m³)
    rho_c=2800.0,      # Crustal density (kg/m³)
    E=70e9,            # Young's modulus (Pa)
    nu=0.25,           # Poisson's ratio
)
```

Or modify the `get_default_eastern_alps_params()` method in `inverse_model.py`.

## Python Module Structure

- **`grd_reader.py`**: GRD file reading/writing (`GridData`, `GRDReader`)
- **`flexure_response.py`**: Flexural response functions (`FlexureResponse`)
- **`convolution.py`**: Convolution operations (`Convolution`)
- **`inverse_model.py`**: Inverse modeling algorithm (`InverseModel`, `InverseParams`, `InverseResult`)
- **`main.py`**: Main program entry point

## Using as a Python Module

You can also import and use the modules in your own Python scripts:

```python
from grd_reader import GRDReader
from inverse_model import InverseModel, InverseParams

# Read data
topography = GRDReader.read_grd("topography.grd")
moho = GRDReader.read_grd("moho.grd")

# Set up parameters
params = InverseParams(Te_min=5.0, Te_max=50.0, Te_step=1.0)

# Run inverse modeling
modeler = InverseModel()
result = modeler.estimate_te(topography, moho, params)

print(f"Best Te: {result.best_Te} km")
```

## Performance Notes

- Uses FFT-based convolution by default (much faster than spatial convolution)
- For grids > 1000×1000, consider increasing Te_step for faster computation
- Processing time scales with: (number of Te values) × (grid size)

## Troubleshooting

### "No module named 'numpy'"
Install dependencies: `pip install -r requirements.txt`

### "Grid dimensions don't match"
Ensure both input GRD files have the same nx and ny values.

### "Cannot open file"
- Check file paths are correct
- Use quotes around paths with spaces
- Verify files exist and are readable

### Poor Results
- Adjust Te search range
- Verify data units (topography in meters, Moho in meters)
- Check density parameters are appropriate for your region

## Differences from C++ Version

- **Faster**: Uses NumPy and SciPy optimized operations
- **FFT Convolution**: Defaults to FFT-based convolution (much faster)
- **Better Error Handling**: More informative error messages
- **Easier to Modify**: Python code is more readable and extensible

## License

[Specify license here]

## References

Braitenberg, C., Ebbing, J., & Götze, H.-J. (2002). Inverse modelling of elastic thickness by convolution method – the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(3-4), 387-404.






<<<<<<< HEAD


=======
>>>>>>> e164d71a87412ccae8297068cbd15b9f29754aad
