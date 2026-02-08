# Usage Guide

## Quick Start

### 1. Compile the Code

#### Using CMake:
```bash
mkdir build
cd build
cmake ../Implementation
make
```

#### Using Visual Studio:
1. Open `Convolution.slnx` in Visual Studio
2. Add all `.cpp` and `.h` files from the `Implementation` folder to your project
3. Build the solution

#### Using g++ directly:
```bash
cd Implementation
g++ -std=c++11 -O2 *.cpp -o inverse_model
```

### 2. Prepare Input Data

You need two GRD files:
- **Topography file**: Gridded elevation data (meters)
- **Moho file**: Observed crust-mantle interface depth (meters)

Example test data location:
- `../Test data_DRK/Single Mount case/1Topo.grd`
- `../Test data_DRK/Single Mount case/Moho_calc.grd`

### 3. Run the Program

```bash
./inverse_model [topography.grd] [moho_observed.grd] [output_dir/]
```

Example:
```bash
./inverse_model \
    "../Test data_DRK/Single Mount case/1Topo.grd" \
    "../Test data_DRK/Single Mount case/Moho_calc.grd" \
    "output/"
```

### 4. View Results

The program will generate:
- `Te_estimated.grd`: Estimated elastic thickness map (km)
- `deflection_modeled.grd`: Modeled Moho deflection (m)
- `misfit.grd`: Spatial misfit between observed and modeled Moho (m²)

View these files using GMT, Python with matplotlib, or other visualization tools.

## Customizing Parameters

Edit `InverseModel.cpp` and modify `getDefaultEasternAlpsParams()`:

```cpp
InverseParams params;
params.Te_min = 5.0;        // Minimum Te (km)
params.Te_max = 50.0;       // Maximum Te (km)
params.Te_step = 1.0;       // Step size (km)
params.rho_m = 3300.0;      // Mantle density (kg/m³)
params.rho_c = 2800.0;      // Crustal density (kg/m³)
params.E = 70e9;            // Young's modulus (Pa)
params.nu = 0.25;           // Poisson's ratio
```

## Program Flow

1. **Read Input Data**: Loads topography and observed Moho depth grids
2. **Parameter Setup**: Configures inverse modeling parameters
3. **Grid Search**: Tests different Te values
4. **Forward Modeling**: For each Te:
   - Computes load from topography
   - Generates flexural response function
   - Convolves to get modeled deflection
5. **Misfit Calculation**: Compares modeled vs. observed Moho
6. **Best Fit Selection**: Chooses Te with minimum misfit
7. **Output**: Saves results to GRD files

## Expected Output

```
========================================
Inverse Modeling of Elastic Thickness
Convolution Method Implementation
Based on Braitenberg et al. (2002)
========================================

Input files:
  Topography: ../Test data_DRK/Single Mount case/1Topo.grd
  Observed Moho: ../Test data_DRK/Single Mount case/Moho_calc.grd
  Output directory: output/

Reading input data...
Topography grid: 200 x 200
  Range: -5000.0 to 8000.0
Moho grid: 200 x 200
  Range: 30000.0 to 45000.0

Inverse modeling parameters:
  Te range: 5 - 50 km
  Te step: 1 km
  ...

Starting inverse modeling...
Testing Te values from 5 to 50 km...
Te = 5.0 km: Misfit = 1234.56, Correlation = 0.85
Te = 6.0 km: Misfit = 987.65, Correlation = 0.87
...
Best Te = 25.0 km
Best misfit = 456.78
Best correlation = 0.92
```

## Troubleshooting

### Error: "Grid dimensions don't match"
- Ensure topography and Moho grids have the same dimensions
- Check that both files use the same coordinate system

### Error: "Cannot open file"
- Verify file paths are correct
- Check file permissions

### Poor Results
- Adjust Te search range
- Check data quality and units
- Verify density parameters are appropriate for your region

## Next Steps

- Implement spatially variable Te estimation
- Add FFT-based convolution for faster computation
- Include uncertainty estimation
- Add support for multiple loading scenarios

