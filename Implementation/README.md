# Inverse Modeling of Elastic Thickness - Implementation

This implementation is based on the method described in:

**Braitenberg, C., Ebbing, J., & Götze, H.-J. (2002).** Inverse modelling of elastic thickness by convolution method – the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(3-4), 387-404.

## Overview

This code implements a convolution-based inverse modeling approach to estimate the effective elastic thickness (Tₑ) of the lithosphere. The method:

1. Uses topography as the load
2. Computes flexural response functions (Green's functions) for different Tₑ values
3. Performs forward modeling via convolution to predict Moho deflection
4. Estimates Tₑ by minimizing the misfit between observed and modeled Moho depths

## Components

### Core Classes

- **GRDReader**: Reads and writes GMT GRD format files
- **FlexureResponse**: Computes flexural response functions (unit impulse responses)
- **Convolution**: Performs 2D convolution operations
- **InverseModel**: Main inverse modeling algorithm

### Key Functions

1. **Forward Modeling**: `InverseModel::forwardModel()`
   - Computes load from topography
   - Generates flexural response function
   - Convolves load with response to get deflection

2. **Inverse Modeling**: `InverseModel::estimateTe()`
   - Grid search over Tₑ values
   - Minimizes misfit between observed and modeled Moho depths

## Building

### Using Visual Studio (Windows)

1. Open `Convolution.slnx` in Visual Studio
2. Add all source files from the Implementation folder to the project
3. Build the project

### Using CMake (Cross-platform)

A CMakeLists.txt file can be created for cross-platform building.

### Using Command Line (g++/clang++)

```bash
g++ -std=c++11 -O2 \
    GRDReader.cpp \
    FlexureResponse.cpp \
    Convolution.cpp \
    InverseModel.cpp \
    main.cpp \
    -o inverse_model
```

## Usage

### Basic Usage

```bash
./inverse_model [topography.grd] [moho_observed.grd] [output_directory/]
```

### Example

```bash
./inverse_model \
    "../Test data_DRK/Single Mount case/1Topo.grd" \
    "../Test data_DRK/Single Mount case/Moho_calc.grd" \
    "output/"
```

## Input Files

- **Topography GRD file**: Elevation data (in meters)
- **Observed Moho GRD file**: Crust-mantle interface depth (in meters)

Both files should be in GMT GRD format (ASCII or binary).

## Output Files

- `Te_estimated.grd`: Estimated effective elastic thickness (km)
- `deflection_modeled.grd`: Modeled Moho deflection (m)
- `misfit.grd`: Spatial misfit map (m²)

## Parameters

Default parameters (Eastern Alps case):

- Tₑ range: 5-50 km
- Tₑ step: 1 km
- Crustal density (ρₐ): 2800 kg/m³
- Mantle density (ρₘ): 3300 kg/m³
- Young's modulus (E): 70 GPa
- Poisson's ratio (ν): 0.25

These can be modified in `InverseModel::getDefaultEasternAlpsParams()` or by creating custom parameter structures.

## Theory

### Flexural Rigidity

\[ D = \frac{E T_e^3}{12(1 - \nu^2)} \]

### Flexural Response

The lithospheric deflection is computed via convolution:

\[ w(x, y) = \iint R(x - x', y - y') \cdot q(x', y') \, dx' \, dy' \]

Where:
- \( w(x, y) \) = deflection
- \( R(x - x', y - y') \) = flexural response function (Green's function)
- \( q(x', y') \) = applied load

### Load Calculation

\[ q = \rho_c \cdot g \cdot h \]

Where:
- \( \rho_c \) = crustal density
- \( g \) = gravitational acceleration
- \( h \) = topography height

## Limitations and Future Improvements

1. **FFT Convolution**: Currently uses spatial convolution, which is slow for large grids. FFT-based convolution should be implemented for better performance.

2. **Spatially Variable Tₑ**: Current implementation estimates a constant Tₑ. Extension to spatially variable Tₑ would require local windowing or regularization.

3. **GRD Format**: The GRD reader assumes ASCII format. Full GMT binary GRD support should be added.

4. **Optimization**: Grid search is used. Gradient-based optimization (e.g., L-BFGS) could be more efficient.

5. **Uncertainty Estimation**: No uncertainty estimates are provided. Bootstrap or Monte Carlo methods could be added.

## References

Braitenberg, C., Ebbing, J., & Götze, H.-J. (2002). Inverse modelling of elastic thickness by convolution method – the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(3-4), 387-404.

## License

[Specify license here]

## Author

Implementation based on Braitenberg et al. (2002) methodology.

