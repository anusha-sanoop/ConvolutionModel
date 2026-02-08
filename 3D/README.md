# 3D Elastic Thickness Inversion

This folder contains code for 3D elastic thickness estimation using the Braitenberg convolution method, extended to work with 3D volumetric data.

## Overview

The 3D implementation extends the 2D method (in `../Code/`) to handle:
- 3D volumetric topography and Moho depth data
- 3D moving window analysis
- Full 3D convolution operations in frequency domain

## Files

- **`data_loader_3d.py`**: Functions to load and prepare 3D volumetric data
  - `read_3d_volume_from_slices()`: Load 3D volume from multiple 2D slices
  - `prepare_3d_data_for_inversion()`: Prepare 3D data (remove mean, apply taper)
  - `check_3d_grid_compatibility()`: Verify coordinate compatibility

- **`elastic_thickness_inversion_3d.py`**: Core 3D inversion class
  - `ElasticThicknessInversion3D`: Main class for 3D elastic thickness inversion
  - Uses 3D FFT for convolution operations
  - Extends Braitenberg method to 3D wavenumber domain

- **`moving_window_analysis_3d.py`**: 3D moving window analysis
  - `MovingWindowAnalysis3D`: Class for 3D spatial analysis
  - Supports 3D windows with configurable depth extent
  - Multiple shift distances in X-Y plane

- **`main_3d.py`**: Main script for 3D inversion workflow
  - Interactive data loading (multiple slices or single file)
  - Global 3D inversion
  - 3D moving window analysis
  - 3D visualization

## Usage

### Running the 3D Inversion

```bash
python main_3d.py
```

### Data Input Options

1. **Multiple 2D Slices**: Provide a file pattern or list of files
   - Example: `topo_slice_*.grd` or list of individual files
   - Files should be sorted by depth/elevation

2. **Single 2D File**: Treat a single 2D grid as one depth slice
   - Useful for testing or when you only have surface data

### Key Differences from 2D Version

- **3D FFT**: Uses `fft.fftn()` and `fft.ifftn()` for 3D frequency domain operations
- **3D Wavenumber**: Calculates 3D wavenumber magnitude `k = sqrt(kx² + ky² + kz²)`
- **Depth Dimension**: Additional depth dimension in all arrays
- **3D Windows**: Moving windows can have configurable depth extent

## Requirements

Same as the 2D version:
- numpy
- scipy
- matplotlib

## Output

The 3D inversion produces:
- `inversion_results_3d.npz`: All results in NumPy format
- `input_data_3d.png`: Visualization of input 3D data
- `te_map_3d_shift_{N}km.png`: Te maps for each shift distance

## Notes

- The 3D method is computationally more intensive than 2D
- Memory requirements scale with volume size (nz × ny × nx)
- For large volumes, consider using smaller window sizes or depth subsets

