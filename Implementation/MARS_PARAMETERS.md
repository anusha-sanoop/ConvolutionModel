# Mars Parameters Configuration

This document describes the Mars-specific parameters for the inverse modeling analysis.

## Mars Physical Parameters

Based on your specifications:

- **Gravity**: 3.72 m/s² (Mars surface gravity)
- **Crustal Density**: 2900 kg/m³
- **Mantle Density**: 3500 kg/m³

## Sliding Window Parameters

- **Window Size**: 1,000,000 m (1000 km)
- **Shift Minimum**: 20,000 m (20 km)
- **Shift Maximum**: 80,000 m (80 km)
- **Shift Step**: 20,000 m (20 km)

## Elastic Thickness Search Range

- **Te Minimum**: 5 km
- **Te Maximum**: 50 km
- **Te Step**: 1 km

## Usage

The Mars parameters are automatically used when running `main.py`. To customize:

```python
from inverse_model import InverseModel, InverseParams

# Use Mars parameters with sliding window
params = InverseModel.get_mars_params(
    window_size=1000000.0,    # 1000 km
    shift_min=20000.0,        # 20 km
    shift_max=80000.0,        # 80 km
    shift_step=20000.0,       # 20 km
    Te_min=5.0,
    Te_max=50.0,
    Te_step=1.0
)

# Or modify individual parameters
params = InverseModel.get_mars_params()
params.Te_min = 10.0  # Change Te minimum
params.window_size = 2000000.0  # Change window size to 2000 km
```

## Sliding Window Method

The sliding window approach:

1. Extracts a window of specified size around each analysis point
2. Estimates Te for that window using grid search
3. Assigns the estimated Te to the center of the window
4. Shifts the window by the specified step size
5. Repeats until the entire grid is analyzed
6. Interpolates the results to create a spatially variable Te map

This allows for:
- Spatial variations in Te
- Local analysis of lithospheric properties
- Better resolution of structural boundaries

## Comparison with Earth Parameters

**Earth (Eastern Alps) Default:**
- Gravity: 9.81 m/s²
- Crustal Density: 2800 kg/m³
- Mantle Density: 3300 kg/m³

**Mars:**
- Gravity: 3.72 m/s² (lower gravity affects flexural response)
- Crustal Density: 2900 kg/m³ (slightly higher)
- Mantle Density: 3500 kg/m³ (higher)

The lower gravity on Mars means that for the same load, the flexural response will be different, requiring different Te values to match observed Moho depths.






<<<<<<< HEAD


=======
>>>>>>> e164d71a87412ccae8297068cbd15b9f29754aad
