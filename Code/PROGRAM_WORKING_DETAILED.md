# Detailed Program Working Explanation

## Overview

This program implements **Elastic Thickness (Te) Inversion** using the **Braitenberg Convolution Method** (Braitenberg et al., 2002). It estimates the effective elastic thickness of the lithosphere from topography and Moho depth data, which is crucial for understanding the mechanical strength of planetary crusts.

**Main Purpose**: Estimate spatial variations in elastic thickness (Te) across a study region by comparing observed Moho depth with flexure-predicted Moho depth from topography.

---

## 1. INPUTS

### 1.1 Required Input Files

The program takes **two input files** in Surfer GRD format (DSAA ASCII grid):

#### **Topography File** (Required)
- **Format**: `.grd` file (Surfer DSAA format)
- **Content**: Topographic elevation data (meters)
- **Example**: `MarsModel3TwoMountainsFar_Topo.grd`

#### **Moho Depth File OR Gravity File** (One Required)
- **Option A - Direct Moho Depth**:
  - **Format**: `.grd` file
  - **Content**: Absolute Moho depth values (meters, typically negative below surface)
  - **Example**: `TwoMountainsFar_2Mohoflexure.grd`

- **Option B - Bouguer Gravity** (Alternative):
  - **Format**: `.grd` file
  - **Content**: Bouguer gravity anomaly (mGal)
  - **Example**: `MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd`
  - **Note**: If gravity is provided, Moho depth is predicted using Parker-Oldenburg inversion

### 1.2 User-Provided Parameters

During execution, the program prompts for:

1. **Reference Moho Depth** (optional):
   - Input in **kilometers**
   - If not provided, uses mean of Moho depth data
   - Used to calculate Moho undulations (deviations from reference)

2. **Moving Window Parameters**:
   - **Window Size**: Size of analysis window (default: 1000 km)
   - **Shift Range**: Minimum and maximum shift distances (default: 20-80 km)
   - **Shift Step**: Step size for shift distances (default: 20 km)
   - **Te Search Range**: Minimum and maximum Te values to search (default: 5-80 km)

### 1.3 Physical Constants (from `constants.py`)

The program uses Mars-specific physical parameters:

- **Densities**:
  - `RHO_LOAD = 2900 kg/m³` (crustal density)
  - `RHO_MANTLE = 3500 kg/m³` (mantle density)
  - `RHO_INFILL = 2900 kg/m³` (infill density for crustal compensation)
  - `DENSITY_CONTRAST = 600 kg/m³` (for gravity inversion)

- **Elastic Parameters**:
  - `YOUNGS_MODULUS = 1.0e11 Pa` (100 GPa)
  - `POISSONS_RATIO = 0.25`

- **Gravitational Parameters**:
  - `GRAVITY = 3.72 m/s²` (Mars surface gravity)
  - `GRAVITATIONAL_CONSTANT = 6.674e-11 m³/kg/s²`

---

## 2. DATA PROCESSING WORKFLOW

### 2.1 Data Loading (`data_loader.py`)

#### Step 1: Read Grid Files
- **Function**: `read_surfer_grd(filepath)`
- **Process**:
  1. Reads DSAA format header
  2. Extracts grid dimensions (nx, ny)
  3. Reads coordinate ranges (xmin, xmax, ymin, ymax)
  4. Reads data values row-by-row
  5. Creates coordinate grids (X, Y) using `meshgrid`
  6. Calculates grid spacing (dx, dy)

#### Step 2: Grid Compatibility Check
- **Function**: `check_grid_compatibility(X1, Y1, X2, Y2)`
- **Purpose**: Ensures topography and Moho grids have:
  - Same dimensions
  - Matching coordinates (within tolerance)

### 2.2 Moho Prediction from Gravity (if applicable)

If gravity data is provided instead of Moho:

- **Module**: `gravity_moho_inversion.py`
- **Method**: **Parker-Oldenburg Frequency Domain Inversion**
- **Process**:
  1. Convert gravity from mGal to m/s²
  2. Apply FFT to convert to wavenumber domain
  3. Iterative inversion:
     - Forward model: Calculate predicted gravity from current Moho estimate
     - Calculate residual: Observed - Predicted gravity
     - Update Moho depth using residual
     - Repeat until convergence (max 10 iterations)
  4. Output: Absolute Moho depth (meters)

**Key Equation**:
```
G(k) = 2πG Δρ × H(k) × exp(-k×d₀)
```
Where:
- `G(k)` = Gravity in wavenumber domain
- `H(k)` = Moho undulation in wavenumber domain
- `Δρ` = Density contrast
- `d₀` = Reference Moho depth

### 2.3 Data Preparation (`data_loader.py`)

#### Step 1: Calculate Anomalies
- **Function**: `prepare_data_for_inversion()`
- **Process**:
  1. **Topography Anomaly**: Remove mean
     ```python
     topo_anom = topography - mean(topography)
     ```
  2. **Moho Undulation**: Remove reference (user-provided or mean)
     ```python
     moho_undulation = moho_depth - reference_moho
     ```
  3. **Purpose**: Work with deviations from reference, not absolute values

#### Step 2: Apply Taper
- **Function**: `apply_taper(data, alpha=0.1)`
- **Purpose**: Reduce edge effects in FFT calculations
- **Method**: 2D cosine taper (Tukey window) applied to 10% of edges
- **Effect**: Smoothly reduces data to zero at boundaries

#### Step 3: Sign Convention Check
- **Purpose**: Verify correct sign convention
- **Expected**: Positive topography → Negative Moho deflection (deeper)
- **If incorrect**: Program prompts user to flip Moho sign

---

## 3. CORE METHODS

### 3.1 Elastic Thickness Inversion (`elastic_thickness_inversion.py`)

#### 3.1.1 Flexure Theory

The program uses **flexural isostasy** theory:

**Key Concept**: The lithosphere behaves as a thin elastic plate that flexes under loads. The amount of flexure depends on:
- Load magnitude (topography)
- Elastic thickness (Te)
- Material properties (E, ν, densities)

#### 3.1.2 Flexure Filter Function

**Function**: `calculate_flexure_filter(k, Te)`

**Mathematical Basis** (Braitenberg et al., 2002):

In wavenumber domain, the flexure equation is:
```
w(k) = h(k) × [ρ_c / (ρ_m - ρ_infill)] × [1 / (1 + Φ(k))]
```

Where:
- `w(k)` = Moho deflection in frequency domain
- `h(k)` = Topographic load in frequency domain
- `Φ(k) = (D × k⁴) / (g × (ρ_m - ρ_infill))` = Flexural parameter
- `D = E × Te³ / (12 × (1 - ν²))` = Flexural rigidity

**Physical Interpretation**:
- **At long wavelengths (k→0)**: `Φ→0`, response approaches **Airy isostasy**
- **At short wavelengths (k→∞)**: `Φ→∞`, plate is rigid, no flexure
- **At intermediate wavelengths**: Flexural response depends on Te

#### 3.1.3 Forward Model: Predict Moho from Topography

**Function**: `predict_moho_flexure(topography_load, Te)`

**Process**:
1. **FFT of Topography**: Convert to wavenumber domain
   ```python
   load_fft = fft.fft2(topography_load)
   ```

2. **Calculate Flexure Filter**: `F(k) = Airy_ratio / (1 + Φ(k))`

3. **Apply Filter**: Multiply in frequency domain
   ```python
   moho_fft = -flexure_filter * load_fft
   ```
   (Negative sign: positive load → downward deflection)

4. **Inverse FFT**: Convert back to spatial domain
   ```python
   moho_pred = np.real(fft.ifft2(moho_fft))
   ```

**Output**: Predicted Moho undulation (zero-mean)

#### 3.1.4 Inversion: Find Best Te

**Function**: `invert_elastic_thickness(topography_load, moho_obs, Te_range)`

**Method**: **Minimize RMS Misfit**

1. **Misfit Function**:
   ```python
   moho_pred = predict_moho_flexure(topography_load, Te)
   residual = moho_obs - moho_pred
   rms = sqrt(mean(residual²))
   ```

2. **Optimization**:
   - **Method**: Bounded optimization (`scipy.optimize.minimize_scalar`)
   - **Search Range**: User-specified (default: 5-80 km)
   - **Goal**: Find Te that minimizes RMS misfit

3. **Output**:
   - `Te_best`: Best-fit elastic thickness (meters)
   - `rms_best`: Minimum RMS misfit (meters)
   - `moho_pred`: Predicted Moho undulation

### 3.2 Moving Window Analysis (`moving_window_analysis.py`)

#### 3.2.1 Purpose

Estimate **spatial variations** in Te across the study area, rather than a single global value.

#### 3.2.2 Process

**Function**: `analyze(topography_anom, moho_undulation, window_size, shift_distance)`

**Steps**:

1. **Define Window Grid**:
   - Calculate window size in pixels: `window_pixels = window_size / dx`
   - Calculate shift in pixels: `shift_pixels = shift_distance / dx`
   - Generate window positions: `x_positions = arange(0, nx-window_pixels, shift_pixels)`

2. **For Each Window**:
   ```python
   for each window position:
       # Extract data
       topo_window = topography_anom[y_start:y_end, x_start:x_end]
       moho_window = moho_undulation[y_start:y_end, x_start:x_end]
       
       # Ensure zero-mean
       topo_window = topo_window - mean(topo_window)
       moho_window = moho_window - mean(moho_window)
       
       # Check data quality
       if std(topo_window) > min_std and std(moho_window) > min_std:
           # Perform inversion
           result = inverter.invert_elastic_thickness(
               topo_window, moho_window, Te_range
           )
           
           # Store results
           Te_map[i, j] = result['Te_best']
           rms_map[i, j] = result['rms_best']
   ```

3. **Output**:
   - `Te_map`: 2D array of Te values at each window center
   - `rms_map`: 2D array of RMS misfit at each window
   - `x_centers`, `y_centers`: Coordinates of window centers

#### 3.2.3 Multiple Shift Distances

**Function**: `analyze_multiple_shifts()`

- Performs moving window analysis for **multiple shift distances**
- **Purpose**: Test sensitivity to window overlap
- **Typical Range**: 20-80 km shift, 20 km steps
- **Output**: Dictionary with results for each shift distance

**Why Multiple Shifts?**:
- Smaller shifts: Higher spatial resolution, more windows
- Larger shifts: Faster computation, less overlap
- Comparison helps assess robustness of results

---

## 4. OUTPUTS

### 4.1 Output Folder Structure

```
Output_YYYYMMDD_HHMMSS/
├── terminal_output.txt          # Complete terminal log
├── input_data.png               # Topography and Moho visualization
├── inversion_results.npz        # All data arrays (NumPy archive)
├── te_map_shift_XXkm.png        # Te map for each shift (2D)
├── rms_map_shift_XXkm.png       # RMS misfit map for each shift (2D)
└── 3D/
    ├── topography_moho_combined_3d.png
    ├── te_map_shift_XXkm_3d.png
    └── rms_map_shift_XXkm_3d.png
```

### 4.2 Output Files

#### 4.2.1 Visualization Files

1. **`input_data.png`**:
   - **Layout**: 1×2 grid (2 panels)
   - **Panel 1**: Topography (m) - color-coded elevation
   - **Panel 2**: Moho depth (m) - color-coded depth
   - **Purpose**: Visualize input data

2. **`te_map_shift_XXkm.png`**:
   - **Content**: Spatial map of elastic thickness (km)
   - **Color Scale**: Jet colormap
   - **One file per shift distance**: e.g., `te_map_shift_20km.png`, `te_map_shift_40km.png`
   - **Purpose**: Show spatial variations in Te

3. **`rms_map_shift_XXkm.png`**:
   - **Content**: RMS misfit map (km)
   - **Color Scale**: RdYlGn_r (green = good fit, red = poor fit)
   - **Purpose**: Assess quality of Te estimates

4. **3D Visualizations**:
   - **`topography_moho_combined_3d.png`**: 3D surface plot of topography and Moho
   - **`te_map_shift_XXkm_3d.png`**: 3D surface plot of Te map
   - **`rms_map_shift_XXkm_3d.png`**: 3D surface plot of RMS map

#### 4.2.2 Data Files

1. **`inversion_results.npz`** (NumPy archive):
   Contains all arrays:
   - `topography`: Original topography (m)
   - `topography_anomaly`: Zero-mean topography (m)
   - `moho_depth`: Original Moho depth (m)
   - `moho_undulation`: Zero-mean Moho (m)
   - `X`, `Y`: Coordinate grids (m)
   - `stats`: Statistics dictionary
   - `Te_map_shift_XXkm`: Te map for each shift (m)
   - `rms_map_shift_XXkm`: RMS map for each shift (m)
   - `x_centers_shift_XXkm`, `y_centers_shift_XXkm`: Window center coordinates

2. **`terminal_output.txt`**:
   - Complete log of all terminal output
   - Includes diagnostic checks, statistics, warnings

### 4.3 Output Statistics

The program prints comprehensive statistics:

- **Data Statistics**:
  - Topography mean, std (original and anomaly)
  - Moho mean, std (original and undulation)
  - Taper status

- **Diagnostic Checks**:
  - Topo-Moho correlation
  - Sign convention verification
  - Compensation ratio (actual vs. Airy)

- **Moving Window Results**:
  - Number of valid windows
  - Te range, mean, median, std
  - Warnings if hitting search bounds

---

## 5. KEY ALGORITHMS

### 5.1 Frequency Domain Convolution

**Why FFT?**
- Convolution in spatial domain: O(N²) operations
- Multiplication in frequency domain: O(N log N) operations
- **Much faster** for large grids

**Process**:
```
Spatial Domain          Frequency Domain
─────────────────       ─────────────────
h(x,y)  ──FFT──>        H(kx,ky)
                        × F(kx,ky)  [Filter]
w(x,y) <──IFFT──        W(kx,ky)
```

### 5.2 Flexural Rigidity Calculation

**Formula**:
```
D = E × Te³ / (12 × (1 - ν²))
```

**Physical Meaning**:
- **D** = Flexural rigidity (resistance to bending)
- **Te³** relationship: Small changes in Te → Large changes in D
- Example: Te doubles → D increases 8× (very sensitive!)

### 5.3 Airy Isostasy Limit

**When Te → 0**:
- `Φ(k) → 0`
- `F(k) → Airy_ratio = ρ_crust / (ρ_mantle - ρ_infill)`
- **Result**: Pure isostatic compensation (no flexural support)

**For Mars** (typical values):
```
Airy_ratio = 2900 / (3500 - 2900) = 2900 / 600 ≈ 4.83
```
Meaning: 1 m of topography → ~4.8 m of Moho deflection

### 5.4 Optimization Strategy

**Bounded Optimization**:
- Uses `scipy.optimize.minimize_scalar` with `method='bounded'`
- **Advantages**:
  - Fast convergence
  - Guaranteed to stay within bounds
  - No need for initial guess

**Alternative: Grid Search**:
- Test discrete Te values
- Useful for sensitivity analysis
- Slower but more robust

---

## 6. PHYSICAL INTERPRETATION

### 6.1 What is Elastic Thickness (Te)?

**Definition**: The thickness of an equivalent elastic plate that would produce the same flexural response as the actual lithosphere.

**Key Points**:
- **NOT** the actual crustal thickness
- **NOT** the mechanical thickness
- **IS** an effective parameter that captures:
  - Temperature-dependent strength
  - Composition variations
  - Geological history

### 6.2 What Te Values Mean

- **Low Te (< 20 km)**:
  - Weak lithosphere
  - High temperature (close to melting)
  - Recent tectonic activity
  - Example: Young volcanic regions

- **High Te (> 50 km)**:
  - Strong lithosphere
  - Low temperature
  - Old, stable regions
  - Example: Ancient cratons

### 6.3 Spatial Variations

**Why Moving Window?**
- Te varies spatially due to:
  - Temperature gradients
  - Compositional changes
  - Geological history
  - Loading history

**Window Size Considerations**:
- **Too small**: Poor statistics, noisy results
- **Too large**: Smooths out real variations
- **Typical**: 3-10× the flexural wavelength

---

## 7. LIMITATIONS AND CONSIDERATIONS

### 7.1 Assumptions

1. **Thin Plate Theory**: Valid for Te << horizontal scale
2. **Elastic Behavior**: No plastic deformation
3. **Uniform Properties**: E, ν, densities constant
4. **2D Approximation**: Assumes plane strain
5. **No Time Dependence**: Steady-state solution

### 7.2 Edge Effects

- **Problem**: FFT assumes periodic boundaries
- **Solution**: Taper applied to edges (10% of domain)
- **Trade-off**: Some data loss at edges

### 7.3 Data Quality Requirements

- **Minimum Window Size**: Must have sufficient data variation
- **Grid Compatibility**: Topography and Moho must align
- **Signal-to-Noise**: Low noise required for reliable inversion

### 7.4 Non-Uniqueness

- **Problem**: Multiple Te values can fit data similarly
- **Mitigation**: 
  - Use RMS maps to assess confidence
  - Compare multiple shift distances
  - Check for bound-hitting

---

## 8. WORKFLOW SUMMARY

```
START
  │
  ├─> Load Topography (.grd)
  │
  ├─> Load Moho OR Gravity (.grd)
  │   └─> [If Gravity] Predict Moho (Parker-Oldenburg)
  │
  ├─> Check Grid Compatibility
  │
  ├─> Prepare Data:
  │   ├─> Calculate Anomalies (remove mean/reference)
  │   ├─> Apply Taper (reduce edge effects)
  │   └─> Check Sign Convention
  │
  ├─> Diagnostic Check:
  │   ├─> Calculate Correlation
  │   ├─> Check Compensation Ratio
  │   └─> Verify Sign Convention
  │
  ├─> Moving Window Analysis:
  │   ├─> For Each Shift Distance:
  │   │   ├─> For Each Window:
  │   │   │   ├─> Extract Data
  │   │   │   ├─> Invert for Te (minimize RMS)
  │   │   │   └─> Store Te and RMS
  │   │   └─> Create Te and RMS Maps
  │   └─> Generate Visualizations
  │
  ├─> Save Results:
  │   ├─> NumPy Archive (.npz)
  │   ├─> Visualization Images (.png)
  │   └─> Terminal Log (.txt)
  │
  └─> END
```

---

## 9. REFERENCES

1. **Braitenberg, C., Ebbing, J., & Götze, H. J. (2002)**. Inverse modelling of elastic thickness by convolution method—the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(2), 387-404.

2. **Parker, R. L. (1973)**. The rapid calculation of potential anomalies. *Geophysical Journal International*, 31(4), 447-455.

3. **Oldenburg, D. W. (1974)**. The inversion and interpretation of gravity anomalies. *Geophysics*, 39(4), 526-536.

---

## 10. USAGE EXAMPLE

```python
# Run the program
python main.py

# Interactive prompts:
# 1. Enter topography file path
# 2. Choose: Moho file or Gravity file
# 3. Enter reference Moho depth (km) [optional]
# 4. Enter window size (km) [default: 1000]
# 5. Enter shift parameters [defaults: 20-80 km, step 20 km]
# 6. Enter Te search range [default: 5-80 km]

# Output: Results saved in timestamped folder
```

---

## END OF DOCUMENTATION

