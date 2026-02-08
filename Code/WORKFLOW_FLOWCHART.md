# Elastic Thickness Inversion - Workflow & Flowchart
## Braitenberg Convolution Method with Gravity-Based Moho Option

---

## Physical Constants (Mars Parameters)

| Parameter | Symbol | Value | Units |
|-----------|--------|-------|-------|
| Crustal density | ρ_load | 2900 | kg/m³ |
| Mantle density | ρ_m | 3500 | kg/m³ |
| Infill density | ρ_infill | 0 (air) or 2900 (subsurface) | kg/m³ |
| Density contrast | Δρ = ρ_m - ρ_infill | 600 (if air) or 600 (if subsurface) | kg/m³ |
| Gravity | g | 3.72 | m/s² |
| Young's modulus | E | 1.0×10¹¹ | Pa |
| Poisson's ratio | ν | 0.25 | - |
| Gravitational constant | G | 6.674×10⁻¹¹ | m³/kg/s² |

**Note on Infill Density:**
- **ρ_infill = 0**: For topography above sea level (in air)
- **ρ_infill = 2900**: For subsurface loads or crustal compensation
- **ρ_infill = 1000**: For submarine topography (water density)
- This is an **extension** beyond the original Braitenberg (2002) paper, which used (ρ_m - ρ_crust)

---

## Key Equations (As Implemented in Code)

### 1. Flexural Rigidity Factor
\[
D_{factor} = \frac{E}{12(1 - \nu^2)}
\]

### 2. Flexural Rigidity
\[
D = D_{factor} \times T_e^3 = \frac{E \cdot T_e^3}{12(1 - \nu^2)}
\]

**Code implementation:**
```python
D_factor = E / (12 * (1 - nu**2))
D = D_factor * Te**3
```

Where:
- D = Flexural rigidity (N·m)
- E = Young's modulus (Pa) = 1.0×10¹¹ Pa
- T_e = Effective elastic thickness (m)
- ν = Poisson's ratio = 0.25

### 3. Airy Compensation Ratio
\[
Airy_{ratio} = \frac{\rho_{load}}{\rho_m - \rho_{infill}}
\]

**Code implementation:**
```python
airy_ratio = rho_load / (rho_m - rho_infill)
```

Where:
- ρ_load = Load density (kg/m³) = 2900 kg/m³
- ρ_m = Mantle density (kg/m³) = 3500 kg/m³
- ρ_infill = Infill density (kg/m³) = 0 (air) or 2900 (subsurface)

**Example values:**
- For air (ρ_infill = 0): Airy_ratio = 2900 / (3500 - 0) = 0.829
- For subsurface (ρ_infill = 2900): Airy_ratio = 2900 / (3500 - 2900) = 4.833

### 4. Flexural Parameter (Wavenumber Domain)
\[
\Phi(k) = \frac{D \times k^4}{g \times (\rho_m - \rho_{infill})}
\]

**Code implementation:**
```python
Phi_k = (D * k**4) / (g * (rho_m - rho_infill))
```

Where:
- k = Wavenumber magnitude (rad/m) = √(kx² + ky²)
- kx = 2π × fftfreq(nx, dx)
- ky = 2π × fftfreq(ny, dy)
- g = Gravitational acceleration (m/s²) = 3.72 m/s²
- Δρ = ρ_m - ρ_infill = Density contrast (kg/m³)

### 5. Flexure Filter Function
\[
F(k) = \frac{Airy_{ratio}}{1 + \Phi(k)} = \frac{\rho_{load} / (\rho_m - \rho_{infill})}{1 + \frac{D \times k^4}{g \times (\rho_m - \rho_{infill})}}
\]

**Code implementation:**
```python
F_k = airy_ratio / (1.0 + Phi_k)
# At k=0: F(0) = airy_ratio (pure isostatic compensation)
# At high k: F(k) → 0 (plate stiffness dominates)
```

**Physical meaning:**
- At k=0 (long wavelengths): F(0) = Airy_ratio → Full isostatic compensation
- At high k (short wavelengths): F(k) → 0 → Plate is too stiff to flex

### 6. Predicted Moho Undulation (Frequency Domain - Flexure)
\[
W(k) = -F(k) \times H(k)
\]

**Code implementation:**
```python
# FFT of topography
H(k) = FFT2[h(x,y)]

# Apply flexure filter
W(k) = -F(k) * H(k)

# Convert back to spatial domain
w(x,y) = real(IFFT2[W(k)])
```

Where:
- H(k) = Topography load in frequency domain (from FFT)
- W(k) = Predicted Moho undulation in frequency domain
- **Negative sign**: Indicates downward deflection (deeper Moho) under positive loads
- w(x,y) = Predicted Moho undulation in spatial domain (m)

### 7. RMS Misfit Calculation
\[
RMS = \sqrt{\frac{1}{N}\sum_{i=1}^{N} (w_{obs,i} - w_{pred,i})^2}
\]

**Code implementation:**
```python
residual = moho_obs - moho_pred
rms = sqrt(mean(residual**2))
```

Where:
- w_obs = Observed Moho undulation (m)
- w_pred = Predicted Moho undulation (m)
- N = Number of data points
- RMS = Root mean square misfit (m)

### 8. Data Preparation

#### Topography Anomaly
\[
h_{anom}(x,y) = h(x,y) - \bar{h}
\]

#### Moho Undulation
\[
w_{undulation}(x,y) = w_{depth}(x,y) - w_{ref}
\]

**Code implementation:**
```python
topo_mean = mean(topography)
moho_mean = ref_moho (if provided) OR mean(moho_depth)

topo_anom = topography - topo_mean
moho_undulation = moho_depth - moho_mean
```

Where:
- h_anom = Topography anomaly (zero-mean)
- w_undulation = Moho undulation relative to reference
- w_ref = Reference Moho depth (user-provided or mean)

### 9. Taper Function (Edge Effect Reduction)
\[
T(x,y) = T_x(x) \times T_y(y)
\]

Where T_x and T_y are Tukey (cosine taper) windows:

\[
T_x(x) = \begin{cases}
\frac{1}{2}\left[1 + \cos\left(\frac{\pi(x - \alpha L/2)}{\alpha L/2}\right)\right] & \text{if } x < \alpha L/2 \\
1 & \text{if } \alpha L/2 \leq x \leq L(1-\alpha/2) \\
\frac{1}{2}\left[1 + \cos\left(\frac{\pi(L - x - \alpha L/2)}{\alpha L/2}\right)\right] & \text{if } x > L(1-\alpha/2)
\end{cases}
\]

**Code implementation:**
```python
from scipy.signal.windows import tukey
taper_x = tukey(nx, alpha)  # alpha = 0.1 (10% taper)
taper_y = tukey(ny, alpha)
taper_2d = outer(taper_y, taper_x)
data_tapered = data * taper_2d
```

Where:
- α = Taper parameter = 0.1 (10% of domain tapered at edges)
- L = Domain length
- Purpose: Reduces edge effects in FFT processing

### 10. Gravity-to-Moho Inversion (Parker-Oldenburg)

#### Forward Model (Gravity from Moho)
\[
\Delta g(k) = 2\pi G \Delta\rho \times H(k) \times e^{-k d_0}
\]

#### Inverse Model (Moho from Gravity) - Iterative
\[
H_{update}(k) = \frac{\Delta g_{residual}(k)}{2\pi G \Delta\rho \times e^{-k d_0}}
\]

**Code implementation:**
```python
# Convert mGal to m/s²
gravity_si = bouguer_gravity * 1e-5

# Iterative inversion
for iteration in range(max_iterations):
    # Forward model
    moho_fft = FFT2(moho_undulation)
    predicted_gravity_fft = 2*pi*G*density_contrast * moho_fft * exp(-k*d0)
    
    # Residual
    residual_fft = gravity_fft - predicted_gravity_fft
    
    # Update
    update_fft = residual_fft / (2*pi*G*density_contrast * exp(-k*d0))
    moho_depth = moho_depth + real(IFFT2(update_fft))
```

Where:
- Δg(k) = Bouguer gravity anomaly in frequency domain (m/s²)
- G = Gravitational constant = 6.674×10⁻¹¹ m³/kg/s²
- Δρ = Density contrast = 600 kg/m³ (typically)
- d₀ = Reference Moho depth (m)
- H(k) = Moho undulation in frequency domain
- Conversion: 1 mGal = 1×10⁻⁵ m/s²

---

## About Infill Density (ρ_infill)

### Original Braitenberg (2002) Formulation:

The original paper uses:
\[
w(k) = h(k) \times \frac{\rho_c}{\rho_m - \rho_c} \times \frac{1}{1 + \Phi(k)}
\]

Where:
- ρ_c = Crustal density
- ρ_m = Mantle density
- Density contrast = ρ_m - ρ_c

### Extension in This Code:

The code generalizes this to:
\[
w(k) = h(k) \times \frac{\rho_{load}}{\rho_m - \rho_{infill}} \times \frac{1}{1 + \Phi(k)}
\]

Where:
- **ρ_infill** = Density of material filling the space above/below the load
- This allows handling different scenarios:
  - **ρ_infill = 0**: Topography in air (original case)
  - **ρ_infill = 2900**: Subsurface loads or crustal compensation
  - **ρ_infill = 1000**: Submarine topography (water)

**Note:** The infill density parameter is an **extension** beyond the original Braitenberg paper to handle more general loading scenarios. The original paper focused on surface topography in air (effectively ρ_infill = 0).

---

## Complete Workflow Flowchart

```
┌─────────────────────────────────────────────────────────────────┐
│                    START: ELASTIC THICKNESS INVERSION            │
│                      Braitenberg Method                          │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 1: USER INPUT                                             │
│  ─────────────────────────────────────────────────────────────── │
│  • Topography file path (.grd)                                  │
│  • Do you have Bouguer gravity data? (y/n)                      │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ IF YES (Gravity Available):                                │  │
│  │  • Bouguer gravity file path (.grd)                        │  │
│  │  → Will predict Moho from gravity                          │  │
│  └──────────────────────────────────────────────────────────┘  │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ IF NO (No Gravity):                                         │  │
│  │  • Moho depth file path (.grd)                             │  │
│  │  → Will use provided Moho depth                             │  │
│  └──────────────────────────────────────────────────────────┘  │
│  • Reference Moho depth (m) [optional, default: mean]           │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 2: FILE VALIDATION                                        │
│  ─────────────────────────────────────────────────────────────── │
│  Check if files exist                                           │
│  ┌─────────────────┐                                            │
│  │ Files exist?    │                                            │
│  └─────┬───────┬───┘                                            │
│        │       │                                                │
│      YES      NO                                                │
│        │       │                                                │
│        │       └───► ERROR: Exit                               │
│        │                                                        │
└────────┼────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 3: CREATE OUTPUT FOLDER                                   │
│  ─────────────────────────────────────────────────────────────── │
│  Create timestamped output directory                            │
│  Format: Output_YYYYMMDD_HHMMSS                                 │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 4: LOAD DATA                                              │
│  ─────────────────────────────────────────────────────────────── │
│  • Read topography GRD file → X_topo, Y_topo, topography, dx, dy│
│  ┌──────────────────────────────────────────────────────────┐ │
│  │ IF GRAVITY METHOD:                                          │ │
│  │  • Read Bouguer gravity GRD file → X_grav, Y_grav, gravity │ │
│  │  • Check grid compatibility (topo vs gravity)              │ │
│  └──────────────────────────────────────────────────────────┘ │
│  ┌──────────────────────────────────────────────────────────┐ │
│  │ IF OBSERVED MOHO METHOD:                                  │ │
│  │  • Read Moho depth GRD file → X_moho, Y_moho, moho_depth │ │
│  │  • Check grid compatibility (topo vs moho)                │ │
│  └──────────────────────────────────────────────────────────┘ │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 5: MOHO PREDICTION/LOADING                                │
│  ─────────────────────────────────────────────────────────────── │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ IF GRAVITY METHOD:                                        │   │
│  │  ─────────────────────────────────────────────────────── │   │
│  │  5.1: Initialize Gravity Inverter                         │   │
│  │    • density_contrast = 600 kg/m³                         │   │
│  │    • reference_moho = user input or 50 km                 │   │
│  │                                                           │   │
│  │  5.2: Parker-Oldenburg Inversion                          │   │
│  │    For iteration = 1 to max_iterations:                    │   │
│  │    ┌─────────────────────────────────────────────────┐  │   │
│  │    │ 1. Convert gravity: mGal → m/s² (×1e-5)         │  │   │
│  │    │    gravity_si = bouguer_gravity * 1e-5          │  │   │
│  │    │ 2. Create wavenumber grid:                      │  │   │
│  │    │    kx = 2π × fftfreq(nx, dx)                    │  │   │
│  │    │    ky = 2π × fftfreq(ny, dy)                    │  │   │
│  │    │    k = √(kx² + ky²)                             │  │   │
│  │    │ 3. FFT of gravity:                              │  │   │
│  │    │    G(k) = FFT2(gravity_si)                       │  │   │
│  │    │ 4. Forward model:                                │  │   │
│  │    │    G_pred(k) = 2πG Δρ H(k) exp(-k·d₀)            │  │   │
│  │    │ 5. Calculate residual:                           │  │   │
│  │    │    G_res(k) = G_obs(k) - G_pred(k)               │  │   │
│  │    │ 6. Update Moho:                                   │  │   │
│  │    │    H_update(k) = G_res(k) /                      │  │   │
│  │    │                   (2πG Δρ exp(-k·d₀))            │  │   │
│  │    │    moho_depth = moho_depth + IFFT2(H_update)      │  │   │
│  │    │ 7. Check convergence                              │  │   │
│  │    └─────────────────────────────────────────────────┘  │   │
│  │    → Output: Predicted Moho depth from gravity          │   │
│  └──────────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ IF OBSERVED MOHO METHOD:                                  │   │
│  │  → Use loaded Moho depth directly                         │   │
│  └──────────────────────────────────────────────────────────┘   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 6: PREPARE DATA FOR INVERSION                             │
│  ─────────────────────────────────────────────────────────────── │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 6.1: Calculate Reference Levels                          │   │
│  │  • topo_mean = mean(topography)                          │   │
│  │  • moho_mean = ref_moho (if provided) OR mean(moho_depth)│   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 6.2: Calculate Anomalies/Undulations                      │   │
│  │  • topo_anom = topography - topo_mean                      │   │
│  │  • moho_undulation = moho_depth - moho_mean               │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 6.3: Apply Taper (Edge Effect Reduction)                 │   │
│  │  • Taper function: T(x,y) = Tukey(nx,α) × Tukey(ny,α)   │   │
│  │  • α = 0.1 (10% of domain tapered at edges)             │   │
│  │  • topo_anom = topo_anom × T                              │   │
│  │  • moho_undulation = moho_undulation × T                 │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 7: INITIALIZE INVERTER                                    │
│  ─────────────────────────────────────────────────────────────── │
│  Create ElasticThicknessInversion object with:                  │
│  • Grid spacing: dx, dy                                          │
│  • Densities: ρ_load=2900, ρ_m=3500, ρ_infill=0                 │
│  • Elastic parameters: E=1.0×10¹¹ Pa, ν=0.25                    │
│  • Gravity: g=3.72 m/s²                                          │
│                                                                  │
│  Calculate:                                                      │
│  • D_factor = E / [12(1-ν²)]                                    │
│  • Airy_ratio = ρ_load / (ρ_m - ρ_infill)                      │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 8: DIAGNOSTIC CHECK                                       │
│  ─────────────────────────────────────────────────────────────── │
│  • Calculate topo-Moho correlation                              │
│  • Check sign convention                                        │
│  ┌─────────────────┐                                            │
│  │ Correlation > 0? │ (Positive = wrong sign)                   │
│  └─────┬───────┬───┘                                            │
│        │       │                                                │
│      YES      NO                                                │
│        │       │                                                │
│        │       └───► Continue                                  │
│        │                                                        │
│        ▼                                                        │
│  ┌────────────────────────────────────┐                        │
│  │ Ask: Flip Moho sign? (y/n)        │                        │
│  └─────┬──────────────────────────┬───┘                        │
│        │                          │                            │
│       YES                         NO                           │
│        │                          │                            │
│        │                          └───► Continue              │
│        │                                                        │
│        ▼                                                        │
│  moho_undulation = -moho_undulation                            │
│                                                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 9: MOVING WINDOW ANALYSIS                                 │
│  ─────────────────────────────────────────────────────────────── │
│  (Performed for BOTH gravity and observed Moho methods)         │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 9.1: Get User Parameters                                   │   │
│  │  • Window size (km) [default: 1000 km]                    │   │
│  │  • Shift range: min, max, step (km) [default: 20-80 km]   │   │
│  │  • Te search range: min, max (km) [default: 5-80 km]     │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 9.2: For Each Shift Distance                             │   │
│  │  ──────────────────────────────────────────────────────── │   │
│  │                                                             │   │
│  │  ┌─────────────────────────────────────────────────────┐  │   │
│  │  │ 9.2.1: Calculate Window Positions                     │  │   │
│  │  │  • window_pixels = window_size / dx                   │  │   │
│  │  │  • shift_pixels = shift_distance / dx                 │  │   │
│  │  │  • x_positions = [0, shift, 2*shift, ...]            │  │   │
│  │  │  • y_positions = [0, shift, 2*shift, ...]            │  │   │
│  │  └─────────────────────────────────────────────────────┘  │   │
│  │                                                             │   │
│  │  ┌─────────────────────────────────────────────────────┐  │   │
│  │  │ 9.2.2: For Each Window                                │  │   │
│  │  │  ────────────────────────────────────────────────── │  │   │
│  │  │                                                       │  │   │
│  │  │  ┌───────────────────────────────────────────────┐  │  │   │
│  │  │  │ Extract window data:                           │  │   │
│  │  │  │  • topo_window = topo_anom[y:y+win, x:x+win]   │  │   │
│  │  │  │  • moho_window = moho_undulation[y:y+win, ...] │  │   │
│  │  │  │  • Demean each window                           │  │   │
│  │  │  └───────────────────────────────────────────────┘  │  │   │
│  │  │                                                       │  │   │
│  │  │  ┌───────────────────────────────────────────────┐  │  │   │
│  │  │  │ Check data quality:                            │  │   │
│  │  │  │  • std(topo_window) > min_std_topo?            │  │   │
│  │  │  │  • std(moho_window) > min_std_moho?           │  │   │
│  │  │  └─────┬───────────────────────────────────────┬─┘  │  │   │
│  │  │        │                                       │     │  │   │
│  │  │      YES                                      NO     │  │   │
│  │  │        │                                       │     │  │   │
│  │  │        │                                       └───► │  │   │
│  │  │        │                                           Skip│  │   │
│  │  │        │                                           window│  │   │
│  │  │        │                                               │  │   │
│  │  │        ▼                                               │  │   │
│  │  │  ┌───────────────────────────────────────────────┐    │  │   │
│  │  │  │ 9.2.3: Invert for Te                          │    │  │   │
│  │  │  │  ──────────────────────────────────────────── │    │  │   │
│  │  │  │                                                │    │  │   │
│  │  │  │  For each Te in search range:                 │    │  │   │
│  │  │  │  ┌─────────────────────────────────────────┐  │    │  │   │
│  │  │  │  │ 1. Calculate D:                          │  │    │  │   │
│  │  │  │  │    D_factor = E / [12(1-ν²)]            │  │    │  │   │
│  │  │  │  │    D = D_factor × Te³                    │  │    │  │   │
│  │  │  │  │ 2. Create wavenumber grid:               │  │    │  │   │
│  │  │  │  │    kx = 2π × fftfreq(nx, dx)             │  │    │  │   │
│  │  │  │  │    ky = 2π × fftfreq(ny, dy)             │  │    │  │   │
│  │  │  │  │    KX, KY = meshgrid(kx, ky)             │  │    │  │   │
│  │  │  │  │    k = √(KX² + KY²)                       │  │    │  │   │
│  │  │  │  │ 3. Calculate Φ(k):                       │  │    │  │   │
│  │  │  │  │    Phi_k = (D × k⁴) / [g × (ρ_m - ρ_infill)]│  │    │  │   │
│  │  │  │  │ 4. Calculate F(k):                        │  │    │  │   │
│  │  │  │  │    F_k = Airy_ratio / (1 + Phi_k)        │  │    │  │   │
│  │  │  │  │ 5. FFT of topography:                    │  │    │  │   │
│  │  │  │  │    H(k) = FFT2(topo_window)               │  │    │  │   │
│  │  │  │  │ 6. Predicted Moho:                       │  │    │  │   │
│  │  │  │  │    W(k) = -F(k) × H(k)                   │  │    │  │   │
│  │  │  │  │ 7. IFFT:                                 │  │    │  │   │
│  │  │  │  │    w_pred = real(IFFT2[W(k)])            │  │    │  │   │
│  │  │  │  │ 8. Calculate RMS misfit:                  │  │    │  │   │
│  │  │  │  │    residual = moho_window - w_pred       │  │    │  │   │
│  │  │  │  │    RMS = √[mean(residual²)]              │  │    │  │   │
│  │  │  │  └─────────────────────────────────────────┘  │    │  │   │
│  │  │  │                                                │    │  │   │
│  │  │  │  Find Te that minimizes RMS                   │    │  │   │
│  │  │  │  (using scipy.optimize.minimize_scalar)        │    │  │   │
│  │  │  │                                                │    │  │   │
│  │  │  │  Store:                                        │    │  │   │
│  │  │  │  • Te_map[i,j] = Te_best                      │    │  │   │
│  │  │  │  • rms_map[i,j] = RMS_best                    │    │  │   │
│  │  │  │  • x_centers[j] = window center x             │    │  │   │
│  │  │  │  • y_centers[i] = window center y             │    │  │   │
│  │  │  └───────────────────────────────────────────────┘    │  │   │
│  │  │                                                       │  │   │
│  │  └─────────────────────────────────────────────────────┘  │   │
│  │                                                             │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ 9.3: Create Te Maps for Each Shift                       │   │
│  │  • Interpolate Te_map to full grid                        │   │
│  │  • Create 2D and 3D visualization plots                   │   │
│  │  • Save: te_map_shift_XXkm.png                            │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 10: CALCULATE RESIDUAL MOHO                               │
│  ─────────────────────────────────────────────────────────────── │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ IF GRAVITY METHOD:                                        │   │
│  │  • Calculate flexure-predicted Moho using mean Te        │   │
│  │    (from moving window analysis)                          │   │
│  │  • residual_moho = gravity_moho - flexure_moho            │   │
│  └──────────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ IF OBSERVED MOHO METHOD:                                  │   │
│  │  • Calculate flexure-predicted Moho using mean Te        │   │
│  │    (from moving window analysis)                          │   │
│  │  • residual_moho = observed_moho - flexure_moho         │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  • Mask edge regions (10% of domain) to reduce edge artifacts  │
│  • residual_moho[edges] = NaN                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 11: PLOT FINAL RESULTS                                    │
│  ─────────────────────────────────────────────────────────────── │
│  Create 2×2 figure (3 panels, 4th hidden):                     │
│  ┌──────────────┬──────────────┐                               │
│  │ Topography   │ Moho Depth   │                               │
│  │   (m)        │    (m)       │                               │
│  │              │ [from gravity│                               │
│  │              │  if gravity] │                               │
│  ├──────────────┼──────────────┤                               │
│  │ Residual     │ (empty)      │                               │
│  │ Moho (m)     │              │                               │
│  └──────────────┴──────────────┘                               │
│  Save as: input_data.png                                       │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 12: SAVE RESULTS                                          │
│  ─────────────────────────────────────────────────────────────── │
│  Save to .npz file:                                            │
│  • topography, topography_anomaly                              │
│  • moho_depth, moho_undulation                                  │
│  • X, Y coordinates                                            │
│  • Te_map_shift_XXkm (from moving window)                      │
│  • rms_map_shift_XXkm (from moving window)                     │
│  • x_centers, y_centers (from moving window)                   │
│  • Statistics                                                  │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                         END                                     │
│  Output folder contains all results and visualizations         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Two Main Workflows

### Workflow A: Gravity-Based Moho Prediction

```
Topography + Bouguer Gravity
    ↓
[Parker-Oldenburg Inversion]
    ↓
Predicted Moho Depth
    ↓
[Data Preparation]
    ↓
[Moving Window Analysis]
    ↓
Te Map (Elastic Thickness)
    ↓
[Calculate Residual Moho]
    ↓
[Plot 3 Panels: Topo, Moho, Residual]
```

### Workflow B: Observed Moho + Flexure Inversion

```
Topography + Observed Moho
    ↓
[Data Preparation]
    ↓
[Moving Window Analysis]
    ↓
Te Map (Elastic Thickness)
    ↓
[Calculate Residual Moho]
    ↓
[Plot 3 Panels: Topo, Moho, Residual]
```

---

## Forward Model Details

### Flexure Model (Frequency Domain):

1. **Input**: Topography anomaly h(x,y) (zero-mean, tapered)
2. **FFT**: H(k) = FFT2[h(x,y)]
3. **Calculate Flexure Filter**:
   - D_factor = E / [12(1-ν²)]
   - D = D_factor × Te³
   - Φ(k) = (D × k⁴) / [g × (ρ_m - ρ_infill)]
   - F(k) = Airy_ratio / (1 + Φ(k))
4. **Apply Filter**: W(k) = -F(k) × H(k)
5. **IFFT**: w_pred(x,y) = real[IFFT2(W(k))]
6. **Output**: Predicted Moho undulation w_pred(x,y)

### Gravity Inversion Model (Parker-Oldenburg):

1. **Input**: Bouguer gravity anomaly Δg(x,y) (mGal)
2. **Convert**: Δg_si = Δg × 1e-5 (m/s²)
3. **FFT**: G(k) = FFT2[Δg_si]
4. **Iterative Inversion**:
   - Forward: G_pred(k) = 2πG Δρ H(k) exp(-k·d₀)
   - Residual: G_res(k) = G_obs(k) - G_pred(k)
   - Update: H_update(k) = G_res(k) / [2πG Δρ exp(-k·d₀)]
   - Moho: h_new = h_old + real[IFFT2(H_update(k))]
5. **Output**: Predicted Moho depth h(x,y)

---

## Output Files

1. **input_data.png**: 2×2 figure with Topography, Moho Depth, Residual Moho (3 panels)
2. **te_map_shift_XXkm.png**: Te map for each shift distance
3. **te_map_shift_XXkm_3d.png**: 3D Te visualization
4. **inversion_results.npz**: All data arrays and results
5. **3D/**: Folder with 3D visualizations

---

## Key Algorithm References

### Flexure Method:
Braitenberg, C., Ebbing, J., & Götze, H. J. (2002).  
"Inverse modelling of elastic thickness by convolution method—the eastern Alps as a case example."  
Earth and Planetary Science Letters, 202(2), 387-404.

**Note on Infill Density:** The original Braitenberg paper uses (ρ_m - ρ_crust) for density contrast. The infill density parameter (ρ_infill) is an **extension** in this implementation to handle different loading scenarios (air, water, subsurface loads).

### Gravity Method:
Parker, R. L. (1973). The rapid calculation of potential anomalies.  
Geophysical Journal International, 31(4), 447-455.

Oldenburg, D. W. (1974). The inversion and interpretation of gravity anomalies.  
Geophysics, 39(4), 526-536.
