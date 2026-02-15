# Code and Paper Guide: Braitenberg Convolution Method for Elastic Thickness

This document explains **Braitenberg et al. (2002)** and how the code implements it: each main function, the equations used, and how they relate to the paper. Use it to understand both the paper and the code in detail.

---

## Reference

**Braitenberg, C., Ebbing, J., & Götze, H. J. (2002).**  
*Inverse modelling of elastic thickness by convolution method—the eastern Alps as a case example.*  
Earth and Planetary Science Letters, 202(2), 387–404.

---

## Part 1: Braitenberg’s Concept and Paper Summary

### 1.1 Physical idea

- The **lithosphere** is treated as a **thin elastic plate** over a fluid mantle. Its strength is summarized by the **effective elastic thickness** \(T_e\).
- **Topography** (e.g. mountains) acts as a **load**. The plate **flexes** in response, so the **Moho** (crust–mantle boundary) is deflected: downward under loads (roots), upward under unloads.
- **Braitenberg’s goal:** Given **observed** topography and **observed** (or gravity-derived) Moho depth, find the \(T_e\) that makes the **flexure-predicted** Moho best match the observed Moho.

So the method is **inverse**: we do not assume \(T_e\); we **infer** it from the relationship between load (topography) and response (Moho).

### 1.2 Convolution in the paper

In the **spatial domain**, the paper writes the predicted Moho deflection as a **convolution** of the load with a **flexural Green’s function** \(G(\mathbf{r}; T_e)\):

\[
w_{\mathrm{pred}}(\mathbf{r}) = -\iint h'(\mathbf{r}')\, G(\mathbf{r} - \mathbf{r}'; T_e)\, d\mathbf{r}'
\quad \equiv \quad
w_{\mathrm{pred}} = -\bigl( G * h' \bigr).
\]

- \(h'\) = topographic load anomaly (zero-mean).
- \(w_{\mathrm{pred}}\) = predicted Moho **undulation** (relative to a reference depth).
- The **minus sign**: positive load (e.g. mountain) → downward deflection (deeper Moho).
- \(G\) depends on \(T_e\) through the **flexural rigidity** \(D\).

So in the paper: “convolution method” = **convolve load with a flexural kernel** to get predicted Moho; then vary \(T_e\) until that prediction fits the observed Moho.

### 1.3 Why the code uses the frequency domain

**Convolution theorem:** Convolution in space \(\Leftrightarrow\) multiplication in the Fourier (wavenumber) domain:

\[
w_{\mathrm{pred}} = -(G * h')
\quad \Leftrightarrow \quad
W(\mathbf{k}) = -F(k)\, H(\mathbf{k}).
\]

- \(H(\mathbf{k})\) = Fourier transform of \(h'\).
- \(F(k)\) = Fourier transform of \(G\); it is the **flexure transfer function** (depends on \(k = |\mathbf{k}|\) and \(T_e\)).
- \(W(\mathbf{k})\) = Fourier transform of \(w_{\mathrm{pred}}\).

So the **same physical model** as the paper is implemented as: **multiply** \(H\) by \(F(k)\) in the wavenumber domain, then inverse FFT to get \(w_{\mathrm{pred}}\). The code does **not** form the spatial kernel \(G\) explicitly; it uses \(F(k)\) (equivalent and efficient with FFT).

### 1.4 Flexure transfer function (paper Eq. 3 and code)

For a thin elastic plate with flexural rigidity \(D\):

**Flexural parameter:**

\[
\Phi(k) = \frac{D\, k^4}{g\,(\rho_m - \rho_{\mathrm{infill}})}.
\]

**Transfer function (ratio of Moho to load in wavenumber domain):**

\[
F(k) = \frac{\rho_{\mathrm{load}}}{\rho_m - \rho_{\mathrm{infill}}} \cdot \frac{1}{1 + \Phi(k)}
= \frac{\mathrm{Airy\_ratio}}{1 + \Phi(k)}.
\]

- **Airy ratio** = \(\rho_{\mathrm{load}}/(\rho_m - \rho_{\mathrm{infill}})\). At \(k=0\) (infinite wavelength), \(F(0)\) equals this ratio → **local (Airy) isostasy**.
- As \(k\) increases, \(\Phi\) grows (\(k^4\)), so \(F(k) \to 0\) → short wavelengths are **supported by the plate** (flexure), not fully compensated.
- \(D\) depends on \(T_e\): \(D = E\,T_e^3 / [12(1-\nu^2)]\). So **stiffer plate (larger \(T_e\)) → larger \(D\) → \(\Phi\) larger at a given \(k\) → \(F\) smaller** → less Moho deflection per unit topography.

**Paper:** Uses \(\rho_{\mathrm{crust}}\) and \(\rho_m\); density contrast \(\rho_m - \rho_{\mathrm{crust}}\).  
**Code:** Generalizes to **infill density** \(\rho_{\mathrm{infill}}\) (e.g. 0 for air, 2900 for crust, 1000 for water), so the contrast is \(\rho_m - \rho_{\mathrm{infill}}\). Same form of \(F(k)\).

### 1.5 Inversion (finding \(T_e\))

- For a given \(T_e\), the code computes \(w_{\mathrm{pred}}\) from \(h'\) using \(F(k)\) (see below).
- **Misfit:** RMS difference between observed Moho undulation \(w'_{\mathrm{obs}}\) and \(w_{\mathrm{pred}}\):

\[
\mathrm{RMS}(T_e) = \sqrt{\frac{1}{N}\sum_{\mathbf{r}} \bigl( w'_{\mathrm{obs}}(\mathbf{r}) - w_{\mathrm{pred}}(\mathbf{r}; T_e) \bigr)^2}.
\]

- **Inversion:** Find \(T_e^*\) that minimizes \(\mathrm{RMS}(T_e)\) over a user-defined range \([T_{e,\min}, T_{e,\max}]\). The code uses a 1D bounded optimizer (e.g. `minimize_scalar`) or a grid search.

---

## Part 2: Constants and Conventions

**File:** `constants.py`

All physical constants used across the code are defined here (Mars defaults).

| Symbol / name        | Typical value | Unit   | Meaning                         |
|----------------------|---------------|--------|---------------------------------|
| \(\rho_{\mathrm{load}}\)  | 2900          | kg/m³  | Crustal (load) density          |
| \(\rho_m\)            | 3500          | kg/m³  | Mantle density                  |
| \(\rho_{\mathrm{infill}}\) | 0 or 2900   | kg/m³  | Infill (air / crust / water)    |
| \(E\)                 | \(10^{11}\)   | Pa     | Young’s modulus                |
| \(\nu\)               | 0.25          | –      | Poisson’s ratio                 |
| \(g\)                 | 3.72          | m/s²   | Surface gravity (Mars)         |
| \(G\)                 | \(6.674\times10^{-11}\) | m³/(kg·s²) | Gravitational constant |
| DEFAULT_REFERENCE_MOHO_DEPTH | 50000 | m      | Fallback reference Moho (50 km) |
| DENSITY_CONTRAST      | \(\rho_m - \rho_{\mathrm{load}}\) | kg/m³ | For gravity inversion |
| DEFAULT_DX, DEFAULT_DY | 1000        | m      | Default grid spacing            |

**Flexural rigidity:**

\[
D = \frac{E\,T_e^3}{12(1-\nu^2)}.
\]

**Airy ratio (used in flexure filter):**

\[
\mathrm{Airy\_ratio} = \frac{\rho_{\mathrm{load}}}{\rho_m - \rho_{\mathrm{infill}}}.
\]

---

## Part 3: Data Loader (`data_loader.py`)

### 3.1 `read_surfer_grd(filepath)`

**Purpose:** Read a Surfer ASCII grid (DSAA format).

**Format:** Header `DSAA`, then lines: `nx ny`, `xmin xmax`, `ymin ymax`, `zmin zmax`, then data rows. Coordinates and values in file are typically in **meters** (topography, Moho); gravity in **mGal**.

**Returns:**  
`X`, `Y` (2D coordinate grids, m), `data` (2D array), `dx`, `dy` (m), `nx`, `ny`, `xmin`, `xmax`, `ymin`, `ymax`.

**Code:** Parses header, builds 1D x/y with `np.linspace`, meshes with `np.meshgrid`, reshapes data to `(ny, nx)` (row = constant y).

---

### 3.2 `write_surfer_grd(filepath, X, Y, data, blank_nan=True)`

**Purpose:** Write a 2D grid to Surfer DSAA. Supports 1D or 2D X/Y; NaN can be written as Surfer blank value.

**Code:** Infers `nx`, `ny` from X/Y; writes header and one row of data per line.

---

### 3.3 `check_grid_compatibility(X1, Y1, X2, Y2)`

**Purpose:** Check that two grids (e.g. topography and gravity, or topography and Moho) have the same shape and coordinates (within tolerance).

**Returns:** `(compatible: bool, message: str)`.

---

### 3.4 `apply_taper(data, alpha=0.1)`

**Purpose:** Apply a 2D **Tukey (cosine) taper** to reduce FFT edge effects.

**Equation:** \(T(x,y) = T_x(x)\,T_y(y)\), with Tukey windows in x and y. `alpha` is the fraction of the domain tapered at each edge (e.g. 0.1 = 10%).

**Code:** `scipy.signal.windows.tukey` for 1D; outer product gives 2D taper; returns `data * taper_2d`.

---

### 3.5 `prepare_data_for_inversion(topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1, ref_moho=None)`

**Purpose:** Convert raw topography and Moho depth into **anomaly** and **undulation** suitable for the convolution/inversion (zero-mean; optional taper).

**Equations:**

- **Topography anomaly (zero-mean load):**
  \[
  h'(\mathbf{r}) = h(\mathbf{r}) - \overline{h}.
  \]

- **Moho undulation (relative to reference):**
  \[
  w'(\mathbf{r}) = w(\mathbf{r}) - w_{\mathrm{ref}},
  \]
  where \(w_{\mathrm{ref}}\) is user-provided reference depth or \(\overline{w}\).

**Code:** Computes means; subtracts to get `topo_anom` and `moho_undulation`; optionally applies `apply_taper` to both; returns these plus a `stats` dict (means, stds, etc.).

---

## Part 4: Elastic Thickness Inversion (`elastic_thickness_inversion.py`)

This module implements the **Braitenberg convolution method** in the **frequency domain**: same physics as the paper, via \(F(k)\) and FFT.

### 4.1 Class `ElasticThicknessInversion.__init__(dx, dy, rho_load, rho_m, rho_infill, E, nu, g)`

**Purpose:** Set grid spacing and physical parameters; precompute **D_factor** and **Airy ratio** used in every filter evaluation.

**Quantities:**

- \(D_{\mathrm{factor}} = E / [12(1-\nu^2)]\) so that \(D = D_{\mathrm{factor}}\,T_e^3\).
- \(\mathrm{Airy\_ratio} = \rho_{\mathrm{load}} / (\rho_m - \rho_{\mathrm{infill}})\).

---

### 4.2 `calculate_flexure_filter(self, k, Te)`

**Purpose:** Compute the **flexure transfer function** \(F(k)\) for a given wavenumber magnitude \(k\) and elastic thickness \(T_e\). This is the frequency-domain equivalent of the paper’s convolution kernel.

**Equations (Braitenberg Eq. 3, equivalent form):**

\[
D = \frac{E\,T_e^3}{12(1-\nu^2)},
\qquad
\Phi(k) = \frac{D\, k^4}{g\,(\rho_m - \rho_{\mathrm{infill}})},
\qquad
F(k) = \frac{\mathrm{Airy\_ratio}}{1 + \Phi(k)}.
\]

**Code:**  
`D = self.D_factor * Te**3`; `Phi_k = (D * k**4) / (self.g * (self.rho_m - self.rho_infill))`; `F_k = self.airy_ratio / (1.0 + Phi_k)`. Returns `F_k` (same shape as `k`).

**Behaviour:** At \(k=0\), \(F(0) = \mathrm{Airy\_ratio}\) (full compensation); at large \(k\), \(F \to 0\) (rigid plate).

---

### 4.3 `predict_moho_flexure(self, topography_load, Te)`

**Purpose:** **Forward model**: from a topographic load anomaly and \(T_e\), compute the predicted Moho undulation (same as paper’s convolution result, via frequency domain).

**Steps (equivalent to \(w_{\mathrm{pred}} = -(G * h')\)):**

1. Demean load: \(h' = \mathrm{topography\_load} - \mathrm{mean}(\mathrm{topography\_load})\).
2. Wavenumber grid: \(k_x = 2\pi\,\mathrm{fftfreq}(n_x, dx)\), \(k_y = 2\pi\,\mathrm{fftfreq}(n_y, dy)\), \(k = \sqrt{k_x^2 + k_y^2}\).
3. \(H(\mathbf{k}) = \mathrm{FFT2}(h')\).
4. \(F(k) = \mathrm{calculate\_flexure\_filter}(k, T_e)\); set \(F(0,0) = \mathrm{Airy\_ratio}\) for stability.
5. \(W(\mathbf{k}) = -F(k)\,H(\mathbf{k})\) (minus sign = downward deflection under positive load).
6. \(w_{\mathrm{pred}}(\mathbf{r}) = \mathrm{real}(\mathrm{IFFT2}(W))\), then demean.

**Returns:** 2D array `moho_pred` (m), zero-mean. Negative values = deeper Moho.

---

### 4.4 `misfit_function(self, Te, topography_load, moho_obs, mask=None)`

**Purpose:** Single scalar misfit for a given \(T_e\): RMS difference between observed and predicted Moho undulation. Used by the optimizer.

**Equation:**

\[
\mathrm{RMS} = \sqrt{\frac{1}{N}\sum \bigl( w'_{\mathrm{obs}} - w_{\mathrm{pred}}(T_e) \bigr)^2}.
\]

**Code:** Calls `predict_moho_flexure(topography_load, Te)`; computes residual `moho_obs - moho_pred` (optionally masked); returns `np.sqrt(np.mean(residual**2))` in meters.

---

### 4.5 `invert_elastic_thickness(self, topography_load, moho_obs, Te_range, mask=None, method='bounded')`

**Purpose:** Find \(T_e^*\) that minimizes RMS misfit over \([T_{e,\min}, T_{e,\max}]\).

**Equation:**

\[
T_e^* = \arg\min_{T_e \in [T_{e,\min},\, T_{e,\max}]} \mathrm{RMS}(T_e).
\]

**Code:**  
- `method='bounded'`: `scipy.optimize.minimize_scalar(..., bounds=Te_range, method='bounded')` on `misfit_function`.  
- `method='grid_search'`: loop over a linear grid of \(T_e\) values, compute RMS for each, take argmin.

**Returns:** Dict with `Te_best` (m), `rms_best` (m), `moho_pred`, `method`, `Te_range`.

---

### 4.6 `sensitivity_analysis(self, topography_load, moho_obs, Te_range, n_points=50, mask=None)`

**Purpose:** Evaluate RMS for many \(T_e\) values to plot Te–RMS curve and assess uniqueness/confidence.

**Returns:** `(Te_values, rms_values)` arrays.

---

### 4.7 `diagnostic_check(self, topography_anom, moho_undulation)`

**Purpose:** Compare observed compensation to Airy expectation; check topo–Moho correlation (sign convention); give a rough interpretation (weak vs strong lithosphere).

**Quantities:**  
- Expected Moho std if Airy: \(\sigma_{w,\mathrm{Airy}} = \mathrm{Airy\_ratio} \times \sigma_h\).  
- Actual ratio: \(\sigma_w / \sigma_h\).  
- Correlation of topo vs Moho (negative expected: high topo → deep Moho).

---

## Part 5: Gravity–Moho Inversion (`gravity_moho_inversion.py`)

**Purpose:** When the user provides **Bouguer gravity** instead of Moho, estimate **Moho depth** from gravity so it can be used as “observed Moho” in the Te inversion. Based on **Parker–Oldenburg** in the frequency domain.

**References:** Parker (1973), Oldenburg (1974).

### 5.1 Class `GravityMohoInversion.__init__(dx, dy, density_contrast, reference_moho, G)`

**Purpose:** Set grid spacing, crust–mantle density contrast \(\Delta\rho\), reference Moho depth \(d_0\), and \(G\). Used in the forward and inverse steps below.

---

### 5.2 `predict_moho_from_gravity(self, bouguer_gravity, max_iterations=10, tolerance=1.0)`

**Purpose:** Recover Moho depth from Bouguer anomaly by iterative Parker–Oldenburg inversion.

**Forward (gravity from Moho):** In wavenumber domain,

\[
\Delta g(\mathbf{k}) = 2\pi G\,\Delta\rho\; \widetilde{w}(\mathbf{k})\; e^{-k\,d_0},
\]

where \(\widetilde{w}\) is the Moho undulation in the frequency domain, \(d_0\) is reference depth.

**Inverse (iterative):**  
1. Start from Moho = reference (flat).  
2. Forward: compute predicted gravity from current Moho.  
3. Residual: \(\Delta g_{\mathrm{obs}} - \Delta g_{\mathrm{pred}}\).  
4. Update: solve for Moho update from residual (divide by \(2\pi G\,\Delta\rho\, e^{-k d_0}\) in k-domain), add to Moho.  
5. Repeat until max change in Moho < `tolerance` or `max_iterations` reached.

**Code:** Converts mGal to m/s² (\(1\,\mathrm{mGal}=10^{-5}\,\mathrm{m/s}^2\)); FFT of gravity; loop with forward FFT of Moho undulation, residual in k-space, update by division in k-space and IFFT; returns 2D Moho depth (m).

---

### 5.3 `predict_moho_simple(self, bouguer_gravity)`

**Purpose:** Quick linear estimate: \(\Delta g \approx 2\pi G\,\Delta\rho\, h\) ⇒ \(h \approx \Delta g / (2\pi G\,\Delta\rho)\). No iteration; useful for rough checks.

---

## Part 6: Moving Window Analysis (`moving_window_analysis.py`)

**Purpose:** Obtain **spatial maps** of \(T_e\) and RMS by running the Braitenberg inversion in many **windows** over the grid (and optionally for several **shift** distances).

### 6.1 Class `MovingWindowAnalysis.__init__(dx, dy, rho_load, rho_m, rho_infill, g)`

**Purpose:** Store grid and density/gravity parameters for creating the inverter in each window.

---

### 6.2 `analyze(self, topography_anom, moho_undulation, window_size, shift_distance, Te_range, min_std_topo, min_std_moho)`

**Purpose:** One **moving-window run**: fixed window size and shift; for each window position, run Te inversion and store one \(T_e\) and one RMS at the window center.

**Steps:**

1. **Window geometry:**  
   `window_pixels = window_size / dx`, `shift_pixels = shift_distance / dx`.  
   Window positions: `x_positions = 0, shift_pixels, 2*shift_pixels, ...` (and same for y) so the full window stays inside the grid.

2. **Per window:**  
   - Extract `topo_window`, `moho_window`; demean each.  
   - Optionally skip if std(topo) < min_std_topo or std(moho) < min_std_moho.  
   - Call `ElasticThicknessInversion().invert_elastic_thickness(topo_window, moho_window, Te_range, method='bounded')`.  
   - Store `Te_best` and `rms_best` in 2D arrays at index \((i,j)\); store window center coordinates in meters in `x_centers`, `y_centers`.

**Returns:** Dict with `Te_map`, `rms_map` (2D, NaN where inversion failed), `x_centers`, `y_centers`, `window_size`, `shift_distance`, `n_windows`, `n_valid`.

---

### 6.3 `analyze_multiple_shifts(self, ...)`

**Purpose:** Run `analyze` for several shift distances (e.g. 20, 40, 60, 80 km). Each shift gives one Te map and one RMS map.

**Returns:** Dict mapping `shift_dist` → result dict from `analyze(...)`.

---

## Part 7: Main Program Flow (`main.py`)

**Purpose:** Orchestrate inputs, data loading, optional gravity–Moho step, data preparation, optional global inversion, then **single-window** or **moving-window** analysis; create figures and save outputs (all displayed in km; Te/RMS maps trimmed at borders as implemented).

**High-level flow:**

1. **Inputs:** Topography .grd path; choice: gravity .grd (y) or Moho .grd (n); reference Moho depth (km) or Enter for mean.
2. **Load:** `read_surfer_grd` for topography; if gravity, load gravity and optionally run `GravityMohoInversion.predict_moho_from_gravity` to get Moho; else load Moho .grd. `check_grid_compatibility` when needed.
3. **Prepare:** `prepare_data_for_inversion` → topo anomaly, Moho undulation, stats.
4. **Optional global inversion:** One Te over the full grid (currently off by default).
5. **Mode:** User chooses single window (1) or moving window (2).
6. **Single window:** Ask window size and Te range; interactive plot to place window; one Te (or Te map within that region); save .npz, .grd, PNGs.
7. **Moving window:** Ask window size, shift range/step, Te range; `MovingWindowAnalysis.analyze_multiple_shifts`; trim border cells from Te/RMS maps; for each shift, create Te and RMS figures (2D and 3D), save .grd and .npz; RMS color scale from full (untrimmed) map.
8. **Final plots:** `plot_final_results` (topography and Moho in km); save all outputs and `terminal_output.txt`.

**Key main.py functions (summary):**

- `create_output_folder()`: timestamped `Output_YYYYMMDD_HHMMSS`.
- `create_filename_prefix(rho_load, ref_moho_km, window_size_km)`: prefix for filenames.
- `plot_final_results(...)`: 1×2 figure Topography (km), Moho depth (km).
- `plot_input_data_3d(...)`: 3D view of topo and Moho in km.
- Te/RMS map figure helpers (e.g. `create_te_map_figure`, `create_rms_map_figure`) use trimmed results and (for RMS) fixed color scale from full map.

---

## Part 8: Equation Summary Table

| Step / quantity        | Equation / expression |
|------------------------|------------------------|
| Flexural rigidity      | \(D = \dfrac{E\,T_e^3}{12(1-\nu^2)}\) |
| Paper (spatial)        | \(w_{\mathrm{pred}} = -(G * h')\) |
| Flexural parameter     | \(\Phi(k) = \dfrac{D\,k^4}{g\,(\rho_m - \rho_{\mathrm{infill}})}\) |
| Transfer function      | \(F(k) = \dfrac{\rho_{\mathrm{load}}/(\rho_m - \rho_{\mathrm{infill}})}{1 + \Phi(k)}\) |
| Code (frequency)       | \(W(\mathbf{k}) = -F(k)\,H(\mathbf{k})\), \(w_{\mathrm{pred}} = \mathrm{real}(\mathrm{IFFT}(W))\) |
| RMS misfit             | \(\mathrm{RMS} = \sqrt{\dfrac{1}{N}\sum (w'_{\mathrm{obs}} - w_{\mathrm{pred}})^2}\) |
| Inversion              | \(T_e^* = \arg\min_{T_e} \mathrm{RMS}(T_e)\) |
| Topo anomaly           | \(h' = h - \overline{h}\) |
| Moho undulation        | \(w' = w - w_{\mathrm{ref}}\) |
| Gravity forward (Parker) | \(\Delta g(\mathbf{k}) = 2\pi G\,\Delta\rho\; \widetilde{w}(\mathbf{k})\; e^{-k\,d_0}\) |

---

## Part 9: Reading the Paper with the Code

- **Paper’s “convolution”** → code’s **frequency-domain multiplication** by \(F(k)\) in `calculate_flexure_filter` and `predict_moho_flexure`.
- **Paper’s “flexural response” / Green’s function** → code’s **transfer function** \(F(k)\) (Eq. 3 in the paper).
- **Paper’s inversion for \(T_e\)** → code’s `misfit_function` + `invert_elastic_thickness` (minimize RMS over \(T_e\)).
- **Spatial moving windows** in the paper → `MovingWindowAnalysis.analyze`: same inversion per window, results on a grid of centers.
- **Density contrast** in the paper (\(\rho_m - \rho_c\)) → code’s generalization with \(\rho_{\mathrm{infill}}\) in \(F(k)\) and in constants.

Using this guide, you can go through the paper equation by equation and see exactly where and how each part is implemented in the code.
