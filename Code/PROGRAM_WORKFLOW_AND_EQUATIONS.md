# Elastic Thickness Inversion Program — Workflow and Equations

This document describes how the program runs, its workflow, and the main equations used for estimating elastic thickness (Te) from topography and Moho depth.

---

## 1. Overview

The program estimates **effective elastic thickness (Te)** of the lithosphere using the **convolution method** (Braitenberg et al., 2002). It compares observed Moho depth (or Moho predicted from gravity) with the Moho predicted by a thin-plate flexure model driven by topography, and finds the Te that gives the best fit.

**Main reference:**  
Braitenberg, C., Ebbing, J., & Götze, H. J. (2002). Inverse modelling of elastic thickness by convolution method—the eastern Alps as a case example. *Earth and Planetary Science Letters*, 202(2), 387–404.

**Spatial vs frequency domain:** Braitenberg et al. formulate the method in the **spatial domain** (convolution of the load with a flexural Green’s function). This program implements the **same relationship** in the **frequency domain** (convolution theorem: multiplication by the transfer function \(F(k)\) in Fourier space), which is mathematically equivalent and numerically efficient (FFT-based).

**Two analysis modes:**
- **Single window:** One selected region → one Te (or a Te map within that region).
- **Moving window:** Many overlapping windows over the full grid → Te map over the whole domain.

---

## 2. High-Level Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  START                                                                       │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  1. USER INPUTS                                                              │
│     • Topography file (.grd)                                                 │
│     • Moho file (.grd) OR Bouguer gravity file (.grd)                        │
│     • Reference Moho depth (km) or use mean                                  │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  2. LOAD DATA (data_loader.py)                                               │
│     • Read Surfer DSAA .grd → X, Y, topography, dx, dy                        │
│     • If gravity: predict Moho from Bouguer (Parker–Oldenburg)                │
│     • Check grid compatibility                                                │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  3. PREPARE DATA FOR INVERSION                                                │
│     • Topography anomaly: h' = topography − mean(topography)                 │
│     • Moho undulation: w' = moho_depth − reference_moho                       │
│     • Optional: 2D taper (Tukey) to reduce FFT edge effects                   │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  4. CHOOSE ANALYSIS MODE                                                     │
│     • Press 1 → Single window                                                │
│     • Press 2 → Moving window                                                │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ├─── Single window ────────────────────────────────────────────────────┤
        │    • Set window size (km)                                             │
        │    • Set Te range (min, max km)                                       │
        │    • Interactive plot: move cursor, left-click to select region      │
        │    • Option A: One Te for whole region (one inversion)                │
        │    • Option B: Te map within region (moving window inside selection)  │
        │    • Save: .npz, .grd, PNG figures                                    │
        │                                                                       │
        └─── Moving window ────────────────────────────────────────────────────┤
             • Set window size, shift range, Te range (km)                      │
             • Run moving window for several shift distances                    │
             • For each shift: Te map, RMS map → PNG + .grd                     │
             • Save: .npz, .grd, PNG figures                                    │
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  5. FINAL PLOTS & OUTPUTS                                                    │
│     • Topography + Moho summary figure (input_data.png)                       │
│     • All outputs in timestamped folder: Output_YYYYMMDD_HHMMSS/              │
└─────────────────────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  END                                                                         │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Data Loading and Preparation

### 3.1 Surfer grid format (DSAA)

Grids are read in **Surfer ASCII DSAA** format:

- Header: `DSAA`
- Line 2: `nx ny`
- Line 3: `xmin xmax`
- Line 4: `ymin ymax`
- Line 5: `zmin zmax`
- Following lines: data rows (row index = Y from bottom to top in the code convention)

Coordinates and values are in **meters** (topography and Moho depth). Gravity is in **mGal**.

### 3.2 Preparation for inversion

- **Topography anomaly** (zero-mean load):
  \[
  h'(\mathbf{r}) = h(\mathbf{r}) - \overline{h}
  \]
- **Moho undulation** (relative to reference):
  \[
  w'(\mathbf{r}) = w(\mathbf{r}) - w_{\mathrm{ref}}
  \]
  where \(w_{\mathrm{ref}}\) is either the user-given reference depth or the mean Moho depth.

- **Taper:** A 2D Tukey (cosine) taper is applied to \(h'\) and \(w'\) at the edges (default 10% of the domain) to reduce FFT edge effects.

---

## 4. Elastic Thickness Inversion (Convolution Method)

Braitenberg et al. (2002) formulate the method in the **spatial domain** as a convolution of the topographic load with a flexural **response function** (Green’s function). This code implements the **same physical model** in the **frequency domain** (convolution theorem: convolution in space = multiplication in Fourier space), which is numerically efficient and gives identical results.

### 4.1 Thin-plate flexure

The lithosphere is treated as a thin elastic plate. Flexural rigidity is

\[
D = \frac{E\,T_e^3}{12(1-\nu^2)}
\]

where:
- \(E\) = Young’s modulus (Pa)  
- \(\nu\) = Poisson’s ratio  
- \(T_e\) = elastic thickness (m)

### 4.2 Spatial domain (Braitenberg et al., 2002)

In the **spatial domain**, the predicted Moho undulation is the **convolution** of the load with a flexural response function \(G(\mathbf{r}; T_e)\):

\[
w_{\mathrm{pred}}(\mathbf{r}) = -\iint h'(\mathbf{r}')\, G(\mathbf{r} - \mathbf{r}'; T_e)\, d\mathbf{r}'
\quad \equiv \quad
w_{\mathrm{pred}} = -\bigl( G * h' \bigr).
\]

The minus sign gives downward deflection under a positive load. The kernel \(G\) depends on \(T_e\) (through \(D\)) and on the density contrast and gravity; its Fourier transform is the transfer function \(F(k)\) used below. So the paper’s formulation is **convolution in space**; no FFT is required in principle.

### 4.3 Frequency domain (implementation used here)

By the **convolution theorem**, convolution in space is multiplication in the wavenumber domain. So the same relation is implemented as:

\[
W(\mathbf{k}) = -F(k)\,H(\mathbf{k}),
\]

where \(H(\mathbf{k})\) is the 2D Fourier transform of \(h'(\mathbf{r})\), and \(F(k)\) is the **flexure transfer function** (Fourier transform of the spatial response function). For the thin-plate model, \(F(k)\) has the form (e.g. Braitenberg et al., 2002, Eq. 3 in equivalent form):

\[
\Phi(k) = \frac{D\,k^4}{g\,(\rho_m - \rho_{\mathrm{infill}})},
\qquad
F(k) = \frac{\rho_{\mathrm{load}}}{\rho_m - \rho_{\mathrm{infill}}} \cdot \frac{1}{1 + \Phi(k)} = \frac{\mathrm{Airy\_ratio}}{1 + \Phi(k)}.
\]

Here \(k = \sqrt{k_x^2 + k_y^2}\) (rad/m); \(g\), \(\rho_{\mathrm{load}}\), \(\rho_m\), \(\rho_{\mathrm{infill}}\) are as above. At \(k=0\), \(F(0)\) is the Airy ratio; at large \(k\), \(\Phi \gg 1\) so \(F \to 0\).

Predicted Moho in space is then:

\[
w_{\mathrm{pred}}(\mathbf{r}) = \mathrm{IFFT}\bigl[ W(\mathbf{k}) \bigr] \quad \text{(real part, then demean)}.
\]

So: **Braitenberg’s spatial-domain convolution** and **this code’s frequency-domain multiplication** describe the same forward model; the code uses the frequency domain for speed (FFT is \(O(N\log N)\) vs \(O(N^2)\) for a naive spatial convolution).

### 4.4 Misfit and inversion

**RMS misfit** between observed and predicted Moho undulation:

\[
\mathrm{RMS}(T_e) = \sqrt{\frac{1}{N}\sum_{\mathbf{r}} \bigl( w'_{\mathrm{obs}}(\mathbf{r}) - w_{\mathrm{pred}}(\mathbf{r}; T_e) \bigr)^2}
\]

(in meters). The program finds

\[
T_e^* = \arg\min_{T_e \in [T_{e,\min},\, T_{e,\max}]} \mathrm{RMS}(T_e)
\]

using a **bounded 1D optimizer** (e.g. `minimize_scalar`, method `bounded`) over the user-defined range \([T_{e,\min}, T_{e,\max}]\) (in m).

---

## 5. Gravity–Moho Step (Optional)

If the user provides **Bouguer gravity** instead of Moho, Moho depth is first estimated with the **Parker–Oldenburg** frequency-domain method.

- **Forward (gravity from Moho):**  
  \(G(\mathbf{k}) = 2\pi G\,\Delta\rho\, \widetilde{w}(\mathbf{k})\, e^{-k\,d_0}\)  
  where \(\widetilde{w}\) is the Moho undulation in the frequency domain, \(d_0\) is reference depth, \(\Delta\rho\) is crust–mantle density contrast.

- **Inversion:** Iteratively update Moho so that predicted gravity matches observed Bouguer anomaly (see `gravity_moho_inversion.py`). The result is a **Moho depth grid** used as “observed” Moho in the Te inversion.

---

## 6. Single-Window Analysis

1. User sets **window size** (km) and **Te range** (min, max in km).  
2. An **interactive map** shows topography and Moho; a rectangle of the given window size follows the cursor.  
3. User **left-clicks** to fix the window.  
4. **Option (1) – One Te for whole region:**  
   - Extract topography anomaly and Moho undulation in that window.  
   - Run **one** Te inversion (minimize RMS over Te in the given range).  
   - Output: one \(T_e^*\), one RMS, plus optional “at upper bound” warning.  
5. **Option (2) – Te map within region:**  
   - Inside the same selected region, run a **moving window** with a smaller analysis window and step (shift).  
   - Each sub-window gives one Te; result is a **Te map** (and RMS map) over the selected region only.  
6. No shifting is involved in option (1); option (2) uses shifting only inside the selected area.

---

## 7. Moving-Window Analysis

1. User sets **window size** (km), **shift range** (min, max, step in km), and **Te range** (km).  
2. For each **shift distance** \(s\) in that range:
   - Place a grid of windows of the given size, spaced by \(s\) (in pixels: step = round(\(s\)/dx)).  
   - For each window: extract \(h'\) and \(w'\), run one Te inversion → one \(T_e^*\) and one RMS.  
   - Store results on the window-center grid → **Te map** and **RMS map** for that shift.  
3. Outputs: for each shift, PNG figures and Surfer `.grd` files for Te and RMS maps; 3D views optional.

---

## 8. Constants (Default: Mars)

| Symbol / name       | Typical value | Unit   | Meaning                    |
|---------------------|---------------|--------|----------------------------|
| \(\rho_{\mathrm{load}}\) | 2900          | kg/m³  | Crustal (load) density     |
| \(\rho_m\)          | 3500          | kg/m³  | Mantle density             |
| \(\rho_{\mathrm{infill}}\) | 2900 or 0   | kg/m³  | Infill (e.g. crust or air) |
| \(E\)               | \(10^{11}\)   | Pa     | Young’s modulus            |
| \(\nu\)             | 0.25          | –      | Poisson’s ratio            |
| \(g\)               | 3.72          | m/s²   | Surface gravity            |
| \(G\)               | \(6.674\times10^{-11}\) | m³/(kg·s²) | Gravitational constant |

Flexural rigidity: \(D = E\,T_e^3 / [12(1-\nu^2)]\).

---

## 9. Outputs

- **Folder:** `Output_YYYYMMDD_HHMMSS/` (and subfolder `3D/` if used).  
- **Data:**  
  - `inversion_results.npz`: arrays (topography, Moho, topo anomaly, Moho undulation, Te map(s), RMS map(s), coordinates, etc.).  
  - Surfer `.grd`: topography, Moho, anomalies, and (for single or moving window) Te and RMS grids.  
- **Figures (PNG):**  
  - `input_data.png`: topography and Moho (km axes).  
  - Single window: `single_window_location.png`, and either `single_window_result.png` (one Te) or `single_region_Te_map.png` / `single_region_rms_map.png` (Te map in region).  
  - Moving window: for each shift, `te_map_shift_XXkm.png`, `rms_map_shift_XXkm.png`, and optional 3D versions.  
- **Log:** `terminal_output.txt` in the same output folder.

---

## 10. Summary of Equations Used

| Step              | Equation / quantity |
|-------------------|----------------------|
| Flexural rigidity | \(D = \dfrac{E\,T_e^3}{12(1-\nu^2)}\) |
| **Spatial (Braitenberg)** | \(w_{\mathrm{pred}} = -(G * h')\) — convolution with flexural Green’s function |
| Flexure parameter | \(\Phi(k) = \dfrac{D\,k^4}{g\,(\rho_m - \rho_{\mathrm{infill}})}\) |
| Flexure filter (Fourier form of \(G\)) | \(F(k) = \dfrac{\rho_{\mathrm{load}}/(\rho_m - \rho_{\mathrm{infill}})}{1 + \Phi(k)}\) |
| Predicted Moho (this code: frequency domain) | \(W(\mathbf{k}) = -F(k)\,H(\mathbf{k})\), then \(w_{\mathrm{pred}} = \mathrm{IFFT}[W]\) |
| RMS misfit        | \(\mathrm{RMS} = \sqrt{\dfrac{1}{N}\sum (w'_{\mathrm{obs}} - w_{\mathrm{pred}})^2}\) |
| Inversion         | \(T_e^* = \arg\min_{T_e} \mathrm{RMS}(T_e)\) over \([T_{e,\min}, T_{e,\max}]\) |

**Note:** Braitenberg et al. (2002) use the **spatial-domain** convolution; this program uses the **equivalent frequency-domain** formulation (same result, efficient FFT implementation).
