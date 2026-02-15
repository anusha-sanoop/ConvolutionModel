# Changelog – main.py and related updates

**Keep this file updated:** Add a new numbered section (and date if useful) whenever code or behaviour changes are made. Summarise what was changed and why.

---

Summary of changes made during development (conversation session).

---

## 1. Display all outputs in km scale (no meters)

- **Plot labels and data**
  - `plot_final_results`: Topography and Moho depth plotted in km (÷1000); titles and colorbars set to "Topography (km)", "Moho depth (km)".
  - `plot_input_data_3d`: X, Y, elevation/depth converted to km; axis and colorbar labels use "(km)".
  - Single-window selection and location figures: data in km; titles "Topography (km)", "Moho depth (km)".
  - Te and RMS map figures (2D/3D): extent and axes in km; axis labels "X (km)", "Y (km)"; Te and RMS already or converted to km.
- **Print / terminal output**
  - Reference Moho, Moho depth range, data statistics (topo mean/std, Moho mean/std, undulation std), sign-convention check (max topo, Moho undulation), global and single-window RMS misfit, confidence-interval threshold, and single-window result figure text all switched from meters to km.
- Internal calculations remain in meters; only displayed values and labels use km.

---

## 2. Fix undefined variables in plot_input_data_3d

- **Issue:** `moho_undulation_km` and `topo_anom_km` were used (mean reference surface, z-limits) but not defined.
- **Change:** At the start of `plot_input_data_3d`, added:
  - `x_km`, `y_km` (coordinates in km),
  - `Xg`, `Yg` from meshgrid of those,
  - `topo_anom_km = topo_anom / 1000.0`,
  - `moho_undulation_km = moho_undulation / 1000.0`.
- First two `plot_surface` calls now use `topo_anom_km` and `moho_undulation_km`.

---

## 3. plot_final_results: ensure topo and Moho are shown in km

- **Issue:** Input-data figure still showed "Topography (m)" and "Moho depth (m)" with meter-scaled colorbars.
- **Change:** In `plot_final_results`:
  - Build `topo_km = topography / 1000.0`, `moho_km = moho_depth / 1000.0`.
  - Use `topo_km` and `moho_km` in `imshow` with percentile-based limits in km.
  - Titles and colorbar labels set to "Topography (km)" and "Moho depth (km)" (and "[from gravity]" variant for Moho when applicable).

---

## 4. Trim borders of Te and RMS maps (unreliable edges)

- **Goal:** Remove outer cells of Te/RMS maps where estimates are less reliable (window near data edge).
- **Implementation:**
  - `TRIM_BORDER_CELLS = 2`: number of cells removed from each side of Te and RMS grids.
  - Helper `trim_te_result(res, trim)`: returns a copy of the result with `Te_map`, `rms_map`, `x_centers`, `y_centers` sliced by `trim` (no trim if grid too small).
  - `mw_results_trimmed`: built from `mw_results_dict` by applying `trim_te_result(..., TRIM_BORDER_CELLS)`.
- **Usage:** All Te and RMS map figures, saved arrays in the .npz, and Surfer .grd exports for Te and RMS use `mw_results_trimmed` instead of `mw_results_dict`.
- Print message added: "Trimmed N border cells from each side of Te/RMS maps (edge estimates excluded)."

---

## 5. RMS color scale after trimming

- **Issue:** After trimming, the RMS colorbar range dropped (e.g. max from ~20 km to ~7 km) because high-RMS edge cells were removed.
- **Change:**
  - `create_rms_map_figure` now accepts optional `rms_vmin` and `rms_vmax` (in km).
  - When provided, the figure uses these as the color scale limits so the bar stays comparable (e.g. 0–20 km).
  - Call site: before calling `create_rms_map_figure`, compute 5th and 95th percentiles of the **untrimmed** RMS map (in km) for that shift and pass them as `rms_vmin`, `rms_vmax`.
- Result: Trimmed RMS map is still plotted (smaller area), but the colorbar range matches the full-map scale.

---

## 6. Remove explanatory comments (per user request)

- Removed the long comment block in `create_rms_map_figure` that described how RMS is plotted and how the color scale works (including trimmed vs full map).
- Removed duplicate/orphaned code that had been introduced during an earlier edit; left a single, correct flow for RMS figure creation and the following "Creating RMS maps." loop.

---

## 7. Restore Te map figure creation (Te maps not saved)

- **Issue:** "Creating Te maps." was printed but no Te map PNGs were saved; only RMS map PNGs appeared in the output folder.
- **Cause:** The loop that calls `create_te_map_figure` for each shift had been removed when cleaning up the duplicate block.
- **Change:** Re-added the loop after the `create_rms_map_figure` definition and before "Creating RMS maps.":
  - For each `(shift_dist, result)` in `mw_results_trimmed.items()`, call `create_te_map_figure(result, shift_dist, X_topo, Y_topo, Te_min, Te_max, output_folder, output_folder_3d, filename_prefix)`.
- Result: For each shift, 2D and 3D Te map PNGs are saved (e.g. `..._te_map_shift_60km.png` and `..._te_map_shift_60km_3d.png`) in the main output folder and in the `3D` subfolder.

---

## 8. Workflow documentation updated to match program

- **Code/PROGRAM_WORKFLOW_AND_EQUATIONS.md**
  - Analysis flow: optional global inversion; then choice Single (1) vs Moving (2). Single = interactive click; Moving = Te + RMS maps with border trimming.
  - High-level workflow: reference Moho in km; optional global step; trim border cells and RMS scale from full map noted for moving window.
  - Data loading: clarified that grid data are in meters, all displayed outputs in km.
  - Outputs (Section 9): display units km; trimmed Te/RMS maps; list of figures (input_data.png, single-window outputs, moving-window te_map_shift_XXkm, rms_map_shift_XXkm, 3D subfolder); RMS color scale and trimming noted.
- **Code/WORKFLOW_FLOWCHART.md**
  - Step 1: Reference Moho depth (m) → (km).
  - Step 9: Replaced single “Moving window” path with “CHOOSE WINDOW MODE” (1 = single, 2 = moving); added short Single vs Moving branches; Step 9b (moving only) includes trim, Te + RMS maps (2D and 3D), RMS scale from full map.
  - Removed “STEP 10: CALCULATE RESIDUAL MOHO”; final plot is 1×2 (Topography and Moho in km), not 2×2 with residual.
  - Step 11 (Save): .npz and .grd content; trimmed Te/RMS; terminal_output.txt.
  - Output files list: input_data.png 1×2 in km; te_map and rms_map PNGs (trimmed); 3D subfolder; .npz and .grd; terminal_output.txt.

---

## 9. Code and paper guide added

- **Code/CODE_AND_PAPER_GUIDE.md** – New detailed guide to understand Braitenberg’s paper and the code together:
  - Braitenberg et al. (2002) concept: convolution method, thin-plate flexure, inverse for \(T_e\).
  - Paper’s spatial convolution and equivalence to frequency-domain \(W = -F(k)H(k)\); equations for \(\Phi(k)\), \(F(k)\), \(D\), Airy ratio.
  - Module-by-module: **constants** (physical parameters, \(D\), Airy ratio); **data_loader** (`read_surfer_grd`, `write_surfer_grd`, `check_grid_compatibility`, `apply_taper`, `prepare_data_for_inversion`) with equations; **elastic_thickness_inversion** (each method: `calculate_flexure_filter`, `predict_moho_flexure`, `misfit_function`, `invert_elastic_thickness`, `sensitivity_analysis`, `diagnostic_check`) with equations; **gravity_moho_inversion** (Parker–Oldenburg, forward/inverse); **moving_window_analysis** (`analyze`, `analyze_multiple_shifts`); **main.py** flow.
  - Equation summary table and “Reading the paper with the code” section.

---

## File modified

- **Code/main.py** – all of the above code changes (sections 1–7).
- **Code/PROGRAM_WORKFLOW_AND_EQUATIONS.md**, **Code/WORKFLOW_FLOWCHART.md** – workflow docs updated (section 8).
- **Code/CODE_AND_PAPER_GUIDE.md** – new (section 9).

---

*Last updated: 2026-02-14 (session summary).*
