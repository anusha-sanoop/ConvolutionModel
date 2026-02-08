"""
3D Moving Window Analysis for Elastic Thickness Estimation

This module implements the moving window technique for 3D volumetric data
to estimate spatial variations in elastic thickness.
"""

import numpy as np
import time
import os
import sys

# Add current directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from elastic_thickness_inversion_3d import ElasticThicknessInversion3D


class MovingWindowAnalysis3D:
    """
    Class for performing 3D moving window analysis of elastic thickness
    """

    def __init__(self, dx, dy, dz, rho_load=2900, rho_m=3500, rho_infill=2900, g=3.72):
        """
        Initialize 3D moving window analysis

        Parameters:
        -----------
        dx, dy, dz : float
            Grid spacing in meters
        rho_load : float
            Load density (kg/m³) - default 2900
        rho_m : float
            Mantle density (kg/m³) - default 3500
        rho_infill : float
            Infill density (kg/m³) - default 2900 (crustal compensation)
        g : float
            Gravity (m/s²) - default 3.72 (Mars)
        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rho_load = rho_load
        self.rho_m = rho_m
        self.rho_infill = rho_infill
        self.g = g

        print(f"\n3D Moving Window Analyzer Initialized:")
        print(
            f"  Using density parameters: ρ_load={rho_load}, ρ_m={rho_m}, ρ_infill={rho_infill}"
        )

    def analyze_3d(
        self,
        topography_anom_3d,
        moho_undulation_3d,
        window_size=1000000,
        window_depth=None,
        shift_distance=20000,
        shift_depth=None,
        Te_range=(5000, 80000),
        min_std_topo=100,
        min_std_moho=100,
    ):
        """
        Perform 3D moving window analysis

        Parameters:
        -----------
        topography_anom_3d : 3D array
            Topography anomaly volume (nz, ny, nx) - zero-mean
        moho_undulation_3d : 3D array
            Moho undulation volume (nz, ny, nx) - zero-mean
        window_size : float
            Size of moving window in X-Y plane (meters)
        window_depth : float, optional
            Depth extent of window (meters). If None, uses full depth
        shift_distance : float
            Distance to shift window in X-Y plane (meters)
        shift_depth : float, optional
            Depth shift (meters). If None, no depth shifting
        Te_range : tuple
            Range of Te values to search (m)
        min_std_topo : float
            Minimum standard deviation in topography for valid window
        min_std_moho : float
            Minimum standard deviation in Moho for valid window

        Returns:
        --------
        results : dict
            Dictionary containing Te_map, rms_map, and coordinates
        """
        nz, ny, nx = topography_anom_3d.shape

        # Calculate window positions
        window_pixels = int(window_size / self.dx)
        shift_pixels = int(shift_distance / self.dx)
        
        if window_depth is None:
            window_depth_pixels = nz
        else:
            window_depth_pixels = int(window_depth / self.dz)
        
        if shift_depth is None:
            shift_depth_pixels = nz  # No shifting
        else:
            shift_depth_pixels = int(shift_depth / self.dz)

        # Ensure window size is reasonable
        if window_pixels >= nx or window_pixels >= ny:
            raise ValueError(
                f"Window size ({window_pixels} pixels) is larger than grid size ({nx}x{ny})"
            )
        if window_depth_pixels > nz:
            raise ValueError(
                f"Window depth ({window_depth_pixels} pixels) is larger than depth dimension ({nz})"
            )

        # Starting positions
        x_positions = np.arange(0, nx - window_pixels, shift_pixels)
        y_positions = np.arange(0, ny - window_pixels, shift_pixels)
        
        if shift_depth_pixels < nz:
            z_positions = np.arange(0, nz - window_depth_pixels, shift_depth_pixels)
        else:
            z_positions = [0]  # Use full depth

        n_windows = len(x_positions) * len(y_positions) * len(z_positions)
        print("\n3D Moving Window Analysis:")
        print(f"  Window size (X-Y): {window_size / 1000:.0f} km ({window_pixels} pixels)")
        print(f"  Window depth: {window_depth_pixels} pixels")
        print(
            f"  Shift distance (X-Y): {shift_distance / 1000:.0f} km ({shift_pixels} pixels)"
        )
        print(f"  Number of windows: {n_windows}")
        print(
            f"  Te search range: {Te_range[0] / 1000:.0f}-{Te_range[1] / 1000:.0f} km"
        )

        # Initialize result arrays (using first z position for 2D map output)
        Te_map = np.full((len(y_positions), len(x_positions)), np.nan)
        rms_map = np.full((len(y_positions), len(x_positions)), np.nan)
        x_centers = np.zeros(len(x_positions))
        y_centers = np.zeros(len(y_positions))

        # Initialize inverter
        inverter = ElasticThicknessInversion3D(
            dx=self.dx,
            dy=self.dy,
            dz=self.dz,
            rho_load=self.rho_load,
            rho_m=self.rho_m,
            rho_infill=self.rho_infill,
            g=self.g,
        )

        # Process each window
        window_count = 0
        valid_count = 0
        start_time = time.time()

        for z_idx, z_start in enumerate(z_positions):
            z_end = z_start + window_depth_pixels
            
            for i, y_start in enumerate(y_positions):
                for j, x_start in enumerate(x_positions):
                    window_count += 1

                    # Extract 3D window data
                    y_end = y_start + window_pixels
                    x_end = x_start + window_pixels

                    topo_window_3d = topography_anom_3d[z_start:z_end, y_start:y_end, x_start:x_end]
                    moho_window_3d = moho_undulation_3d[z_start:z_end, y_start:y_end, x_start:x_end]

                    # Ensure windows are zero-mean
                    for z in range(topo_window_3d.shape[0]):
                        topo_window_3d[z, :, :] = topo_window_3d[z, :, :] - np.mean(topo_window_3d[z, :, :])
                        moho_window_3d[z, :, :] = moho_window_3d[z, :, :] - np.mean(moho_window_3d[z, :, :])

                    # Check if window has sufficient data variation
                    topo_std = np.std(topo_window_3d)
                    moho_std = np.std(moho_window_3d)
                    
                    if topo_std > min_std_topo and moho_std > min_std_moho:
                        try:
                            # Perform inversion
                            result = inverter.invert_elastic_thickness_3d(
                                topo_window_3d,
                                moho_window_3d,
                                Te_range=Te_range,
                                method="bounded",
                            )

                            # Store result (use first z position for 2D map)
                            if z_idx == 0:
                                Te_map[i, j] = result["Te_best"]
                                rms_map[i, j] = result["rms_best"]
                                valid_count += 1

                        except Exception as e:
                            if window_count % 100 == 0:
                                print(f"    Warning: Window {window_count} failed: {e}")
                            continue

                    # Store window center coordinates (only for first z)
                    if z_idx == 0:
                        x_centers[j] = (x_start + window_pixels // 2) * self.dx
                        y_centers[i] = (y_start + window_pixels // 2) * self.dy

                    # Progress update
                    if window_count % max(1, n_windows // 10) == 0:
                        elapsed = time.time() - start_time
                        progress = window_count / n_windows * 100
                        print(
                            f"    Progress: {progress:.1f}% ({window_count}/{n_windows}) - {elapsed:.1f}s"
                        )

        # Calculate statistics
        valid_Te = Te_map[~np.isnan(Te_map)]
        if len(valid_Te) > 0:
            print(
                f"\n  Results: {len(valid_Te)} valid windows ({valid_count / n_windows * 100:.1f}%)"
            )
            print(
                f"    Te range: {valid_Te.min() / 1000:.1f} - {valid_Te.max() / 1000:.1f} km"
            )
            print(
                f"    Te mean: {valid_Te.mean() / 1000:.1f} ± {valid_Te.std() / 1000:.1f} km"
            )
            print(f"    Te median: {np.median(valid_Te) / 1000:.1f} km")
        else:
            print("\n  Warning: No valid windows found!")

        return {
            "Te_map": Te_map,
            "rms_map": rms_map,
            "x_centers": x_centers,
            "y_centers": y_centers,
            "n_windows": n_windows,
            "n_valid": valid_count,
            "window_size": window_size,
            "window_depth": window_depth_pixels * self.dz if window_depth is not None else None,
            "shift_distance": shift_distance,
        }

    def analyze_multiple_shifts_3d(
        self,
        topography_anom_3d,
        moho_undulation_3d,
        window_size=1000000,
        window_depth=None,
        shift_min=20000,
        shift_max=80000,
        shift_step=20000,
        Te_range=(5000, 80000),
        min_std_topo=100,
        min_std_moho=100,
    ):
        """
        Perform 3D moving window analysis with multiple shift distances

        Parameters:
        -----------
        topography_anom_3d : 3D array
            Topography anomaly volume
        moho_undulation_3d : 3D array
            Moho undulation volume
        window_size : float
            Size of moving window in X-Y plane (meters)
        window_depth : float, optional
            Depth extent of window (meters)
        shift_min : float
            Minimum shift distance in X-Y plane (meters)
        shift_max : float
            Maximum shift distance in X-Y plane (meters)
        shift_step : float
            Step size for shift distance (meters)
        Te_range : tuple
            Range of Te values to search (m)
        min_std_topo : float
            Minimum standard deviation in topography for valid window
        min_std_moho : float
            Minimum standard deviation in Moho for valid window

        Returns:
        --------
        results : dict
            Dictionary with results for each shift distance
        """
        shift_distances = np.arange(shift_min, shift_max + shift_step, shift_step)
        all_results = {}

        for shift_dist in shift_distances:
            print("\n" + "=" * 60)
            print(f"Analyzing with shift distance: {shift_dist / 1000:.0f} km")
            print("=" * 60)

            result = self.analyze_3d(
                topography_anom_3d,
                moho_undulation_3d,
                window_size=window_size,
                window_depth=window_depth,
                shift_distance=shift_dist,
                Te_range=Te_range,
                min_std_topo=min_std_topo,
                min_std_moho=min_std_moho,
            )
            all_results[shift_dist] = result

        return all_results

