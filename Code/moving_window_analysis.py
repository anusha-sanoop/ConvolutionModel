"""
Moving Window Analysis for Elastic Thickness Estimation

This module implements the moving window technique to estimate spatial
variations in elastic thickness across the study area.
"""

import numpy as np
import time
from elastic_thickness_inversion import ElasticThicknessInversion
from constants import (
    RHO_LOAD,
    RHO_MANTLE,
    RHO_INFILL,
    GRAVITY,
)


class MovingWindowAnalysis:
    """
    Class for performing moving window analysis of elastic thickness
    """

    def __init__(self, dx, dy, rho_load=None, rho_m=None, rho_infill=None, g=None):
        """
        Initialize moving window analysis

        Parameters:
        -----------
        dx, dy : float
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
        # Use constants from constants.py if not provided
        self.rho_load = rho_load if rho_load is not None else RHO_LOAD
        self.rho_m = rho_m if rho_m is not None else RHO_MANTLE
        self.rho_infill = rho_infill if rho_infill is not None else RHO_INFILL
        self.g = g if g is not None else GRAVITY

    def analyze(
        self,
        topography_anom,
        moho_undulation,
        window_size=1000000,
        shift_distance=20000,
        Te_range=(5000, 80000),
        min_std_topo=0,
        min_std_moho=0,
    ):
        """
        Perform moving window analysis

        Parameters:
        -----------
        topography_anom : 2D array
            Topography anomaly (m) - zero-mean
        moho_undulation : 2D array
            Moho undulation (m) - zero-mean
        window_size : float
            Size of moving window in meters (default: 1000 km)
        shift_distance : float
            Distance to shift window in meters (default: 20 km)
        Te_range : tuple
            Range of Te values to search (m)
        min_std_topo : float
            Minimum standard deviation in topography for valid window (default 0 = all windows)
        min_std_moho : float
            Minimum standard deviation in Moho for valid window (default 0 = all windows)

        Returns:
        --------
        results : dict
            Dictionary containing Te_map, rms_map, and coordinates
        """
        ny, nx = topography_anom.shape

        # Calculate window positions
        window_pixels = int(window_size / self.dx)
        shift_pixels = int(shift_distance / self.dx)

        # Ensure minimum values (at least 1 pixel)
        if window_pixels < 1:
            window_pixels = 1
        
        if shift_pixels < 1:
            shift_pixels = 1

        # Ensure window size is reasonable
        if window_pixels >= nx or window_pixels >= ny:
            raise ValueError(
                f"Window size ({window_pixels} pixels) is larger than grid size ({nx}x{ny})"
            )

        # Starting positions
        x_positions = np.arange(0, nx - window_pixels, shift_pixels)
        y_positions = np.arange(0, ny - window_pixels, shift_pixels)

        n_windows = len(x_positions) * len(y_positions)

        # Initialize result arrays
        Te_map = np.full((len(y_positions), len(x_positions)), np.nan)
        rms_map = np.full((len(y_positions), len(x_positions)), np.nan)
        x_centers = np.zeros(len(x_positions))
        y_centers = np.zeros(len(y_positions))

        # Initialize inverter
        inverter = ElasticThicknessInversion(
            dx=self.dx,
            dy=self.dy,
            rho_load=self.rho_load,
            rho_m=self.rho_m,
            rho_infill=self.rho_infill,
            g=self.g,
        )

        # Process each window
        window_count = 0
        valid_count = 0
        start_time = time.time()

        for i, y_start in enumerate(y_positions):
            for j, x_start in enumerate(x_positions):
                window_count += 1

                # Extract window data
                y_end = y_start + window_pixels
                x_end = x_start + window_pixels

                topo_window = topography_anom[y_start:y_end, x_start:x_end]
                moho_window = moho_undulation[y_start:y_end, x_start:x_end]

                # Ensure windows are zero-mean
                topo_window = topo_window - np.mean(topo_window)
                moho_window = moho_window - np.mean(moho_window)

                # Check if window has sufficient data variation
                if (
                    np.std(topo_window) > min_std_topo
                    and np.std(moho_window) > min_std_moho
                ):
                    try:
                        # Perform inversion
                        result = inverter.invert_elastic_thickness(
                            topo_window,
                            moho_window,
                            Te_range=Te_range,
                            method="bounded",
                        )

                        Te_map[i, j] = result["Te_best"]
                        rms_map[i, j] = result["rms_best"]
                        valid_count += 1

                    except Exception:
                        continue

                # Store window center coordinates IN METERS (not pixel indices)
                x_centers[j] = (x_start + window_pixels // 2) * self.dx
                y_centers[i] = (y_start + window_pixels // 2) * self.dy

                # Progress update (no terminal output; timing kept for potential future use)
                if window_count % max(1, n_windows // 10) == 0:
                    _elapsed = time.time() - start_time

        # Calculate statistics
        valid_Te = Te_map[~np.isnan(Te_map)]

        return {
            "Te_map": Te_map,
            "rms_map": rms_map,
            "x_centers": x_centers,
            "y_centers": y_centers,
            "n_windows": n_windows,
            "n_valid": valid_count,
            "window_size": window_size,
            "shift_distance": shift_distance,
        }

    def analyze_multiple_shifts(
        self,
        topography_anom,
        moho_undulation,
        window_size=1000000,
        shift_min=20000,
        shift_max=80000,
        shift_step=20000,
        Te_range=(5000, 80000),
        min_std_topo=0,
        min_std_moho=0,
    ):
        """
        Perform moving window analysis with multiple shift distances

        Parameters:
        -----------
        topography_anom : 2D array
            Topography anomaly (m) - zero-mean
        moho_undulation : 2D array
            Moho undulation (m) - zero-mean
        window_size : float
            Size of moving window in meters (default: 1000 km)
        shift_min : float
            Minimum shift distance in meters (default: 20 km)
        shift_max : float
            Maximum shift distance in meters (default: 80 km)
        shift_step : float
            Step size for shift distance in meters (default: 20 km)
        Te_range : tuple
            Range of Te values to search (m)
        min_std_topo : float
            Minimum standard deviation in topography for valid window (default 0 = all windows)
        min_std_moho : float
            Minimum standard deviation in Moho for valid window (default 0 = all windows)

        Returns:
        --------
        results : dict
            Dictionary with results for each shift distance
        """
        shift_distances = np.arange(shift_min, shift_max + shift_step, shift_step)
        all_results = {}

        for shift_dist in shift_distances:
            result = self.analyze(
                topography_anom,
                moho_undulation,
                window_size=window_size,
                shift_distance=shift_dist,
                Te_range=Te_range,
                min_std_topo=min_std_topo,
                min_std_moho=min_std_moho,
            )
            all_results[shift_dist] = result

        return all_results
