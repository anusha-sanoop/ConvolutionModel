"""
3D Elastic Thickness Inversion using Braitenberg Convolution Method

This module extends the 2D convolution method to 3D volumetric data.
Works with 3D topography and Moho depth volumes.
"""

import numpy as np
from scipy import fft
from scipy.optimize import minimize_scalar


class ElasticThicknessInversion3D:
    """
    Class for 3D inverse modeling of elastic thickness using convolution method
    Based on Braitenberg et al. (2002), extended to 3D
    """

    def __init__(
        self,
        dx=1000,
        dy=1000,
        dz=1000,
        rho_load=2900,
        rho_m=3500,
        rho_infill=0,
        E=1.0e11,
        nu=0.25,
        g=3.72,
    ):
        """
        Initialize the 3D inversion class

        Parameters:
        -----------
        dx, dy, dz : float
            Grid spacing in meters (default: 1000m)
        rho_load : float
            Load density (kg/m³) - default 2900 (crustal density for Mars)
        rho_m : float
            Mantle density (kg/m³) - default 3500 (Mars mantle)
        rho_infill : float
            Infill density (kg/m³) - default 0 (air) for topography above datum
        E : float
            Young's modulus (Pa) - default 1.0e11
        nu : float
            Poisson's ratio - default 0.25
        g : float
            Gravitational acceleration (m/s²) - default 3.72 (Mars)
        """
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.rho_load = rho_load
        self.rho_m = rho_m
        self.rho_infill = rho_infill
        self.E = E
        self.nu = nu
        self.g = g

        # Flexural rigidity factor (D = D_factor * Te**3)
        self.D_factor = self.E / (12 * (1 - self.nu**2))

        # Calculate expected Airy compensation ratio
        self.airy_ratio = self.rho_load / (self.rho_m - self.rho_infill)

        print("\n3D Elastic Thickness Inverter Initialized:")
        print(f"  Grid spacing: dx={dx}m, dy={dy}m, dz={dz}m")
        print(f"  Load density: {self.rho_load} kg/m³")
        print(f"  Mantle density: {self.rho_m} kg/m³")
        print(f"  Infill density: {self.rho_infill} kg/m³")
        print(f"  Density contrast: {self.rho_m - self.rho_infill} kg/m³")
        print(f"  Airy compensation ratio: {self.airy_ratio:.3f}")
        print(f"  Young's modulus: {self.E:.2e} Pa")
        print(f"  Gravity: {self.g} m/s²")

    def calculate_flexure_filter_3d(self, k, Te):
        """
        Calculate 3D flexure filter function F(k) in wavenumber domain
        
        Parameters:
        -----------
        k : array
            Wavenumber magnitude (rad/m) - 3D array
        Te : float
            Elastic thickness (m)
        
        Returns:
        --------
        F_k : array
            Flexure filter function (3D array)
        """
        if Te < 1e-3:
            Te = 0.0

        D = self.D_factor * Te**3  # Flexural rigidity

        # Flexural parameter Φ(k) = (D k^4) / (g Δρ)
        Phi_k = (D * k**4) / (self.g * (self.rho_m - self.rho_infill))

        # Transfer function: F(k) = Airy_ratio / (1 + Φ(k))
        F_k = self.airy_ratio / (1.0 + Phi_k)

        return F_k

    def predict_moho_flexure_3d(self, topography_load_3d, Te):
        """
        Forward model: calculate predicted Moho flexure from 3D load and Te
        
        Parameters:
        -----------
        topography_load_3d : 3D array
            Topographic load volume (nz, ny, nx) - should be zero-mean
        Te : float
            Elastic thickness (m)
        
        Returns:
        --------
        moho_pred_3d : 3D array
            Predicted Moho undulation (nz, ny, nx) - zero-mean
        """
        nz, ny, nx = topography_load_3d.shape

        # Ensure zero-mean for each depth level
        topography_load_3d = topography_load_3d.copy()
        for z in range(nz):
            topography_load_3d[z, :, :] = topography_load_3d[z, :, :] - np.mean(
                topography_load_3d[z, :, :]
            )

        # Create 3D wavenumber grids
        kx = 2 * np.pi * fft.fftfreq(nx, self.dx)
        ky = 2 * np.pi * fft.fftfreq(ny, self.dy)
        kz = 2 * np.pi * fft.fftfreq(nz, self.dz)
        
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        k = np.sqrt(KX**2 + KY**2 + KZ**2)

        # 3D FFT of topographic load
        load_fft_3d = fft.fftn(topography_load_3d)

        # Calculate flexure filter
        flexure_filter_3d = self.calculate_flexure_filter_3d(k, Te)

        # Set DC component explicitly
        flexure_filter_3d[0, 0, 0] = self.airy_ratio

        # Predicted Moho in frequency domain
        moho_fft_3d = -flexure_filter_3d * load_fft_3d

        # Convert back to spatial domain
        moho_pred_3d = np.real(fft.ifftn(moho_fft_3d))

        # Ensure zero-mean for each depth level
        for z in range(nz):
            moho_pred_3d[z, :, :] = moho_pred_3d[z, :, :] - np.mean(moho_pred_3d[z, :, :])

        return moho_pred_3d

    def misfit_function_3d(self, Te, topography_load_3d, moho_obs_3d, mask_3d=None):
        """
        Calculate misfit between observed and predicted 3D Moho undulations
        
        Parameters:
        -----------
        Te : float
            Elastic thickness (m)
        topography_load_3d : 3D array
            Topography volume (nz, ny, nx) - zero-mean
        moho_obs_3d : 3D array
            Observed Moho undulations (nz, ny, nx) - zero-mean
        mask_3d : 3D array, optional
            Mask for valid data points
        
        Returns:
        --------
        rms : float
            RMS misfit (in meters of Moho)
        """
        moho_pred_3d = self.predict_moho_flexure_3d(topography_load_3d, Te)

        if mask_3d is not None:
            residual = (moho_obs_3d - moho_pred_3d)[mask_3d]
        else:
            residual = moho_obs_3d - moho_pred_3d

        rms = np.sqrt(np.mean(residual**2))
        return rms

    def invert_elastic_thickness_3d(
        self,
        topography_load_3d,
        moho_obs_3d,
        Te_range=(1000, 100000),
        mask_3d=None,
        method="bounded",
    ):
        """
        Invert for elastic thickness from 3D data
        
        Parameters:
        -----------
        topography_load_3d : 3D array
            Topography volume (nz, ny, nx) - should be zero-mean
        moho_obs_3d : 3D array
            Observed Moho undulations (nz, ny, nx) - should be zero-mean
        Te_range : tuple
            Range of Te values to search (m)
        mask_3d : 3D array, optional
            Mask for valid data points
        method : str
            Optimization method ('bounded' or 'grid_search')
        
        Returns:
        --------
        result : dict
            Inversion results containing Te_best, rms_best, and predicted moho
        """
        if method == "bounded":
            result = minimize_scalar(
                self.misfit_function_3d,
                bounds=Te_range,
                args=(topography_load_3d, moho_obs_3d, mask_3d),
                method="bounded",
            )

            Te_best = result.x
            rms_best = result.fun

        elif method == "grid_search":
            Te_values = np.linspace(Te_range[0], Te_range[1], 50)
            rms_values = []

            for Te in Te_values:
                rms = self.misfit_function_3d(Te, topography_load_3d, moho_obs_3d, mask_3d)
                rms_values.append(rms)

            rms_values = np.array(rms_values)
            best_idx = np.argmin(rms_values)
            Te_best = Te_values[best_idx]
            rms_best = rms_values[best_idx]

        else:
            raise ValueError("Method must be 'bounded' or 'grid_search'")

        # Calculate final predicted Moho
        moho_pred_3d = self.predict_moho_flexure_3d(topography_load_3d, Te_best)

        return {
            "Te_best": Te_best,
            "rms_best": rms_best,
            "moho_pred_3d": moho_pred_3d,
            "method": method,
            "Te_range": Te_range,
        }

    def sensitivity_analysis_3d(
        self, topography_load_3d, moho_obs_3d, Te_range=(1000, 100000), n_points=50, mask_3d=None
    ):
        """
        Perform sensitivity analysis for 3D data
        
        Parameters:
        -----------
        topography_load_3d : 3D array
            Topography volume
        moho_obs_3d : 3D array
            Observed Moho undulations
        Te_range : tuple
            Range of Te values to test (m)
        n_points : int
            Number of Te values to test
        mask_3d : 3D array, optional
            Mask for valid data points
        
        Returns:
        --------
        Te_values : array
            Tested Te values (m)
        rms_values : array
            Corresponding RMS misfits (meters)
        """
        Te_values = np.linspace(Te_range[0], Te_range[1], n_points)
        rms_values = []

        print(
            f"Testing {n_points} Te values from {Te_range[0] / 1000:.1f} to {Te_range[1] / 1000:.1f} km..."
        )

        for i, Te in enumerate(Te_values):
            if i % 10 == 0:
                print(f"Progress: {i + 1}/{n_points}")

            rms = self.misfit_function_3d(Te, topography_load_3d, moho_obs_3d, mask_3d)
            rms_values.append(rms)

        return np.array(Te_values), np.array(rms_values)

