"""
Elastic Thickness Inversion using Braitenberg Convolution Method

Based on:
Braitenberg, C., Ebbing, J., & Götze, H. J. (2002).
Inverse modelling of elastic thickness by convolution method—the eastern Alps as a case example.
Earth and Planetary Science Letters, 202(2), 387-404.

This module implements the convolution method for estimating effective elastic thickness (Te)
of the lithosphere from topography and Moho depth data.
"""

import numpy as np
from scipy import fft
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from constants import (
    RHO_LOAD,
    RHO_MANTLE,
    RHO_INFILL,
    YOUNGS_MODULUS,
    POISSONS_RATIO,
    GRAVITY,
    DEFAULT_DX,
    DEFAULT_DY,
)


class ElasticThicknessInversion:
    """
    Class for inverse modeling of elastic thickness using convolution method
    Based on Braitenberg et al. (2002)
    """

    def __init__(
        self,
        dx=None,
        dy=None,
        rho_load=None,
        rho_m=None,
        rho_infill=None,
        E=None,
        nu=None,
        g=None,
    ):
        """
        Initialize the inversion class

        Parameters:
        -----------
        dx, dy : float
            Grid spacing in meters (default: 1000m)
        rho_load : float
            Load density (kg/m³) - default 2900 (crustal density for Mars)
        rho_m : float
            Mantle density (kg/m³) - default 3500 (Mars mantle)
        rho_infill : float
            Infill density (kg/m³) - default 0 (air) for topography above datum
            Use 2900 for subsurface loads or rho_water for submarine topography
        E : float
            Young's modulus (Pa) - default 1.0e11
        nu : float
            Poisson's ratio - default 0.25
        g : float
            Gravitational acceleration (m/s²) - default 3.72 (Mars)
        """
        # Use constants from constants.py if not provided
        self.dx = dx if dx is not None else DEFAULT_DX
        self.dy = dy if dy is not None else DEFAULT_DY
        self.rho_load = rho_load if rho_load is not None else RHO_LOAD
        self.rho_m = rho_m if rho_m is not None else RHO_MANTLE
        self.rho_infill = rho_infill if rho_infill is not None else RHO_INFILL
        self.E = E if E is not None else YOUNGS_MODULUS
        self.nu = nu if nu is not None else POISSONS_RATIO
        self.g = g if g is not None else GRAVITY

        # Flexural rigidity factor (D = D_factor * Te**3)
        self.D_factor = self.E / (12 * (1 - self.nu**2))

        # Calculate expected Airy compensation ratio
        # For topography in air: w/h = rho_crust / (rho_m - rho_infill)
        self.airy_ratio = self.rho_load / (self.rho_m - self.rho_infill)

        print("\nElastic Thickness Inverter Initialized:")
        print(f"  Load density: {self.rho_load} kg/m³")
        print(f"  Mantle density: {self.rho_m} kg/m³")
        print(f"  Infill density: {self.rho_infill} kg/m³")
        print(f"  Density contrast: {self.rho_m - self.rho_infill} kg/m³")
        print(f"  Airy compensation ratio: {self.airy_ratio:.3f}")
        print(f"  Young's modulus: {self.E:.2e} Pa")
        print(f"  Gravity: {self.g} m/s²")

    def calculate_flexure_filter(self, k, Te):
        """
        Calculate flexure filter function F(k) in wavenumber domain
        Based on Braitenberg et al. (2002), Eq. 3

        The flexure equation in wavenumber domain is:
        w(k) = h(k) × [ρ_c / (ρ_m - ρ_c)] × [1 / (1 + Φ(k))]

        where Φ(k) = (D × k^4) / (g × (ρ_m - ρ_c))

        For our case:
        - ρ_c = rho_load (crustal density)
        - ρ_m = rho_m (mantle density)
        - Compensating material = rho_infill

        Parameters:
        -----------
        k : array
            Wavenumber magnitude (rad/m)
        Te : float
            Elastic thickness (m)

        Returns:
        --------
        F_k : array
            Flexure filter function (ratio of Moho/Topo in frequency domain)
        """
        if Te < 1e-3:  # Treat Te=0 as pure Airy isostasy
            Te = 0.0

        D = self.D_factor * Te**3  # Flexural rigidity

        # Flexural parameter Φ(k) = (D k^4) / (g Δρ)
        Phi_k = (D * k**4) / (self.g * (self.rho_m - self.rho_infill))

        # Transfer function: F(k) = Airy_ratio / (1 + Φ(k))
        # At k=0: F(0) = Airy_ratio (pure isostatic compensation)
        # At high k: F→0 (plate stiffness dominates)

        F_k = self.airy_ratio / (1.0 + Phi_k)

        return F_k

    def predict_moho_flexure(self, topography_load, Te):
        """
        Forward model: calculate predicted Moho flexure from load and Te
        This follows the convolution method

        Parameters:
        -----------
        topography_load : 2D array
            Topographic load anomaly (m) - should be zero-mean
        Te : float
            Elastic thickness (m)

        Returns:
        --------
        moho_pred : 2D array
            Predicted Moho undulation (m) - zero-mean
            NEGATIVE values = downward deflection (deeper Moho)
        """
        ny, nx = topography_load.shape

        # Ensure zero-mean (critical for FFT)
        topography_load = topography_load - np.mean(topography_load)

        # Create wavenumber grids
        kx = 2 * np.pi * fft.fftfreq(nx, self.dx)
        ky = 2 * np.pi * fft.fftfreq(ny, self.dy)
        KX, KY = np.meshgrid(kx, ky)
        k = np.sqrt(KX**2 + KY**2)

        # FFT of topographic load
        load_fft = fft.fft2(topography_load)

        # Calculate flexure filter
        flexure_filter = self.calculate_flexure_filter(k, Te)

        # Set DC component (k=0) explicitly for numerical stability
        # At k=0, response should be Airy
        flexure_filter[0, 0] = self.airy_ratio

        # Predicted Moho in frequency domain W(k) = F(k) * H(k)
        # CRITICAL: Add negative sign for downward deflection under positive loads
        moho_fft = -flexure_filter * load_fft

        # Convert back to spatial domain w(r)
        moho_pred = np.real(fft.ifft2(moho_fft))

        # Ensure zero-mean output
        moho_pred = moho_pred - np.mean(moho_pred)

        return moho_pred

    def misfit_function(self, Te, topography_load, moho_obs, mask=None):
        """
        Calculate misfit between observed and predicted Moho undulations

        Parameters:
        -----------
        Te : float
            Elastic thickness (m)
        topography_load : 2D array
            Topography anomaly (m) - zero-mean
        moho_obs : 2D array
            Observed Moho undulations (m) - zero-mean
        mask : 2D array, optional
            Mask for valid data points

        Returns:
        --------
        rms : float
            RMS misfit (in meters of Moho)
        """
        moho_pred = self.predict_moho_flexure(topography_load, Te)

        if mask is not None:
            residual = (moho_obs - moho_pred)[mask]
        else:
            residual = moho_obs - moho_pred

        # DO NOT demean - we want to penalize amplitude mismatch!
        # The inputs are already zero-mean, so the residual measures
        # the actual difference in compensation amplitude

        rms = np.sqrt(np.mean(residual**2))
        return rms

    def invert_elastic_thickness(
        self,
        topography_load,
        moho_obs,
        Te_range=(1000, 100000),
        mask=None,
        method="bounded",
    ):
        """
        Invert for elastic thickness by minimizing Moho misfit

        Parameters:
        -----------
        topography_load : 2D array
            Topography anomaly (m) - should be zero-mean
        moho_obs : 2D array
            Observed Moho undulations (m) - should be zero-mean
        Te_range : tuple
            Range of Te values to search (m) - default (1, 100) km
        mask : 2D array, optional
            Mask for valid data points
        method : str
            Optimization method ('bounded' or 'grid_search')

        Returns:
        --------
        result : dict
            Inversion results containing Te_best, rms_best (in meters),
            and predicted moho.
        """
        if method == "bounded":
            # Bounded optimization
            result = minimize_scalar(
                self.misfit_function,
                bounds=Te_range,
                args=(topography_load, moho_obs, mask),
                method="bounded",
            )

            Te_best = result.x
            rms_best = result.fun

        elif method == "grid_search":
            # Grid search approach
            Te_values = np.linspace(Te_range[0], Te_range[1], 50)
            rms_values = []

            for Te in Te_values:
                rms = self.misfit_function(Te, topography_load, moho_obs, mask)
                rms_values.append(rms)

            rms_values = np.array(rms_values)
            best_idx = np.argmin(rms_values)
            Te_best = Te_values[best_idx]
            rms_best = rms_values[best_idx]

        else:
            raise ValueError("Method must be 'bounded' or 'grid_search'")

        # Calculate final predicted Moho
        moho_pred = self.predict_moho_flexure(topography_load, Te_best)

        return {
            "Te_best": Te_best,
            "rms_best": rms_best,
            "moho_pred": moho_pred,
            "method": method,
            "Te_range": Te_range,
        }

    def sensitivity_analysis(
        self, topography_load, moho_obs, Te_range=(1000, 100000), n_points=50, mask=None
    ):
        """
        Perform sensitivity analysis by testing different Te values

        Parameters:
        -----------
        topography_load : 2D array
            Topography anomaly (m)
        moho_obs : 2D array
            Observed Moho undulations (m)
        Te_range : tuple
            Range of Te values to test (m)
        n_points : int
            Number of Te values to test
        mask : 2D array, optional
            Mask for valid data points

        Returns:
        --------
        Te_values : array
            Tested Te values (m)
        rms_values : array
            Corresponding RMS misfits (in meters)
        """
        Te_values = np.linspace(Te_range[0], Te_range[1], n_points)
        rms_values = []

        print(
            f"Testing {n_points} Te values from {Te_range[0] / 1000:.1f} to {Te_range[1] / 1000:.1f} km..."
        )

        for i, Te in enumerate(Te_values):
            if i % 10 == 0:
                print(f"Progress: {i + 1}/{n_points}")

            rms = self.misfit_function(Te, topography_load, moho_obs, mask)
            rms_values.append(rms)

        return np.array(Te_values), np.array(rms_values)

    def diagnostic_check(self, topography_anom, moho_undulation):
        """
        Perform diagnostic checks on input data

        Parameters:
        -----------
        topography_anom : 2D array
            Topography anomaly (m)
        moho_undulation : 2D array
            Moho undulation (m)

        Returns:
        --------
        diagnostics : dict
            Diagnostic information
        """
        topo_std = np.std(topography_anom)
        moho_std = np.std(moho_undulation)

        # Expected Airy compensation
        expected_moho = topo_std * self.airy_ratio

        # Actual compensation ratio
        actual_ratio = moho_std / topo_std if topo_std > 0 else 0

        # Check sign correlation
        topo_flat = topography_anom.flatten()
        moho_flat = moho_undulation.flatten()
        correlation = np.corrcoef(topo_flat, moho_flat)[0, 1]

        # Rough Te estimate from compensation ratio
        # Lower ratio suggests higher Te (more flexural support)
        compensation_percentage = (actual_ratio / self.airy_ratio) * 100

        diagnostics = {
            "topo_std": topo_std,
            "moho_std": moho_std,
            "expected_moho_std": expected_moho,
            "actual_compensation_ratio": actual_ratio,
            "airy_compensation_ratio": self.airy_ratio,
            "compensation_percentage": compensation_percentage,
            "correlation": correlation,
        }

        print("\n DIAGNOSTIC CHECK ")
        print(f"Topography std: {topo_std:.1f} m")
        print(f"Moho undulation std: {moho_std:.1f} m")
        print(f"Expected Airy compensation: {expected_moho:.1f} m")
        print(f"Actual/Expected ratio: {actual_ratio / self.airy_ratio:.3f}")
        print(f"Compensation: {compensation_percentage:.1f}% of Airy")
        print(f"Topo-Moho correlation: {correlation:.3f}")

        # Check sign convention
        print("\n--- SIGN CONVENTION CHECK ---")
        if correlation > 0:
            print("⚠️  WARNING: Positive correlation detected!")
            print("   For loads (positive topo), Moho should deflect DOWN (negative)")
            print("   Your data shows: High topo → Deep Moho (both positive)")
            print("   → This is CORRECT if Moho is stored as negative depth")
        elif correlation < -0.5:
            print("✓ Negative correlation: High topo → Shallow Moho")
            print("  This suggests INVERTED sign convention")
            print("  → Try flipping Moho sign: moho_undulation *= -1")
        else:
            print("⚠️  Weak correlation - check data quality")

        if compensation_percentage > 90:
            print("\n→ Suggests weak/thin lithosphere (low Te)")
        elif compensation_percentage < 50:
            print("\n→ Suggests strong/thick lithosphere (high Te)")
        else:
            print("\n→ Suggests intermediate flexural support")

        print("\n")

        return diagnostics

    def plot_flexure_response(self, Te_values=[5000, 20000, 50000, 100000]):
        """
        Plot flexure filter response for different Te values

        Parameters:
        -----------
        Te_values : list
            List of Te values to plot (in meters)
        """
        # Create wavelength array
        wavelengths = np.logspace(3, 6, 100)  # 1 km to 1000 km
        k = 2 * np.pi / wavelengths

        fig, ax = plt.subplots(figsize=(10, 6))

        for Te in Te_values:
            F_k = self.calculate_flexure_filter(k, Te)
            ax.loglog(
                wavelengths / 1000, F_k, linewidth=2, label=f"Te = {Te / 1000:.0f} km"
            )

        # Add Airy isostasy line
        airy_response = np.ones_like(k) * self.airy_ratio
        ax.loglog(
            wavelengths / 1000, airy_response, "k--", linewidth=2, label="Airy (Te=0)"
        )

        ax.set_xlabel("Wavelength (km)", fontsize=12)
        ax.set_ylabel("Flexure Response (m Moho / m Topo)", fontsize=12)
        ax.set_title("Flexural Response Function", fontsize=14, fontweight="bold")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10)
        ax.set_ylim([0.1, 10])

        plt.tight_layout()
        return fig

    def diagnostic_check(self, topography_anom, moho_undulation):
        """
        Perform diagnostic checks on input data

        Parameters:
        -----------
        topography_anom : 2D array
            Topography anomaly (m)
        moho_undulation : 2D array
            Moho undulation (m)

        Returns:
        --------
        diagnostics : dict
            Diagnostic information
        """
        topo_std = np.std(topography_anom)
        moho_std = np.std(moho_undulation)

        # Expected Airy compensation
        expected_moho = topo_std * self.airy_ratio

        # Actual compensation ratio
        actual_ratio = moho_std / topo_std if topo_std > 0 else 0

        # Rough Te estimate from compensation ratio
        # Lower ratio suggests higher Te (more flexural support)
        compensation_percentage = (actual_ratio / self.airy_ratio) * 100

        diagnostics = {
            "topo_std": topo_std,
            "moho_std": moho_std,
            "expected_moho_std": expected_moho,
            "actual_compensation_ratio": actual_ratio,
            "airy_compensation_ratio": self.airy_ratio,
            "compensation_percentage": compensation_percentage,
        }

        print("\n DIAGNOSTIC CHECK ")
        print(f"Topography std: {topo_std:.1f} m")
        print(f"Moho undulation std: {moho_std:.1f} m")
        print(f"Expected Airy compensation: {expected_moho:.1f} m")
        print(f"Actual/Expected ratio: {actual_ratio / self.airy_ratio:.3f}")
        print(f"Compensation: {compensation_percentage:.1f}% of Airy")

        if compensation_percentage > 90:
            print("→ Suggests weak/thin lithosphere (low Te)")
        elif compensation_percentage < 50:
            print("→ Suggests strong/thick lithosphere (high Te)")
        else:
            print("→ Suggests intermediate flexural support")

        print("\n")

        return diagnostics
