"""
Gravity-based Moho Depth Estimation using Parker-Oldenburg Method

Based on:
Parker, R. L. (1973). The rapid calculation of potential anomalies.
Geophysical Journal International, 31(4), 447-455.

Oldenburg, D. W. (1974). The inversion and interpretation of gravity anomalies.
Geophysics, 39(4), 526-536.
"""

import numpy as np
from scipy import fft
from constants import (
    DENSITY_CONTRAST,
    DEFAULT_REFERENCE_MOHO_DEPTH,
    GRAVITATIONAL_CONSTANT,
    DEFAULT_DX,
    DEFAULT_DY,
)


class GravityMohoInversion:
    """
    Class for estimating Moho depth from Bouguer gravity anomalies
    using Parker-Oldenburg frequency domain method
    """

    def __init__(
        self,
        dx=None,
        dy=None,
        density_contrast=None,
        reference_moho=None,
        G=None,
    ):
        """
        Initialize gravity inversion

        Parameters:
        -----------
        dx, dy : float
            Grid spacing in meters
        density_contrast : float
            Density contrast between crust and mantle (kg/m³)
            Default: 600 kg/m³ (typical for Mars: 3500 - 2900)
        reference_moho : float
            Reference Moho depth in meters (default: 50 km)
        G : float
            Gravitational constant (m³/kg/s²)
        """
        # Use constants from constants.py if not provided
        self.dx = dx if dx is not None else DEFAULT_DX
        self.dy = dy if dy is not None else DEFAULT_DY
        self.density_contrast = density_contrast if density_contrast is not None else DENSITY_CONTRAST
        # reference_moho should be provided by user or calculated from data, not from constants
        # DEFAULT_REFERENCE_MOHO_DEPTH is only a fallback default
        self.reference_moho = reference_moho if reference_moho is not None else DEFAULT_REFERENCE_MOHO_DEPTH
        self.G = G if G is not None else GRAVITATIONAL_CONSTANT

        print("\nGravity Moho Inverter Initialized:")
        print(f"  Density contrast: {density_contrast} kg/m³")
        print(f"  Reference Moho depth: {reference_moho / 1000:.1f} km")
        print(f"  Gravitational constant: {G:.3e} m³/kg/s²")

    def predict_moho_from_gravity(self, bouguer_gravity, max_iterations=10, tolerance=1.0):
        """
        Predict Moho depth from Bouguer gravity anomaly using Parker-Oldenburg method

        Parameters:
        -----------
        bouguer_gravity : 2D array
            Bouguer gravity anomaly (mGal)
        max_iterations : int
            Maximum number of iterations for convergence
        tolerance : float
            Convergence tolerance (m)

        Returns:
        --------
        moho_depth : 2D array
            Predicted Moho depth (m, negative = below surface)
        """
        ny, nx = bouguer_gravity.shape

        # Ensure zero-mean
        bouguer_gravity = bouguer_gravity - np.mean(bouguer_gravity)

        # Convert mGal to m/s² (1 mGal = 1e-5 m/s²)
        gravity_si = bouguer_gravity * 1e-5

        # Create wavenumber grids
        kx = 2 * np.pi * fft.fftfreq(nx, self.dx)
        ky = 2 * np.pi * fft.fftfreq(ny, self.dy)
        KX, KY = np.meshgrid(kx, ky)
        k = np.sqrt(KX**2 + KY**2)

        # Avoid division by zero at k=0
        k[k == 0] = 1e-10

        # FFT of gravity anomaly
        gravity_fft = fft.fft2(gravity_si)

        # Initialize Moho depth with reference depth
        moho_depth = np.full((ny, nx), -self.reference_moho)

        # Iterative Parker-Oldenburg inversion
        for iteration in range(max_iterations):
            # Calculate gravity from current Moho estimate
            moho_undulation = moho_depth - (-self.reference_moho)
            moho_fft = fft.fft2(moho_undulation)

            # Forward model: gravity from Moho
            # G(k) = 2πG Δρ * H(k) * exp(-k*d0)
            predicted_gravity_fft = (
                2 * np.pi * self.G * self.density_contrast * moho_fft * np.exp(-k * self.reference_moho)
            )

            # Calculate residual
            residual_fft = gravity_fft - predicted_gravity_fft

            # Update Moho depth
            # H(k) = G_residual(k) / (2πG Δρ * exp(-k*d0))
            update_fft = residual_fft / (
                2 * np.pi * self.G * self.density_contrast * np.exp(-k * self.reference_moho)
            )
            update = np.real(fft.ifft2(update_fft))

            # Update Moho depth
            moho_depth = moho_depth + update

            # Check convergence
            max_change = np.max(np.abs(update))
            if max_change < tolerance:
                print(f"  Converged after {iteration + 1} iterations (max change: {max_change:.2f} m)")
                break

        if iteration == max_iterations - 1:
            print(f"  Warning: Reached maximum iterations ({max_iterations})")

        return moho_depth

    def predict_moho_simple(self, bouguer_gravity):
        """
        Simplified gravity-to-Moho conversion (linear approximation)
        For quick estimates without iteration

        Parameters:
        -----------
        bouguer_gravity : 2D array
            Bouguer gravity anomaly (mGal)

        Returns:
        --------
        moho_depth : 2D array
            Predicted Moho depth (m)
        """
        # Convert mGal to m/s²
        gravity_si = bouguer_gravity * 1e-5

        # Simple relationship: Δg ≈ 2πG Δρ h
        # For small undulations: h ≈ Δg / (2πG Δρ)
        moho_undulation = gravity_si / (2 * np.pi * self.G * self.density_contrast)

        # Convert to absolute depth
        moho_depth = -self.reference_moho + moho_undulation

        return moho_depth
