"""
Flexural Response Function Calculations
Implements the unit impulse response for thin elastic plate flexure
"""

import numpy as np
from scipy import special
from typing import Optional
from grd_reader import GridData


class FlexureResponse:
    """Class for computing flexural response functions"""
    
    @staticmethod
    def compute_flexural_rigidity(Te: float, E: float = 70e9, nu: float = 0.25) -> float:
        """
        Compute flexural rigidity from elastic thickness
        
        D = E * Te^3 / (12 * (1 - nu^2))
        
        Args:
            Te: Effective elastic thickness (km)
            E: Young's modulus (Pa), default 70 GPa
            nu: Poisson's ratio, default 0.25
            
        Returns:
            Flexural rigidity D (Pa * m^3)
        """
        Te_m = Te * 1000.0  # Convert km to m
        return E * Te_m**3 / (12.0 * (1.0 - nu**2))
    
    @staticmethod
    def compute_flexural_wavelength(D: float, rho_m: float = 3300.0, 
                                   rho_c: float = 2800.0, g: float = 9.81) -> float:
        """
        Compute flexural wavelength
        
        alpha = (4*pi*D / (g * (rho_m - rho_c)))^(1/4)
        wavelength = 2*pi*alpha
        
        Args:
            D: Flexural rigidity (Pa * m^3)
            rho_m: Mantle density (kg/m^3)
            rho_c: Crustal density (kg/m^3)
            g: Gravitational acceleration (m/s^2)
            
        Returns:
            Flexural wavelength (km)
        """
        alpha = (4.0 * np.pi * D / (g * (rho_m - rho_c)))**0.25
        wavelength = 2.0 * np.pi * alpha / 1000.0  # Convert to km
        return wavelength
    
    @staticmethod
    def compute_load_from_topography(topography: GridData, rho_c: float = 2800.0, 
                                     g: float = 9.81) -> GridData:
        """
        Compute load from topography
        
        q = rho_c * g * h
        
        Args:
            topography: GridData containing topography (m)
            rho_c: Crustal density (kg/m^3)
            g: Gravitational acceleration (m/s^2)
            
        Returns:
            GridData containing load (Pa)
        """
        load = GridData()
        load.nx = topography.nx
        load.ny = topography.ny
        load.dx = topography.dx
        load.dy = topography.dy
        load.xmin = topography.xmin
        load.xmax = topography.xmax
        load.ymin = topography.ymin
        load.ymax = topography.ymax
        
        # q = rho_c * g * h
        load.data = rho_c * g * topography.data
        load.update_range()
        
        return load
    
    @staticmethod
    def compute_response_function(Te: float, nx: int, ny: int, dx: float, dy: float,
                                  rho_m: float = 3300.0, rho_c: float = 2800.0,
                                  E: float = 70e9, nu: float = 0.25) -> GridData:
        """
        Compute the flexural response function (Green's function) for a given elastic thickness
        
        Args:
            Te: Effective elastic thickness (km)
            nx, ny: Grid dimensions
            dx, dy: Grid spacing (km)
            rho_m: Mantle density (kg/m^3)
            rho_c: Crustal density (kg/m^3)
            E: Young's modulus (Pa)
            nu: Poisson's ratio
            
        Returns:
            GridData containing the response function
        """
        D = FlexureResponse.compute_flexural_rigidity(Te, E, nu)
        return FlexureResponse._compute_response_spectral(D, nx, ny, dx, dy, rho_m, rho_c)
    
    @staticmethod
    def _compute_response_spectral(D: float, nx: int, ny: int, dx: float, dy: float,
                                   rho_m: float, rho_c: float, g: float = 9.81) -> GridData:
        """
        Compute 2D flexural response using spectral method
        
        Uses the analytical solution for a point load on an infinite elastic plate:
        w(r) = (load / (2*pi*D*alpha^2)) * kei(r/alpha)
        where kei is the Kelvin kei function
        """
        response = GridData()
        response.nx = nx
        response.ny = ny
        response.dx = dx
        response.dy = dy
        response.xmin = -(nx / 2) * dx
        response.xmax = (nx / 2) * dx
        response.ymin = -(ny / 2) * dy
        response.ymax = (ny / 2) * dy
        
        # Flexural parameter (characteristic length) in km
        alpha = (4.0 * D / (g * (rho_m - rho_c)))**0.25 / 1000.0
        
        # Scale factor
        scale = 1.0 / (2.0 * np.pi * D * (alpha * 1000.0)**2)  # Convert alpha to m
        
        # Create coordinate arrays
        y_coords = np.linspace(response.ymin, response.ymax - response.dy, ny)
        x_coords = np.linspace(response.xmin, response.xmax - response.dx, nx)
        X, Y = np.meshgrid(x_coords, y_coords)
        
        # Compute distance from center
        R = np.sqrt(X**2 + Y**2)
        
        # Compute response using Kelvin kei function
        R_alpha = R / alpha
        R_alpha_safe = np.where(R_alpha < 1e-10, 1e-10, R_alpha)
        
        # Use approximation for Kelvin kei function
        # For large r, use asymptotic approximation
        kei_values = np.zeros_like(R_alpha_safe)
        
        # For small r, use series expansion
        mask_small = R_alpha_safe < 2.0
        r_small = R_alpha_safe[mask_small]
        if np.any(mask_small):
            # Series expansion for kei(x) ≈ -pi/4 + x^2/4 * ln(x/2) - x^2/4 + ...
            r2 = r_small**2
            r4 = r2**2
            log_term = np.log(r_small / 2.0)
            kei_values[mask_small] = -np.pi / 4.0 + log_term * (1.0 - r2 / 4.0 + r4 / 64.0) - r2 / 4.0 * (1.0 - r2 / 8.0)
        
        # For large r, use asymptotic expansion
        mask_large = ~mask_small
        r_large = R_alpha_safe[mask_large]
        if np.any(mask_large):
            sqrt_r = np.sqrt(r_large)
            exp_term = np.exp(-r_large / np.sqrt(2.0))
            cos_term = np.cos(r_large / np.sqrt(2.0) - np.pi / 8.0)
            kei_values[mask_large] = -np.sqrt(np.pi / (2.0 * r_large)) * exp_term * cos_term
        
        # At center (r=0), use limiting value
        mask_zero = R < 1e-10
        if np.any(mask_zero):
            kei_values[mask_zero] = -np.pi / 4.0
        
        response.data = scale * kei_values
        response.update_range()
        
        return response

