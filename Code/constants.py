"""
Physical Constants for Elastic Thickness Inversion

This module centralizes all physical constants used throughout the codebase.
Modify values here to update them across all modules.

Constants are defined for Mars by default, but can be adjusted for other planetary bodies.
"""


# DENSITY PARAMETERS (kg/m³)


# Load density (crustal density)
RHO_LOAD = 2900.0  # kg/m³ (Mars crustal density)

# Mantle density
RHO_MANTLE = 3500.0  # kg/m³ (Mars mantle density)

# Infill density
# - 0 for topography in air (above datum)
# - 2900 for subsurface loads or crustal compensation
# - ~1000 for submarine topography (water)
RHO_INFILL = 2900.0  # kg/m³ (default: crustal compensation)

# Density contrast (for gravity inversion)
# Typically: RHO_MANTLE - RHO_LOAD

#DENSITY_CONTRAST = 600.0  # kg/m³ (3500 - 2900 for Mars)
DENSITY_CONTRAST = RHO_MANTLE - RHO_LOAD



# ELASTIC PARAMETERS


# Young's modulus
YOUNGS_MODULUS = 1.0e11  # Pa (100 GPa)

# Poisson's ratio
POISSONS_RATIO = 0.25  # dimensionless



# GRAVITATIONAL PARAMETERS


# Gravitational acceleration
GRAVITY = 3.72  # m/s² (Mars surface gravity)

# Universal gravitational constant
GRAVITATIONAL_CONSTANT = 6.674e-11  # m³/kg/s²



# MOHO PARAMETERS


# Default reference Moho depth (for gravity inversion when user doesn't provide one)
# This is a fallback default, not a constant - should be provided by user or calculated from data
DEFAULT_REFERENCE_MOHO_DEPTH = 50000.0  # m (50 km default fallback)



# GRID PARAMETERS (defaults)


# Default grid spacing
DEFAULT_DX = 1000.0  # m
DEFAULT_DY = 1000.0  # m



# CONSTANT SETS FOR DIFFERENT SCENARIOS


# For topography in air (above datum)
CONSTANTS_AIR = {
    'rho_load': RHO_LOAD,
    'rho_m': RHO_MANTLE,
    'rho_infill': 0.0,  # Air
    'E': YOUNGS_MODULUS,
    'nu': POISSONS_RATIO,
    'g': GRAVITY,
}

# For subsurface loads or crustal compensation
CONSTANTS_SUBSURFACE = {
    'rho_load': RHO_LOAD,
    'rho_m': RHO_MANTLE,
    'rho_infill': RHO_INFILL,  # Crustal density
    'E': YOUNGS_MODULUS,
    'nu': POISSONS_RATIO,
    'g': GRAVITY,
}

# For gravity inversion
CONSTANTS_GRAVITY = {
    'density_contrast': DENSITY_CONTRAST,
    'reference_moho': DEFAULT_REFERENCE_MOHO_DEPTH,  # Fallback default only
    'G': GRAVITATIONAL_CONSTANT,
}
