"""
Data Loader for Surfer Grid Files (.grd)

This module provides functions to read Surfer ASCII grid files (DSAA format)
containing topography and Moho depth data.
"""

import numpy as np


def read_surfer_grd(filepath):
    """
    Read Surfer ASCII grid file (.grd) in DSAA format

    Parameters:
    -----------
    filepath : str
        Path to the .grd file

    Returns:
    --------
    X : 2D array
        X coordinate grid (meters)
    Y : 2D array
        Y coordinate grid (meters)
    data : 2D array
        Data values (same units as in file)
    dx : float
        Grid spacing in X direction (meters)
    dy : float
        Grid spacing in Y direction (meters)
    nx : int
        Number of columns
    ny : int
        Number of rows
    xmin, xmax : float
        X coordinate range
    ymin, ymax : float
        Y coordinate range
    """
    with open(filepath, "r") as f:
        # Read header
        header = f.readline().strip()
        if header != "DSAA":
            raise ValueError(f"File {filepath} is not in DSAA format. Found: {header}")

        # Read grid dimensions
        nx, ny = map(int, f.readline().split())

        # Read coordinate ranges
        xmin, xmax = map(float, f.readline().split())
        ymin, ymax = map(float, f.readline().split())
        zmin, zmax = map(float, f.readline().split())

        # Read all data values
        data_values = []
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                data_values.extend(map(float, line.split()))

        # Reshape data array (Surfer stores data row by row, starting from top)
        data = np.array(data_values).reshape(ny, nx)

        # Create coordinate grids
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        X, Y = np.meshgrid(x, y)

        # Calculate grid spacing
        dx = (xmax - xmin) / (nx - 1) if nx > 1 else 0
        dy = (ymax - ymin) / (ny - 1) if ny > 1 else 0

        print(f"Loaded {filepath}:")
        print(f"  Grid size: {nx} x {ny}")
        print(f"  X range: {xmin / 1000:.1f} to {xmax / 1000:.1f} km")
        print(f"  Y range: {ymin / 1000:.1f} to {ymax / 1000:.1f} km")
        print(f"  Grid spacing: {dx / 1000:.2f} x {dy / 1000:.2f} km")
        print(f"  Data range: {np.min(data):.2f} to {np.max(data):.2f}")
        print(f"  Data mean: {np.mean(data):.2f}, std: {np.std(data):.2f}")

        return X, Y, data, dx, dy, nx, ny, xmin, xmax, ymin, ymax


# Surfer no-data value (used when writing NaN so Surfer can blank cells)
SURFER_BLANK = 1.70141e38


def write_surfer_grd(filepath, X, Y, data, blank_nan=True):
    """
    Write a 2D grid to a Surfer ASCII grid file (.grd) in DSAA format.

    Parameters
    ----------
    filepath : str
        Path to the output .grd file.
    X : 2D array or 1D array
        X coordinates. If 2D (from meshgrid), first row is used; if 1D, used as-is.
    Y : 2D array or 1D array
        Y coordinates. If 2D (from meshgrid), first column is used; if 1D, used as-is.
    data : 2D array
        Data values, shape (ny, nx). Row 0 = smallest Y (consistent with read_surfer_grd).
    blank_nan : bool, optional
        If True (default), replace NaN with Surfer blank value (1.70141e38) so Surfer can blank cells.
    """
    # Get 1D axes from X, Y
    if np.ndim(X) == 2:
        x = X[0, :]
        y = Y[:, 0]
    else:
        x = np.asarray(X).ravel()
        y = np.asarray(Y).ravel()

    nx, ny = len(x), len(y)
    if data.shape != (ny, nx):
        raise ValueError(f"data shape {data.shape} does not match grid (ny={ny}, nx={nx})")

    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    z = np.asarray(data, dtype=float).copy()
    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if not np.isfinite(zmin):
        zmin, zmax = 0.0, 0.0
    if blank_nan and np.any(np.isnan(z)):
        z = np.where(np.isnan(z), SURFER_BLANK, z)

    with open(filepath, "w") as f:
        f.write("DSAA\n")
        f.write(f"{nx} {ny}\n")
        f.write(f"{xmin} {xmax}\n")
        f.write(f"{ymin} {ymax}\n")
        f.write(f"{zmin} {zmax}\n")
        # Write data row by row (first row = smallest y)
        for i in range(ny):
            line = " ".join(f"{z[i, j]:.6g}" for j in range(nx))
            f.write(line + "\n")


def check_grid_compatibility(X1, Y1, X2, Y2):
    """
    Check if two grids are compatible (same dimensions and coordinates)

    Parameters:
    -----------
    X1, Y1 : 2D arrays
        First grid coordinates
    X2, Y2 : 2D arrays
        Second grid coordinates

    Returns:
    --------
    compatible : bool
        True if grids are compatible
    message : str
        Description of compatibility status
    """
    if X1.shape != X2.shape or Y1.shape != Y2.shape:
        return False, f"Grid dimensions don't match: {X1.shape} vs {X2.shape}"

    if not np.allclose(X1, X2, rtol=1e-5):
        return False, "X coordinates don't match"

    if not np.allclose(Y1, Y2, rtol=1e-5):
        return False, "Y coordinates don't match"

    return True, "Grids are compatible"


def apply_taper(data, alpha=0.1):
    """
    Apply 2D cosine taper to reduce edge effects in FFT

    Parameters:
    -----------
    data : 2D array
        Input data
    alpha : float
        Taper parameter (0-1), fraction of data to taper at edges

    Returns:
    --------
    tapered_data : 2D array
        Tapered data
    """
    from scipy.signal.windows import tukey

    ny, nx = data.shape
    taper_x = tukey(nx, alpha)
    taper_y = tukey(ny, alpha)
    taper_2d = np.outer(taper_y, taper_x)

    return data * taper_2d


def prepare_data_for_inversion(
    topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1, ref_moho=None
):
    """
    Prepare topography and Moho data for inversion by:
    1. Removing mean (working with anomalies/undulations)
    2. Optionally applying taper to reduce edge effects

    Parameters:
    -----------
    topography : 2D array
        Topography data (m)
    moho_depth : 2D array
        Absolute Moho depth data (m)
    apply_taper_flag : bool
        Whether to apply taper (default: True)
    taper_alpha : float
        Taper parameter (default: 0.1)
    ref_moho : float, optional
        Reference Moho depth (m). If None, uses mean of moho_depth.

    Returns:
    --------
    topo_anom : 2D array
        Topography anomaly (m)
    moho_undulation : 2D array
        Moho undulation relative to reference (m)
    stats : dict
        Statistics about the data
    """
    # Remove mean to work with anomalies
    topo_mean = np.mean(topography)
    
    # Use provided reference or calculate mean
    if ref_moho is not None:
        moho_mean = ref_moho
    else:
        moho_mean = np.mean(moho_depth)

    topo_anom = topography - topo_mean
    moho_undulation = moho_depth - moho_mean

    # Apply taper if requested
    if apply_taper_flag:
        topo_anom = apply_taper(topo_anom, alpha=taper_alpha)
        moho_undulation = apply_taper(moho_undulation, alpha=taper_alpha)

    # Calculate statistics
    stats = {
        "topo_mean": topo_mean,
        "topo_std": np.std(topography),
        "topo_anom_std": np.std(topo_anom),
        "moho_mean": moho_mean,
        "moho_std": np.std(moho_depth),
        "moho_undulation_std": np.std(moho_undulation),
        "tapered": apply_taper_flag,
    }

    return topo_anom, moho_undulation, stats
