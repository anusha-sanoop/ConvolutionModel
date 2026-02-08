"""
3D Data Loader for Volumetric Data

This module provides functions to read and process 3D volumetric data
for elastic thickness inversion. Supports multiple 2D slices or true 3D volumes.
"""

import numpy as np


def read_3d_volume_from_slices(file_pattern, z_values=None):
    """
    Read 3D volume from multiple 2D grid files (slices)
    
    Parameters:
    -----------
    file_pattern : str or list
        Pattern for file names (e.g., "data_z{}.grd") or list of file paths
    z_values : array, optional
        Z coordinate values for each slice (depth or elevation in meters)
        If None, will use indices
    
    Returns:
    --------
    X : 3D array
        X coordinate grid (meters)
    Y : 3D array
        Y coordinate grid (meters)
    Z : 3D array
        Z coordinate grid (meters)
    data : 3D array
        Data values (shape: nz, ny, nx)
    dx, dy, dz : float
        Grid spacing in each direction (meters)
    """
    import glob
    
    if isinstance(file_pattern, str):
        # Find all matching files
        files = sorted(glob.glob(file_pattern))
    else:
        files = sorted(file_pattern)
    
    if len(files) == 0:
        raise ValueError(f"No files found matching pattern: {file_pattern}")
    
    # Read first file to get dimensions
    from Code.data_loader import read_surfer_grd
    X_2d, Y_2d, first_slice, dx, dy, nx, ny, xmin, xmax, ymin, ymax = read_surfer_grd(files[0])
    
    nz = len(files)
    data = np.zeros((nz, ny, nx))
    
    # Read all slices
    for i, filepath in enumerate(files):
        _, _, slice_data, _, _, _, _, _, _, _, _ = read_surfer_grd(filepath)
        data[i, :, :] = slice_data
    
    # Create 3D coordinate grids
    if z_values is None:
        z_values = np.arange(nz) * 1000  # Default 1 km spacing
    
    Z_coords = np.zeros((nz, ny, nx))
    X_coords = np.zeros((nz, ny, nx))
    Y_coords = np.zeros((nz, ny, nx))
    
    for i, z_val in enumerate(z_values):
        Z_coords[i, :, :] = z_val
        X_coords[i, :, :] = X_2d
        Y_coords[i, :, :] = Y_2d
    
    dz = np.mean(np.diff(z_values)) if len(z_values) > 1 else 1000.0
    
    return X_coords, Y_coords, Z_coords, data, dx, dy, dz


def read_3d_surfer_grd(filepath):
    """
    Read 3D Surfer grid file (if format supports 3D)
    For now, this is a placeholder - most .grd files are 2D
    
    Parameters:
    -----------
    filepath : str
        Path to the 3D .grd file
    
    Returns:
    --------
    X, Y, Z : 3D arrays
        Coordinate grids
    data : 3D array
        Data values
    dx, dy, dz : float
        Grid spacing
    """
    # This would need to be implemented based on specific 3D format
    # For now, raise NotImplementedError
    raise NotImplementedError(
        "3D .grd format not yet implemented. "
        "Use read_3d_volume_from_slices() to load multiple 2D slices."
    )


def prepare_3d_data_for_inversion(
    topography_3d, moho_3d, apply_taper_flag=True, taper_alpha=0.1
):
    """
    Prepare 3D data for inversion by removing mean and applying taper
    
    Parameters:
    -----------
    topography_3d : 3D array
        Topography volume (nz, ny, nx)
    moho_3d : 3D array
        Moho depth volume (nz, ny, nx)
    apply_taper_flag : bool
        Whether to apply cosine taper to edges
    taper_alpha : float
        Taper width as fraction of domain (0-1)
    
    Returns:
    --------
    topo_anom_3d : 3D array
        Topography anomaly (zero-mean)
    moho_undulation_3d : 3D array
        Moho undulation (zero-mean)
    stats : dict
        Statistics dictionary
    """
    # Remove mean from each slice
    topo_anom_3d = topography_3d.copy()
    moho_undulation_3d = moho_3d.copy()
    
    nz, ny, nx = topography_3d.shape
    
    # Calculate mean for each depth level
    for z in range(nz):
        topo_anom_3d[z, :, :] = topo_anom_3d[z, :, :] - np.mean(topo_anom_3d[z, :, :])
        moho_undulation_3d[z, :, :] = moho_undulation_3d[z, :, :] - np.mean(moho_undulation_3d[z, :, :])
    
    # Apply taper if requested
    if apply_taper_flag:
        # Create 2D taper mask
        taper_x = int(nx * taper_alpha)
        taper_y = int(ny * taper_alpha)
        
        taper_mask = np.ones((ny, nx))
        
        # X-direction taper
        for i in range(taper_x):
            taper_val = 0.5 * (1 - np.cos(np.pi * i / taper_x))
            taper_mask[:, i] *= taper_val
            taper_mask[:, nx - 1 - i] *= taper_val
        
        # Y-direction taper
        for i in range(taper_y):
            taper_val = 0.5 * (1 - np.cos(np.pi * i / taper_y))
            taper_mask[i, :] *= taper_val
            taper_mask[ny - 1 - i, :] *= taper_val
        
        # Apply taper to all depth levels
        for z in range(nz):
            topo_anom_3d[z, :, :] *= taper_mask
            moho_undulation_3d[z, :, :] *= taper_mask
    
    # Calculate statistics
    stats = {
        "topo_mean": np.mean(topography_3d),
        "topo_std": np.std(topography_3d),
        "topo_anom_std": np.std(topo_anom_3d),
        "moho_mean": np.mean(moho_3d),
        "moho_std": np.std(moho_3d),
        "moho_undulation_std": np.std(moho_undulation_3d),
        "tapered": apply_taper_flag,
        "shape": (nz, ny, nx),
    }
    
    return topo_anom_3d, moho_undulation_3d, stats


def check_3d_grid_compatibility(X1, Y1, Z1, X2, Y2, Z2, tol=1.0):
    """
    Check if two 3D grids are compatible (same coordinates)
    
    Parameters:
    -----------
    X1, Y1, Z1 : 3D arrays
        First grid coordinates
    X2, Y2, Z2 : 3D arrays
        Second grid coordinates
    tol : float
        Tolerance in meters
    
    Returns:
    --------
    compatible : bool
        Whether grids are compatible
    message : str
        Status message
    """
    # Check shape compatibility
    if X1.shape != X2.shape:
        return False, f"Grid shapes differ: {X1.shape} vs {X2.shape}"
    
    # Check coordinate compatibility
    x_diff = np.max(np.abs(X1 - X2))
    y_diff = np.max(np.abs(Y1 - Y2))
    z_diff = np.max(np.abs(Z1 - Z2))
    
    if x_diff > tol or y_diff > tol or z_diff > tol:
        return False, f"Coordinate mismatch: X={x_diff:.1f}m, Y={y_diff:.1f}m, Z={z_diff:.1f}m"
    
    return True, "Grids are compatible"

