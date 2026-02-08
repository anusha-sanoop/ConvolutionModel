"""
Inverse Modeling of Elastic Thickness
Implements the Braitenberg et al. (2002) convolution method
"""

import numpy as np
from typing import Dict, List, Optional
from dataclasses import dataclass, field
from grd_reader import GridData, GRDReader
from flexure_response import FlexureResponse
from convolution import Convolution


@dataclass
class InverseParams:
    """Structure to hold inverse modeling parameters"""
    Te_min: float = 5.0           # Minimum Te to search (km)
    Te_max: float = 50.0          # Maximum Te to search (km)
    Te_step: float = 1.0          # Step size for Te search (km)
    rho_m: float = 3300.0         # Mantle density (kg/m^3)
    rho_c: float = 2800.0         # Crustal density (kg/m^3)
    E: float = 70e9               # Young's modulus (Pa)
    nu: float = 0.25              # Poisson's ratio
    g: float = 9.81               # Gravitational acceleration (m/s^2)
    max_iterations: int = 100     # Maximum iterations for optimization
    tolerance: float = 0.001      # Convergence tolerance
    # Window-based analysis parameters
    window_size: Optional[float] = None     # Window size in meters (for sliding window analysis)
    shift_min: Optional[float] = None       # Minimum shift in meters
    shift_max: Optional[float] = None       # Maximum shift in meters
    shift_step: Optional[float] = None      # Shift step in meters
    use_sliding_window: bool = False  # Enable sliding window analysis


@dataclass
class InverseResult:
    """Structure to hold inverse modeling results"""
    Te_map: GridData              # Spatial map of estimated Te
    modeled_deflection: GridData  # Modeled deflection
    misfit_map: GridData          # Spatial misfit map
    best_misfit: float            # Best overall misfit
    best_correlation: float       # Best correlation coefficient
    best_Te: float                # Best Te value
    Te_values: List[float]        # Te values tested
    misfits: List[float]          # Corresponding misfits


class InverseModel:
    """Class for inverse modeling of elastic thickness"""
    
    def __init__(self):
        """Initialize the inverse modeler"""
        pass
    
    @staticmethod
    def get_default_eastern_alps_params() -> InverseParams:
        """Get default parameters for Eastern Alps (from paper)"""
        return InverseParams(
            Te_min=5.0,
            Te_max=50.0,
            Te_step=1.0,
            rho_m=3300.0,
            rho_c=2800.0,
            E=70e9,
            nu=0.25,
            g=9.81,
            max_iterations=100,
            tolerance=0.001,
            use_sliding_window=False
        )
    
    @staticmethod
    def get_mars_params(window_size: float = 1000000.0, shift_min: float = 20000.0,
                       shift_max: float = 80000.0, shift_step: float = 20000.0,
                       Te_min: float = 5.0, Te_max: float = 50.0, Te_step: float = 1.0) -> InverseParams:
        """
        Get parameters for Mars analysis
        
        Args:
            window_size: Window size in meters (default: 1000 km = 1,000,000 m)
            shift_min: Minimum shift in meters (default: 20 km = 20,000 m)
            shift_max: Maximum shift in meters (default: 80 km = 80,000 m)
            shift_step: Shift step in meters (default: 20 km = 20,000 m)
            Te_min: Minimum Te to search (km)
            Te_max: Maximum Te to search (km)
            Te_step: Te step size (km)
        """
        return InverseParams(
            Te_min=Te_min,
            Te_max=Te_max,
            Te_step=Te_step,
            rho_m=3500.0,          # Mars mantle density
            rho_c=2900.0,          # Mars crustal density
            E=70e9,
            nu=0.25,
            g=3.72,                # Mars gravity
            max_iterations=100,
            tolerance=0.001,
            window_size=window_size,
            shift_min=shift_min,
            shift_max=shift_max,
            shift_step=shift_step,
            use_sliding_window=True
        )
    
    def forward_model(self, topography: GridData, Te: float, 
                     params: InverseParams) -> GridData:
        """
        Perform forward modeling for a given Te
        
        Args:
            topography: GridData of topography
            Te: Effective elastic thickness (km)
            params: Model parameters
            
        Returns:
            GridData of modeled deflection/Moho depth
        """
        # 1. Compute load from topography
        load = FlexureResponse.compute_load_from_topography(
            topography, params.rho_c, params.g
        )
        
        # 2. Compute flexural response function
        # Convert dx, dy from meters to km for response function calculation
        dx_km = topography.dx / 1000.0  # Convert m to km
        dy_km = topography.dy / 1000.0  # Convert m to km
        
        response = FlexureResponse.compute_response_function(
            Te, topography.nx, topography.ny,
            dx_km, dy_km,  # Response function expects km
            params.rho_m, params.rho_c,
            params.E, params.nu
        )
        
        # 3. Convolve load with response to get deflection
        deflection = Convolution.convolve(load, response, method='fft')
        
        return deflection
    
    def estimate_te(self, topography: GridData, observed_moho: GridData,
                   params: InverseParams) -> InverseResult:
        """
        Perform inverse modeling to estimate Te from observed Moho depth
        
        Args:
            topography: GridData of topography
            observed_moho: GridData of observed Moho depth
            params: Inverse modeling parameters
            
        Returns:
            InverseResult containing the estimated Te and modeling results
        """
        result = InverseResult(
            Te_map=GridData(),
            modeled_deflection=GridData(),
            misfit_map=GridData(),
            best_misfit=float('inf'),
            best_correlation=0.0,
            best_Te=params.Te_min,
            Te_values=[],
            misfits=[]
        )
        
        # Check grid compatibility
        if topography.nx != observed_moho.nx or topography.ny != observed_moho.ny:
            raise ValueError(
                f"Topography and Moho grids must have same dimensions. "
                f"Got {topography.nx}x{topography.ny} and {observed_moho.nx}x{observed_moho.ny}"
            )
        
        # Initialize result grids
        result.Te_map.nx = topography.nx
        result.Te_map.ny = topography.ny
        result.Te_map.dx = topography.dx
        result.Te_map.dy = topography.dy
        result.Te_map.xmin = topography.xmin
        result.Te_map.xmax = topography.xmax
        result.Te_map.ymin = topography.ymin
        result.Te_map.ymax = topography.ymax
        
        result.misfit_map = GridData()
        result.misfit_map.nx = topography.nx
        result.misfit_map.ny = topography.ny
        result.misfit_map.dx = topography.dx
        result.misfit_map.dy = topography.dy
        result.misfit_map.xmin = topography.xmin
        result.misfit_map.xmax = topography.xmax
        result.misfit_map.ymin = topography.ymin
        result.misfit_map.ymax = topography.ymax
        
        print(f"Starting inverse modeling...")
        print(f"Grid size: {topography.nx} x {topography.ny}")
        print(f"Testing Te values from {params.Te_min} to {params.Te_max} km (step: {params.Te_step} km)...")
        
        # Grid search for best global Te
        Te_values = np.arange(params.Te_min, params.Te_max + params.Te_step, params.Te_step)
        
        for Te in Te_values:
            result.Te_values.append(Te)
            
            # Forward model
            modeled = self.forward_model(topography, Te, params)
            
            # Compute misfit
            misfit = Convolution.compute_misfit(observed_moho, modeled)
            correlation = Convolution.compute_correlation(observed_moho, modeled)
            
            result.misfits.append(misfit)
            
            print(f"  Te = {Te:6.1f} km: Misfit = {misfit:10.3f}, Correlation = {correlation:6.3f}")
            
            if misfit < result.best_misfit:
                result.best_misfit = misfit
                result.best_Te = Te
                result.modeled_deflection = modeled
                result.best_correlation = correlation
        
        # Create Te map (for now, use constant Te - can be extended to spatially varying)
        result.Te_map.data = np.ones((topography.ny, topography.nx)) * result.best_Te
        
        # Create misfit map
        diff = observed_moho.data - result.modeled_deflection.data
        result.misfit_map.data = diff**2
        
        result.Te_map.update_range()
        result.misfit_map.update_range()
        
        print(f"\nBest Te = {result.best_Te:.2f} km")
        print(f"Best misfit (RMS) = {result.best_misfit:.3f}")
        print(f"Best correlation = {result.best_correlation:.3f}")
        
        return result
    
    def estimate_te_local(self, topography: GridData, observed_moho: GridData,
                         params: InverseParams) -> float:
        """
        Estimate Te at a single location using grid search
        
        Args:
            topography: Local topography
            observed_moho: Observed Moho depth
            params: Inverse modeling parameters
            
        Returns:
            Best Te value
        """
        best_Te = params.Te_min
        best_misfit = float('inf')
        
        Te_values = np.arange(params.Te_min, params.Te_max + params.Te_step, params.Te_step)
        
        # Debug: check if forward model produces different results
        misfits_list = []
        
        for Te in Te_values:
            # Forward model
            try:
                modeled = self.forward_model(topography, Te, params)
                
                # Compute misfit
                misfit = Convolution.compute_misfit(observed_moho, modeled)
                misfits_list.append(misfit)
                
                if misfit < best_misfit:
                    best_misfit = misfit
                    best_Te = Te
            except Exception as e:
                # Skip this Te value if forward modeling fails
                misfits_list.append(float('inf'))
                continue
        
        # Check if we found a valid solution
        if best_misfit == float('inf') or len(misfits_list) == 0:
            # Return middle value as fallback
            best_Te = (params.Te_min + params.Te_max) / 2.0
        
        return best_Te
    
    def _extract_window(self, grid: GridData, center_x: float, center_y: float, 
                       window_size: float) -> GridData:
        """
        Extract a window around a center point
        
        Args:
            grid: Input grid
            center_x: X coordinate of window center (in same units as grid coordinates)
            center_y: Y coordinate of window center (in same units as grid coordinates)
            window_size: Size of window (in meters)
            
        Returns:
            GridData containing the window
        """
        # Validate grid data
        if grid.data is None:
            raise ValueError("Grid data is None - grid not properly loaded")
        
        if grid.nx == 0 or grid.ny == 0:
            raise ValueError(f"Grid has invalid dimensions: {grid.nx} x {grid.ny}")
        
        # Ensure grid.data is 2D
        if grid.data.ndim == 1:
            raise ValueError(f"Grid data is 1D but expected 2D. Shape: {grid.data.shape}, Expected: ({grid.ny}, {grid.nx})")
        
        if grid.data.shape != (grid.ny, grid.nx):
            raise ValueError(f"Grid data shape mismatch: data shape {grid.data.shape} != expected ({grid.ny}, {grid.nx})")
        
        # Validate grid spacing
        if grid.dx <= 0:
            raise ValueError(f"Invalid grid spacing dx: {grid.dx}")
        if grid.dy <= 0:
            raise ValueError(f"Invalid grid spacing dy: {grid.dy}")
        
        # Convert window size from meters to grid coordinate units
        # Assume grid coordinates are in meters (typical for Surfer GRD files)
        # window_size is in meters, grid.dx is in meters, so no conversion needed
        window_half = window_size / (2.0 * grid.dx)  # window_size and dx both in meters
        
        # Find center indices
        center_i = int((center_y - grid.ymin) / grid.dy)
        center_j = int((center_x - grid.xmin) / grid.dx)
        
        # Calculate window bounds in grid indices
        i_start = max(0, int(center_i - window_half))
        i_end = min(grid.ny, int(center_i + window_half))
        j_start = max(0, int(center_j - window_half))
        j_end = min(grid.nx, int(center_j + window_half))
        
        # Ensure we have valid window bounds
        if i_end <= i_start or j_end <= j_start:
            raise ValueError(f"Invalid window bounds: i=[{i_start}:{i_end}], j=[{j_start}:{j_end}]")
        
        # Extract window
        window = GridData()
        window.nx = j_end - j_start
        window.ny = i_end - i_start
        window.dx = grid.dx
        window.dy = grid.dy
        window.xmin = grid.xmin + j_start * grid.dx
        window.ymin = grid.ymin + i_start * grid.dy
        window.xmax = window.xmin + window.nx * window.dx
        window.ymax = window.ymin + window.ny * window.dy
        
        # Extract data with bounds checking
        try:
            window.data = grid.data[i_start:i_end, j_start:j_end].copy()
        except IndexError as e:
            raise IndexError(
                f"Failed to extract window. Grid shape: {grid.data.shape}, "
                f"Requested: [{i_start}:{i_end}, {j_start}:{j_end}]. "
                f"Grid bounds: ny={grid.ny}, nx={grid.nx}. "
                f"Original error: {str(e)}"
            )
        
        window.update_range()
        
        return window
    
    def estimate_te_spatial_variable(self, topography: GridData, observed_moho: GridData,
                                    params: InverseParams) -> InverseResult:
        """
        Estimate spatially variable Te using sliding window approach
        
        Args:
            topography: Topography grid
            observed_moho: Observed Moho grid
            params: Inverse modeling parameters (must have window_size, shift_* set)
            
        Returns:
            InverseResult with spatially variable Te map
        """
        if not params.use_sliding_window or params.window_size is None:
            raise ValueError("Sliding window parameters not set. Use use_sliding_window=True and set window_size.")
        
        result = InverseResult(
            Te_map=GridData(),
            modeled_deflection=GridData(),
            misfit_map=GridData(),
            best_misfit=float('inf'),
            best_correlation=0.0,
            best_Te=params.Te_min,
            Te_values=[],
            misfits=[]
        )
        
        # Initialize result grids
        result.Te_map.nx = topography.nx
        result.Te_map.ny = topography.ny
        result.Te_map.dx = topography.dx
        result.Te_map.dy = topography.dy
        result.Te_map.xmin = topography.xmin
        result.Te_map.xmax = topography.xmax
        result.Te_map.ymin = topography.ymin
        result.Te_map.ymax = topography.ymax
        result.Te_map.data = np.full((topography.ny, topography.nx), np.nan)
        
        result.misfit_map = GridData()
        result.misfit_map.nx = topography.nx
        result.misfit_map.ny = topography.ny
        result.misfit_map.dx = topography.dx
        result.misfit_map.dy = topography.dy
        result.misfit_map.xmin = topography.xmin
        result.misfit_map.xmax = topography.xmax
        result.misfit_map.ymin = topography.ymin
        result.misfit_map.ymax = topography.ymax
        result.misfit_map.data = np.full((topography.ny, topography.nx), np.nan)
        
        result.modeled_deflection = GridData()
        result.modeled_deflection.nx = topography.nx
        result.modeled_deflection.ny = topography.ny
        result.modeled_deflection.dx = topography.dx
        result.modeled_deflection.dy = topography.dy
        result.modeled_deflection.xmin = topography.xmin
        result.modeled_deflection.xmax = topography.xmax
        result.modeled_deflection.ymin = topography.ymin
        result.modeled_deflection.ymax = topography.ymax
        result.modeled_deflection.data = np.zeros((topography.ny, topography.nx))
        
        # Convert shift parameters from meters to grid units
        # Assuming grid coordinates are in meters (typical for Surfer GRD files)
        shift_min_grid = params.shift_min / topography.dx if params.shift_min else 0
        shift_max_grid = params.shift_max / topography.dx if params.shift_max else 0
        shift_step_grid = params.shift_step / topography.dx if params.shift_step else 1
        
        window_half_grid = params.window_size / (2.0 * topography.dx)
        
        # Validate input grids
        if topography.data is None or observed_moho.data is None:
            raise ValueError("Input grids have no data - check that GRD files were read correctly")
        
        if topography.nx == 0 or topography.ny == 0:
            raise ValueError(f"Topography grid has invalid dimensions: {topography.nx} x {topography.ny}")
        
        if observed_moho.nx == 0 or observed_moho.ny == 0:
            raise ValueError(f"Observed Moho grid has invalid dimensions: {observed_moho.nx} x {observed_moho.ny}")
        
        if topography.data.ndim != 2:
            raise ValueError(f"Topography data must be 2D, got {topography.data.ndim}D with shape {topography.data.shape}")
        
        if observed_moho.data.ndim != 2:
            raise ValueError(f"Observed Moho data must be 2D, got {observed_moho.data.ndim}D with shape {observed_moho.data.shape}")
        
        print(f"Starting spatially variable Te estimation...")
        print(f"  Window size: {params.window_size / 1000.0:.1f} km ({params.window_size:.0f} m)")
        print(f"  Grid size: {topography.nx} x {topography.ny}")
        print(f"  Grid spacing: dx={topography.dx:.2f} m, dy={topography.dy:.2f} m")
        print(f"  Grid extent: X=[{topography.xmin:.0f}, {topography.xmax:.0f}] m, Y=[{topography.ymin:.0f}, {topography.ymax:.0f}] m")
        print(f"  Topography data shape: {topography.data.shape}")
        print(f"  Observed Moho data shape: {observed_moho.data.shape}")
        print(f"  Window half-size in grid units: {window_half_grid:.1f}")
        
        # Validate window size is reasonable
        if window_half_grid <= 0:
            raise ValueError(
                f"Window half-size is {window_half_grid:.2f} grid units. "
                f"This suggests a unit mismatch. Window size: {params.window_size} m, "
                f"Grid spacing dx: {topography.dx} m"
            )
        
        if window_half_grid * 2 > min(topography.nx, topography.ny):
            raise ValueError(
                f"Window size ({window_half_grid * 2:.1f} grid units) is larger than "
                f"grid dimensions ({topography.nx} x {topography.ny})"
            )
        
        # Generate shift positions
        # If shift_step is specified, use it as the spacing between window centers
        # Otherwise, use window size as spacing (no overlap)
        if params.shift_step is not None and params.shift_step > 0:
            # Use shift_step as spacing between window centers
            step_grid = params.shift_step / topography.dx
            # Start from window edge (not center), end before grid edge
            start = window_half_grid
            end_x = max(start, topography.nx - window_half_grid)
            end_y = max(start, topography.ny - window_half_grid)
            shift_x_values = np.arange(start, end_x, step_grid)
            shift_y_values = np.arange(start, end_y, step_grid)
        elif params.shift_min is not None and params.shift_max is not None:
            # Use specified shift range
            shift_x_values = np.arange(shift_min_grid, shift_max_grid + shift_step_grid, shift_step_grid)
            shift_y_values = np.arange(shift_min_grid, shift_max_grid + shift_step_grid, shift_step_grid)
            
            # Ensure shifts are reasonable (don't exceed grid bounds)
            max_shift_x = max(0, topography.nx - window_half_grid * 2 - 1)
            max_shift_y = max(0, topography.ny - window_half_grid * 2 - 1)
            shift_x_values = shift_x_values[(shift_x_values >= window_half_grid) & (shift_x_values < max_shift_x)]
            shift_y_values = shift_y_values[(shift_y_values >= window_half_grid) & (shift_y_values < max_shift_y)]
        else:
            # Default: analyze at regular intervals (every window_size, no overlap)
            step = window_half_grid * 2
            if step <= 0:
                step = min(topography.nx, topography.ny) / 10  # Fallback: 10 windows across
            shift_x_values = np.arange(window_half_grid, topography.nx - window_half_grid, step)
            shift_y_values = np.arange(window_half_grid, topography.ny - window_half_grid, step)
        
        # Ensure we have valid shift values
        if len(shift_x_values) == 0:
            shift_x_values = np.array([topography.nx / 2])  # Center point
        if len(shift_y_values) == 0:
            shift_y_values = np.array([topography.ny / 2])  # Center point
        
        print(f"  Number of analysis windows: {len(shift_x_values)} x {len(shift_y_values)} = {len(shift_x_values) * len(shift_y_values)}")
        
        total_windows = len(shift_x_values) * len(shift_y_values)
        window_count = 0
        
        # Process each window position
        for shift_y in shift_y_values:
            for shift_x in shift_x_values:
                window_count += 1
                # shift_x and shift_y are grid indices, convert to coordinates
                center_x = topography.xmin + shift_x * topography.dx
                center_y = topography.ymin + shift_y * topography.dy
                
                # Extract windows around this center point
                try:
                    topo_window = self._extract_window(topography, center_x, center_y, params.window_size)
                    moho_window = self._extract_window(observed_moho, center_x, center_y, params.window_size)
                except (ValueError, IndexError) as e:
                    print(f"  Warning: Skipping window at ({center_x:.1f}, {center_y:.1f}): {str(e)}")
                    continue
                
                if topo_window.nx < 10 or topo_window.ny < 10:
                    continue  # Skip windows that are too small
                
                # Estimate Te for this window
                try:
                    local_Te = self.estimate_te_local(topo_window, moho_window, params)
                    
                    # Find indices in original grid corresponding to window center
                    center_i = int((center_y - topography.ymin) / topography.dy)
                    center_j = int((center_x - topography.xmin) / topography.dx)
                    
                    if 0 <= center_i < topography.ny and 0 <= center_j < topography.nx:
                        result.Te_map.data[center_i, center_j] = local_Te
                        
                        # Forward model for this Te to get deflection
                        modeled_window = self.forward_model(topo_window, local_Te, params)
                        
                        # Compute local misfit (RMS)
                        local_misfit = Convolution.compute_misfit(moho_window, modeled_window)
                        result.misfit_map.data[center_i, center_j] = local_misfit * local_misfit  # Store squared for consistency
                
                except Exception as e:
                    print(f"  Warning: Failed to process window at ({center_x:.1f}, {center_y:.1f}): {e}")
                    continue
                
                if window_count % 10 == 0:
                    print(f"  Processed {window_count}/{total_windows} windows...")
        
        print(f"Completed processing {window_count} windows.")
        
        # Interpolate Te map to fill gaps (simple approach: nearest neighbor)
        # For production, use more sophisticated interpolation
        try:
            from scipy.interpolate import griddata
            use_scipy = True
        except ImportError:
            use_scipy = False
            print("  Warning: scipy.interpolate not available, using simple nearest neighbor fill")
        
        # Find valid points and fill gaps
        valid_mask = ~np.isnan(result.Te_map.data)
        if np.any(valid_mask):
            if use_scipy:
                y_coords, x_coords = np.mgrid[
                    topography.ymin:topography.ymax:topography.ny*1j,
                    topography.xmin:topography.xmax:topography.nx*1j
                ]
                
                valid_y, valid_x = np.where(valid_mask)
                valid_points = np.column_stack([
                    x_coords[valid_y, valid_x],
                    y_coords[valid_y, valid_x]
                ])
                valid_values = result.Te_map.data[valid_mask]
                
                # Interpolate to full grid
                points = np.column_stack([x_coords.ravel(), y_coords.ravel()])
                interpolated = griddata(valid_points, valid_values, points, method='nearest')
                result.Te_map.data = interpolated.reshape((topography.ny, topography.nx))
            else:
                # Simple fill: use mean of valid values for NaN regions
                nan_mask = np.isnan(result.Te_map.data)
                if np.any(nan_mask):
                    mean_Te = np.nanmean(result.Te_map.data)
                    result.Te_map.data[nan_mask] = mean_Te
        
        result.Te_map.update_range()
        
        # Compute global forward model using interpolated Te map (simplified)
        # For each grid point, use local Te to compute deflection
        # This is computationally expensive, so we do a simplified version
        
        result.best_Te = np.nanmean(result.Te_map.data)
        result.best_misfit = np.nanmean(result.misfit_map.data) if np.any(~np.isnan(result.misfit_map.data)) else float('inf')
        
        return result

