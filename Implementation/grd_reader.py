"""
GRD File Reader/Writer for GMT grid format
Handles reading and writing of GRD (GMT grid) files
Supports both ASCII and NetCDF binary formats
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional

# Try to import netCDF4 for binary GRD support
try:
    import netCDF4
    HAS_NETCDF = True
except ImportError:
    HAS_NETCDF = False


class GridData:
    """Structure to hold GRD grid data"""
    def __init__(self):
        self.data: np.ndarray = None  # Grid values as numpy array
        self.nx: int = 0               # Number of columns
        self.ny: int = 0               # Number of rows
        self.xmin: float = 0.0         # Minimum X coordinate
        self.xmax: float = 0.0         # Maximum X coordinate
        self.ymin: float = 0.0         # Minimum Y coordinate
        self.ymax: float = 0.0         # Maximum Y coordinate
        self.dx: float = 0.0           # Grid spacing in X
        self.dy: float = 0.0           # Grid spacing in Y
        self.z_min: float = 0.0        # Minimum data value
        self.z_max: float = 0.0        # Maximum data value
    
    def update_range(self):
        """Update z_min and z_max from data"""
        if self.data is not None and self.data.size > 0:
            # Exclude NaN and invalid values
            valid_data = self.data[np.isfinite(self.data)]
            if valid_data.size > 0:
                self.z_min = float(np.min(valid_data))
                self.z_max = float(np.max(valid_data))


class GRDReader:
    """Class for reading and writing GMT GRD format files"""
    
    @staticmethod
    def file_exists(filename: str) -> bool:
        """Check if a file exists"""
        return Path(filename).exists()
    
    @staticmethod
    def _is_netcdf_file(filename: str) -> bool:
        """Check if file is likely a NetCDF file"""
        try:
            with open(filename, 'rb') as f:
                # Check for NetCDF magic number (CDF\x01 or CDF\x02)
                magic = f.read(4)
                return magic == b'CDF\x01' or magic == b'CDF\x02' or magic == b'\x89HDF'
        except:
            return False
    
    @staticmethod
    def _is_surfer_ascii(filename: str) -> bool:
        """Check if file is a Surfer ASCII GRD file (starts with DSAA)"""
        try:
            with open(filename, 'r') as f:
                first_line = f.readline().strip()
                return first_line == 'DSAA'
        except:
            return False
    
    @staticmethod
    def read_grd(filename: str) -> GridData:
        """
        Read a GRD file (ASCII or NetCDF binary format)
        
        Args:
            filename: Path to the GRD file
            
        Returns:
            GridData object containing the grid information
        """
        grid = GridData()
        
        if not GRDReader.file_exists(filename):
            raise FileNotFoundError(f"Cannot open file: {filename}")
        
        path = Path(filename)
        
        # Check if it's a Surfer ASCII GRD file first (most common)
        if GRDReader._is_surfer_ascii(str(path)):
            try:
                return GRDReader._read_surfer_ascii_grd(str(path))
            except Exception as e:
                raise IOError(f"Failed to read Surfer ASCII GRD file {filename}: {str(e)}")
        
        # Check if it's a NetCDF file
        if GRDReader._is_netcdf_file(str(path)) or (HAS_NETCDF and filename.endswith('.grd')):
            try:
                return GRDReader._read_netcdf_grd(str(path))
            except Exception as e:
                print(f"Warning: Failed to read as NetCDF ({str(e)}), trying GMT ASCII format...")
        
        # Try to read as GMT ASCII GRD
        try:
            with open(path, 'r') as f:
                header = {}
                line = f.readline()
                
                # Read header lines
                while line:
                    line = line.strip()
                    if not line:
                        line = f.readline()
                        continue
                    
                    # Parse header information
                    parts = line.split()
                    if len(parts) >= 2:
                        key = parts[0].lower()
                        value = parts[1]
                        
                        if key == 'ncols':
                            header['ncols'] = int(value)
                        elif key == 'nrows':
                            header['nrows'] = int(value)
                        elif key == 'xllcorner' or key == 'xllcenter':
                            header['xllcorner'] = float(value)
                        elif key == 'yllcorner' or key == 'yllcenter':
                            header['yllcorner'] = float(value)
                        elif key == 'cellsize':
                            header['cellsize'] = float(value)
                        elif key == 'nodata_value' or key == 'nodata':
                            header['nodata'] = float(value)
                            # Header ends here
                            break
                    
                    line = f.readline()
                
                # Read grid dimensions
                grid.nx = header.get('ncols', 0)
                grid.ny = header.get('nrows', 0)
                grid.dx = header.get('cellsize', 1.0)
                grid.dy = header.get('cellsize', 1.0)
                grid.xmin = header.get('xllcorner', 0.0)
                grid.ymin = header.get('yllcorner', 0.0)
                grid.xmax = grid.xmin + grid.nx * grid.dx
                grid.ymax = grid.ymin + grid.ny * grid.dy
                
                # Read data
                data_list = []
                for line in f:
                    line = line.strip()
                    if line:
                        row = [float(x) for x in line.split()]
                        if len(row) == grid.nx:
                            data_list.append(row)
                
                # Convert to numpy array (note: GRD files are row-major, first row is north)
                if len(data_list) == 0:
                    raise ValueError("No data rows read from file")
                
                if grid.ny == 0 or grid.nx == 0:
                    raise ValueError(f"Invalid grid dimensions: {grid.nx} columns × {grid.ny} rows")
                
                if len(data_list) == grid.ny:
                    # Reverse rows to match standard convention (first row at bottom)
                    data_list.reverse()
                    grid.data = np.array(data_list, dtype=np.float64)
                    
                    # Ensure 2D array
                    if grid.data.ndim == 1:
                        # Reshape if needed
                        grid.data = grid.data.reshape((grid.ny, grid.nx))
                    elif grid.data.shape != (grid.ny, grid.nx):
                        # Fix shape if incorrect
                        print(f"Warning: Data shape {grid.data.shape} doesn't match expected ({grid.ny}, {grid.nx})")
                        if grid.data.size == grid.ny * grid.nx:
                            grid.data = grid.data.reshape((grid.ny, grid.nx))
                        else:
                            raise ValueError(
                                f"Data size mismatch: got {grid.data.size} values, "
                                f"expected {grid.ny * grid.nx}"
                            )
                    
                    grid.update_range()
                else:
                    raise ValueError(
                        f"Data size mismatch: expected {grid.ny} rows, got {len(data_list)}. "
                        f"Grid dimensions: {grid.nx} columns × {grid.ny} rows"
                    )
        
        except Exception as e:
            # If ASCII reading fails and we haven't tried NetCDF yet, try NetCDF
            if HAS_NETCDF and not GRDReader._is_netcdf_file(str(path)):
                try:
                    print(f"ASCII reading failed, trying NetCDF format...")
                    grid = GRDReader._read_netcdf_grd(str(path))
                except Exception as e2:
                    raise IOError(
                        f"Failed to read GRD file {filename} as both ASCII and NetCDF. "
                        f"ASCII error: {str(e)}, NetCDF error: {str(e2)}"
                    )
            else:
                raise IOError(
                    f"Failed to read GRD file {filename}: {str(e)}. "
                    f"If this is a binary NetCDF file, install netCDF4: pip install netCDF4"
                )
        
        return grid
    
    @staticmethod
    def _read_surfer_ascii_grd(filename: str) -> GridData:
        """
        Read Surfer ASCII GRD format
        
        Format:
        DSAA
        nx ny
        xlo xhi
        ylo yhi
        zlo zhi
        data values (whitespace-separated, row by row)
        
        Args:
            filename: Path to Surfer ASCII GRD file
            
        Returns:
            GridData object containing the grid information
        """
        grid = GridData()
        
        with open(filename, 'r') as f:
            # Check header
            header = f.readline().strip()
            if header != 'DSAA':
                raise ValueError(f"Not a Surfer ASCII GRD file. Expected 'DSAA', got '{header}'")
            
            # Read dimensions
            nx, ny = map(int, f.readline().split())
            grid.nx = nx
            grid.ny = ny
            
            # Read coordinate ranges
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())
            
            grid.xmin = xlo
            grid.xmax = xhi
            grid.ymin = ylo
            grid.ymax = yhi
            
            # Calculate grid spacing
            if nx > 1:
                grid.dx = (xhi - xlo) / (nx - 1)
            else:
                grid.dx = 1.0
            
            if ny > 1:
                grid.dy = (yhi - ylo) / (ny - 1)
            else:
                grid.dy = 1.0
            
            # Read data values (all remaining lines)
            data_list = []
            for line in f:
                line = line.strip()
                if line:
                    # Split by whitespace and convert to float
                    values = line.split()
                    data_list.extend([float(v) for v in values])
            
            # Convert to numpy array
            data_array = np.array(data_list, dtype=np.float64)
            
            # Validate size
            expected_size = nx * ny
            if data_array.size != expected_size:
                raise ValueError(
                    f"Grid size mismatch: expected {expected_size} values, got {data_array.size}"
                )
            
            # Reshape to (ny, nx) - Surfer format is row-major, first row is at ylo
            grid.data = data_array.reshape((ny, nx))
            
            # Update z range from actual data (more accurate than header)
            grid.update_range()
            
            # Verify z range from header matches data range (with tolerance)
            if abs(grid.z_min - zlo) > 0.1 or abs(grid.z_max - zhi) > 0.1:
                # Header range might be approximate, use actual data range
                pass  # We already updated from data
        
        return grid
    
    @staticmethod
    def _read_netcdf_grd(filename: str) -> GridData:
        """
        Read NetCDF GRD format (GMT binary format)
        
        Args:
            filename: Path to NetCDF GRD file
            
        Returns:
            GridData object containing the grid information
        """
        if not HAS_NETCDF:
            raise ImportError(
                "netCDF4 library is required to read binary GRD files. "
                "Install it with: pip install netCDF4"
            )
        
        grid = GridData()
        
        with netCDF4.Dataset(filename, 'r') as nc:
            # GMT NetCDF grids typically have these variables/dimensions
            # Try common variable names
            var_names = ['z', 'z_values', 'data', 'elevation', 'topography', 'bathymetry']
            data_var = None
            
            for var_name in var_names:
                if var_name in nc.variables:
                    data_var = nc.variables[var_name]
                    break
            
            if data_var is None:
                # Try to find the first 2D variable
                for var_name in nc.variables:
                    var = nc.variables[var_name]
                    if len(var.dimensions) == 2:
                        data_var = var
                        break
            
            if data_var is None:
                raise ValueError("Could not find 2D data variable in NetCDF file")
            
            # Get dimensions
            dims = data_var.dimensions
            if len(dims) != 2:
                raise ValueError(f"Expected 2D data, found {len(dims)}D")
            
            # Get dimension sizes
            x_dim = nc.dimensions[dims[1]]
            y_dim = nc.dimensions[dims[0]]
            grid.nx = len(x_dim)
            grid.ny = len(y_dim)
            
            # Read data
            grid.data = np.array(data_var[:], dtype=np.float64)
            
            # Handle masked arrays
            if np.ma.is_masked(grid.data):
                grid.data = np.ma.filled(grid.data, np.nan)
            
            # Try to get coordinate information
            # Common attribute names
            if 'x_range' in nc.variables or 'x_range' in nc.ncattrs():
                if 'x_range' in nc.variables:
                    x_range = nc.variables['x_range'][:]
                else:
                    x_range = getattr(nc, 'x_range', None)
                if x_range is not None:
                    grid.xmin = float(x_range[0])
                    grid.xmax = float(x_range[1])
                    grid.dx = (grid.xmax - grid.xmin) / (grid.nx - 1) if grid.nx > 1 else 1.0
            elif dims[1] in nc.variables:
                x_var = nc.variables[dims[1]]
                if len(x_var) == grid.nx:
                    grid.xmin = float(x_var[0])
                    grid.xmax = float(x_var[-1])
                    grid.dx = (grid.xmax - grid.xmin) / (grid.nx - 1) if grid.nx > 1 else float(x_var[1] - x_var[0]) if len(x_var) > 1 else 1.0
            
            if 'y_range' in nc.variables or 'y_range' in nc.ncattrs():
                if 'y_range' in nc.variables:
                    y_range = nc.variables['y_range'][:]
                else:
                    y_range = getattr(nc, 'y_range', None)
                if y_range is not None:
                    grid.ymin = float(y_range[0])
                    grid.ymax = float(y_range[1])
                    grid.dy = (grid.ymax - grid.ymin) / (grid.ny - 1) if grid.ny > 1 else 1.0
            elif dims[0] in nc.variables:
                y_var = nc.variables[dims[0]]
                if len(y_var) == grid.ny:
                    grid.ymin = float(y_var[0])
                    grid.ymax = float(y_var[-1])
                    grid.dy = (grid.ymax - grid.ymin) / (grid.ny - 1) if grid.ny > 1 else float(y_var[1] - y_var[0]) if len(y_var) > 1 else 1.0
            
            # If coordinates not found, use defaults
            if grid.dx == 0.0:
                grid.dx = 1.0
                grid.xmin = 0.0
                grid.xmax = grid.nx * grid.dx
            if grid.dy == 0.0:
                grid.dy = 1.0
                grid.ymin = 0.0
                grid.ymax = grid.ny * grid.dy
            
            # Ensure data is 2D and correct orientation
            if grid.data.ndim == 1:
                grid.data = grid.data.reshape((grid.ny, grid.nx))
            elif grid.data.shape != (grid.ny, grid.nx):
                # May need to transpose
                if grid.data.shape == (grid.nx, grid.ny):
                    grid.data = grid.data.T
                else:
                    raise ValueError(
                        f"Data shape {grid.data.shape} doesn't match dimensions "
                        f"({grid.ny}, {grid.nx})"
                    )
            
            grid.update_range()
        
        return grid
    
    @staticmethod
    def write_grd(filename: str, grid: GridData):
        """
        Write a GRD file
        
        Args:
            filename: Output file path
            grid: GridData object to write
        """
        if grid.data is None:
            raise ValueError("Cannot write empty grid")
        
        path = Path(filename)
        
        # Prepare data (reverse rows for GRD format)
        data_to_write = np.flipud(grid.data.copy())
        
        with open(path, 'w') as f:
            # Write header
            f.write(f"ncols         {grid.nx}\n")
            f.write(f"nrows         {grid.ny}\n")
            f.write(f"xllcorner     {grid.xmin:.8f}\n")
            f.write(f"yllcorner     {grid.ymin:.8f}\n")
            f.write(f"cellsize      {grid.dx:.8f}\n")
            f.write(f"NODATA_value  -9999\n")
            
            # Write data
            for i in range(grid.ny):
                row = data_to_write[i, :]
                f.write(" ".join(f"{val:.6e}" for val in row))
                f.write("\n")

