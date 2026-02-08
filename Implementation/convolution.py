"""
Convolution Operations for Flexural Modeling
Implements 2D convolution for flexural modeling
"""

import numpy as np
from scipy import signal
from typing import Tuple
from grd_reader import GridData


class Convolution:
    """Class for performing convolution operations"""
    
    @staticmethod
    def convolve(load: GridData, response: GridData, method: str = 'fft') -> GridData:
        """
        Perform 2D convolution of load with response function
        
        w(x,y) = ∫∫ R(x-x', y-y') * q(x', y') dx' dy'
        
        Args:
            load: GridData representing the load
            response: GridData representing the unit impulse response
            method: 'fft' for FFT-based (fast) or 'spatial' for spatial convolution
            
        Returns:
            GridData containing the deflection
        """
        if method == 'fft':
            return Convolution.convolve_fft(load, response)
        else:
            return Convolution.convolve_spatial(load, response)
    
    @staticmethod
    def convolve_fft(load: GridData, response: GridData) -> GridData:
        """
        Perform 2D convolution using FFT (fast for large grids)
        
        Args:
            load: GridData representing the load
            response: GridData representing the unit impulse response
            
        Returns:
            GridData containing the deflection
        """
        deflection = GridData()
        deflection.nx = load.nx
        deflection.ny = load.ny
        deflection.dx = load.dx
        deflection.dy = load.dy
        deflection.xmin = load.xmin
        deflection.xmax = load.xmax
        deflection.ymin = load.ymin
        deflection.ymax = load.ymax
        
        # Ensure response is centered properly
        # Shift response so center is at (0,0)
        response_shifted = np.fft.fftshift(response.data)
        
        # Perform FFT convolution
        load_fft = np.fft.fft2(load.data)
        response_fft = np.fft.fft2(response_shifted, s=load.data.shape)
        
        # Multiply and inverse FFT
        result_fft = load_fft * response_fft
        deflection.data = np.real(np.fft.ifft2(result_fft))
        
        # Scale by grid spacing
        deflection.data *= load.dx * load.dy * 1e6  # Convert to meters
        
        deflection.update_range()
        return deflection
    
    @staticmethod
    def convolve_spatial(load: GridData, response: GridData) -> GridData:
        """
        Perform spatial convolution (slower but simpler)
        
        Args:
            load: GridData representing the load
            response: GridData representing the unit impulse response
            
        Returns:
            GridData containing the deflection
        """
        deflection = GridData()
        deflection.nx = load.nx
        deflection.ny = load.ny
        deflection.dx = load.dx
        deflection.dy = load.dy
        deflection.xmin = load.xmin
        deflection.xmax = load.xmax
        deflection.ymin = load.ymin
        deflection.ymax = load.ymax
        
        # Use scipy's convolution for efficiency
        # Flip response for convolution (correlation -> convolution)
        response_flipped = np.flipud(np.fliplr(response.data))
        
        # Perform 2D convolution with 'same' mode to preserve size
        deflection.data = signal.convolve2d(
            load.data, response_flipped, mode='same', boundary='fill', fillvalue=0.0
        )
        
        # Scale by grid spacing
        deflection.data *= load.dx * load.dy * 1e6  # Convert to meters
        
        deflection.update_range()
        return deflection
    
    @staticmethod
    def compute_misfit(observed: GridData, modeled: GridData) -> float:
        """
        Compute misfit between observed and modeled data (RMS error)
        
        Args:
            observed: Observed data grid
            modeled: Modeled data grid
            
        Returns:
            Root mean square error
        """
        if observed.nx != modeled.nx or observed.ny != modeled.ny:
            raise ValueError("Grid dimensions don't match for misfit calculation")
        
        # Compute RMS error
        diff = observed.data - modeled.data
        valid_mask = np.isfinite(diff)
        
        if np.sum(valid_mask) == 0:
            return float('inf')
        
        mse = np.mean(diff[valid_mask]**2)
        return np.sqrt(mse)
    
    @staticmethod
    def compute_correlation(grid1: GridData, grid2: GridData) -> float:
        """
        Compute correlation coefficient between two grids
        
        Args:
            grid1: First grid
            grid2: Second grid
            
        Returns:
            Correlation coefficient
        """
        if grid1.nx != grid2.nx or grid1.ny != grid2.ny:
            raise ValueError("Grid dimensions don't match for correlation calculation")
        
        # Flatten arrays
        data1 = grid1.data.flatten()
        data2 = grid2.data.flatten()
        
        # Remove NaN/inf values
        valid_mask = np.isfinite(data1) & np.isfinite(data2)
        
        if np.sum(valid_mask) < 2:
            return 0.0
        
        data1_valid = data1[valid_mask]
        data2_valid = data2[valid_mask]
        
        # Compute correlation
        correlation_matrix = np.corrcoef(data1_valid, data2_valid)
        return float(correlation_matrix[0, 1])

