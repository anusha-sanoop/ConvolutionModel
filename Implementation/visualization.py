"""
Visualization functions for inverse modeling results
Creates and saves plots of Te maps, deflection, misfit, etc.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from grd_reader import GridData


def save_te_map(te_map: GridData, output_path: str, title: str = "Estimated Elastic Thickness (Te)"):
    """
    Save visualization of Te map
    
    Args:
        te_map: GridData containing Te values
        output_path: Path to save the image
        title: Plot title
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create coordinate arrays
    x = np.linspace(te_map.xmin, te_map.xmax, te_map.nx)
    y = np.linspace(te_map.ymin, te_map.ymax, te_map.ny)
    X, Y = np.meshgrid(x, y)
    
    # Create plot
    im = ax.contourf(X / 1000.0, Y / 1000.0, te_map.data, levels=20, cmap='viridis')
    ax.contour(X / 1000.0, Y / 1000.0, te_map.data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Elastic Thickness (km)', rotation=270, labelpad=20)
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title(title)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_deflection_map(deflection: GridData, output_path: str, title: str = "Modeled Moho Deflection"):
    """
    Save visualization of modeled deflection
    
    Args:
        deflection: GridData containing deflection values
        output_path: Path to save the image
        title: Plot title
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create coordinate arrays
    x = np.linspace(deflection.xmin, deflection.xmax, deflection.nx)
    y = np.linspace(deflection.ymin, deflection.ymax, deflection.ny)
    X, Y = np.meshgrid(x, y)
    
    # Create plot
    im = ax.contourf(X / 1000.0, Y / 1000.0, deflection.data, levels=20, cmap='seismic')
    ax.contour(X / 1000.0, Y / 1000.0, deflection.data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Deflection (m)', rotation=270, labelpad=20)
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title(title)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_misfit_map(misfit: GridData, output_path: str, title: str = "Misfit Map"):
    """
    Save visualization of misfit map
    
    Args:
        misfit: GridData containing misfit values
        output_path: Path to save the image
        title: Plot title
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create coordinate arrays
    x = np.linspace(misfit.xmin, misfit.xmax, misfit.nx)
    y = np.linspace(misfit.ymin, misfit.ymax, misfit.ny)
    X, Y = np.meshgrid(x, y)
    
    # Use log scale for misfit if values vary widely
    misfit_data = misfit.data
    if np.nanmax(misfit_data) / np.nanmin(misfit_data[np.isfinite(misfit_data) & (misfit_data > 0)]) > 100:
        # Use log scale
        misfit_plot = np.log10(np.clip(misfit_data, 1e-10, None))
        cbar_label = 'Log10(Misfit)'
    else:
        misfit_plot = misfit_data
        cbar_label = 'Misfit (m²)'
    
    # Create plot
    im = ax.contourf(X / 1000.0, Y / 1000.0, misfit_plot, levels=20, cmap='hot')
    ax.contour(X / 1000.0, Y / 1000.0, misfit_plot, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, rotation=270, labelpad=20)
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title(title)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_comparison_plot(observed: GridData, modeled: GridData, output_path: str, 
                        title: str = "Observed vs Modeled Comparison"):
    """
    Save side-by-side comparison of observed and modeled data
    
    Args:
        observed: GridData containing observed values
        modeled: GridData containing modeled values
        output_path: Path to save the image
        title: Plot title
    """
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    
    # Create coordinate arrays
    x = np.linspace(observed.xmin, observed.xmax, observed.nx)
    y = np.linspace(observed.ymin, observed.ymax, observed.ny)
    X, Y = np.meshgrid(x, y)
    
    # Observed
    im1 = ax1.contourf(X / 1000.0, Y / 1000.0, observed.data, levels=20, cmap='terrain')
    ax1.contour(X / 1000.0, Y / 1000.0, observed.data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax1.set_title('Observed')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1, label='Depth (m)')
    
    # Modeled
    im2 = ax2.contourf(X / 1000.0, Y / 1000.0, modeled.data, levels=20, cmap='terrain')
    ax2.contour(X / 1000.0, Y / 1000.0, modeled.data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax2.set_title('Modeled')
    ax2.set_xlabel('X (km)')
    ax2.set_ylabel('Y (km)')
    ax2.set_aspect('equal')
    plt.colorbar(im2, ax=ax2, label='Depth (m)')
    
    # Difference
    diff = observed.data - modeled.data
    im3 = ax3.contourf(X / 1000.0, Y / 1000.0, diff, levels=20, cmap='RdBu_r')
    ax3.contour(X / 1000.0, Y / 1000.0, diff, levels=10, colors='black', alpha=0.3, linewidths=0.5)
    ax3.set_title('Difference (Observed - Modeled)')
    ax3.set_xlabel('X (km)')
    ax3.set_ylabel('Y (km)')
    ax3.set_aspect('equal')
    plt.colorbar(im3, ax=ax3, label='Difference (m)')
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_te_vs_misfit_plot(te_values: list, misfits: list, output_path: str, 
                           title: str = "Te vs Misfit"):
    """
    Save plot of Te values vs misfit
    
    Args:
        te_values: List of Te values tested
        misfits: List of corresponding misfit values
        output_path: Path to save the image
        title: Plot title
    """
    if not te_values or not misfits:
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(te_values, misfits, 'b-o', linewidth=2, markersize=6)
    ax.set_xlabel('Elastic Thickness (Te, km)')
    ax.set_ylabel('Misfit (RMS)')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Mark minimum
    min_idx = np.argmin(misfits)
    ax.plot(te_values[min_idx], misfits[min_idx], 'ro', markersize=12, label=f'Minimum at Te={te_values[min_idx]:.1f} km')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_shift_comparison_grid(te_maps: list, rms_maps: list, shift_values: list, 
                               output_path: str, title_prefix: str = ""):
    """
    Save a 2x3 grid comparing Te maps and RMS maps for different shift values
    
    Args:
        te_maps: List of GridData objects containing Te maps for each shift
        rms_maps: List of GridData objects containing RMS maps for each shift
        shift_values: List of shift values in km
        output_path: Path to save the image
        title_prefix: Optional prefix for the title
    """
    if len(te_maps) != len(rms_maps) or len(te_maps) != len(shift_values):
        raise ValueError("te_maps, rms_maps, and shift_values must have the same length")
    
    if len(shift_values) == 0:
        return
    
    # Create 2 rows x N columns grid
    n_shifts = len(shift_values)
    fig, axes = plt.subplots(2, n_shifts, figsize=(6*n_shifts, 10))
    
    if n_shifts == 1:
        axes = axes.reshape(2, 1)
    
    # Get coordinate ranges from first map (assuming all have same extent)
    if len(te_maps) > 0:
        x = np.linspace(te_maps[0].xmin, te_maps[0].xmax, te_maps[0].nx)
        y = np.linspace(te_maps[0].ymin, te_maps[0].ymax, te_maps[0].ny)
        X, Y = np.meshgrid(x, y)
    else:
        return
    
    # Top row: Te Maps
    for i, (te_map, shift) in enumerate(zip(te_maps, shift_values)):
        ax = axes[0, i]
        
        # Plot Te map
        im = ax.contourf(X / 1000.0, Y / 1000.0, te_map.data, levels=20, cmap='viridis')
        ax.contour(X / 1000.0, Y / 1000.0, te_map.data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Te (km)', rotation=270, labelpad=15)
        
        ax.set_xlabel('X (km)')
        if i == 0:
            ax.set_ylabel('Y (km)')
        ax.set_title(f'Te Map - Shift {shift:.0f} km')
        ax.set_aspect('equal')
    
    # Bottom row: RMS Maps
    for i, (rms_map, shift) in enumerate(zip(rms_maps, shift_values)):
        ax = axes[1, i]
        
        # RMS data should already be in meters (not squared)
        rms_data = rms_map.data
        
        # Plot RMS map with hot colormap
        im = ax.contourf(X / 1000.0, Y / 1000.0, rms_data, levels=20, cmap='hot')
        ax.contour(X / 1000.0, Y / 1000.0, rms_data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('RMS (m of Moho)', rotation=270, labelpad=15)
        
        ax.set_xlabel('X (km)')
        if i == 0:
            ax.set_ylabel('Y (km)')
        ax.set_title(f'RMS Map - Shift {shift:.0f} km')
        ax.set_aspect('equal')
    
    # Add overall title
    if title_prefix:
        fig.suptitle(title_prefix, fontsize=14, y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96] if title_prefix else None)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

