"""
Script to visualize input data and output results:
- Original Topography (input)
- Original Moho Depth (input)
- Elastic Thickness Te (output)
- Residual Moho (output)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import griddata
import os
import sys

# Add current directory to path to import modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from elastic_thickness_inversion import ElasticThicknessInversion
from constants import (
    RHO_LOAD,
    RHO_MANTLE,
    RHO_INFILL,
    GRAVITY,
)


def load_results(npz_path):
    """Load saved inversion results"""
    data = np.load(npz_path, allow_pickle=True)
    
    # Get basic data
    topography = data["topography"]
    moho_depth = data["moho_depth"]
    topo_anom = data["topography_anomaly"]
    moho_undulation = data["moho_undulation"]
    X = data["X"]
    Y = data["Y"]
    
    # Get Te map (use first available shift)
    te_maps = {}
    for key in data.keys():
        if key.startswith("Te_map_shift_"):
            shift_km = int(key.replace("Te_map_shift_", "").replace("km", ""))
            te_maps[shift_km] = data[key]
    
    if len(te_maps) == 0:
        raise ValueError("No Te maps found in saved data!")
    
    # Use the first available shift
    shift_km = sorted(te_maps.keys())[0]
    te_map = te_maps[shift_km]
    print(f"Using Te map from shift = {shift_km} km")
    
    # Reconstruct center coordinates (for backward compatibility)
    ny_te, nx_te = te_map.shape
    x_m = X[0, :]
    y_m = Y[:, 0]
    te_x_centers = np.linspace(x_m.min(), x_m.max(), nx_te)
    te_y_centers = np.linspace(y_m.min(), y_m.max(), ny_te)
    
    return {
        "topography": topography,
        "moho_depth": moho_depth,
        "topo_anom": topo_anom,
        "moho_undulation": moho_undulation,
        "X": X,
        "Y": Y,
        "te_map": te_map,
        "x_centers": te_x_centers,
        "y_centers": te_y_centers,
        "shift_km": shift_km,
    }


def interpolate_te_to_full_grid(te_map, x_centers, y_centers, X_full, Y_full):
    """Interpolate Te map from window centers to full grid"""
    # Create meshgrid of centers
    X_centers, Y_centers = np.meshgrid(x_centers, y_centers)
    
    # Flatten for interpolation
    points = np.column_stack([X_centers.flatten(), Y_centers.flatten()])
    values = te_map.flatten()
    
    # Remove NaN values
    valid = ~np.isnan(values)
    points_valid = points[valid]
    values_valid = values[valid]
    
    if len(points_valid) == 0:
        raise ValueError("No valid Te values found!")
    
    # Interpolate to full grid
    points_full = np.column_stack([X_full.flatten(), Y_full.flatten()])
    te_interp_flat = griddata(
        points_valid, values_valid, points_full, method="linear", fill_value=np.nan
    )
    
    # Reshape to 2D
    te_interp = te_interp_flat.reshape(X_full.shape)
    
    # Fill NaN with nearest neighbor for edge regions
    if np.any(np.isnan(te_interp)):
        te_interp_nn = griddata(
            points_valid, values_valid, points_full, method="nearest", fill_value=np.nan
        ).reshape(X_full.shape)
        te_interp[np.isnan(te_interp)] = te_interp_nn[np.isnan(te_interp)]
    
    return te_interp


def plot_input_output_comparison(
    X, Y, topography, moho_depth, te_map, te_x_centers, te_y_centers, 
    moho_undulation, topo_anom, output_path, shift_km, dx, dy
):
    """Create 2x2 figure showing inputs and outputs"""
    
    # Calculate extent
    x_m = X[0, :]
    y_m = Y[:, 0]
    extent = [x_m.min(), x_m.max(), y_m.min(), y_m.max()]
    
    # Interpolate Te to full grid
    X_full, Y_full = np.meshgrid(x_m, y_m)
    te_interp = interpolate_te_to_full_grid(te_map, te_x_centers, te_y_centers, X_full, Y_full)
    
    # Calculate predicted and residual Moho
    print("  Calculating predicted Moho using mean Te...")
    inverter = ElasticThicknessInversion(
        dx=dx, dy=dy, rho_load=RHO_LOAD, rho_m=RHO_MANTLE, rho_infill=0, g=GRAVITY
    )
    te_mean = np.nanmean(te_interp)
    print(f"  Mean Te: {te_mean/1000:.1f} km")
    
    # Use topography anomaly for prediction (as done in inversion)
    moho_pred = inverter.predict_moho_flexure(topo_anom, te_mean)
    
    # Residual = observed moho_undulation - predicted moho_undulation
    residual_moho = moho_undulation - moho_pred
    
    print(f"  Residual Moho range: {residual_moho.min():.1f} to {residual_moho.max():.1f} m")
    print(f"  Residual Moho RMS: {np.sqrt(np.mean(residual_moho**2)):.1f} m")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # Formatter
    def format_func(value, tick_number):
        return f"{int(value)}"
    
    # Panel 1: Original Topography
    im1 = axes[0, 0].imshow(
        topography, extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    axes[0, 0].set_title("Topography (m)", fontsize=14, fontweight="bold")
    axes[0, 0].set_xlabel("X (m)", fontsize=12)
    axes[0, 0].set_ylabel("Y (m)", fontsize=12)
    axes[0, 0].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[0, 0].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im1, ax=axes[0, 0], label="Topography (m)")
    
    # Panel 2: Original Moho Depth
    im2 = axes[0, 1].imshow(
        moho_depth, extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    axes[0, 1].set_title("Moho depth (m)", fontsize=14, fontweight="bold")
    axes[0, 1].set_xlabel("X (m)", fontsize=12)
    axes[0, 1].set_ylabel("Y (m)", fontsize=12)
    axes[0, 1].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[0, 1].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im2, ax=axes[0, 1], label="Moho depth (m)")
    
    # Panel 3: Elastic Thickness (Te)
    im3 = axes[1, 0].imshow(
        te_interp / 1000, extent=extent, cmap="jet", origin="lower", aspect="equal",
        vmin=np.nanpercentile(te_interp, 5) / 1000,
        vmax=np.nanpercentile(te_interp, 95) / 1000,
    )
    axes[1, 0].set_title(f"Te (km)", fontsize=14, fontweight="bold")
    axes[1, 0].set_xlabel("X (m)", fontsize=12)
    axes[1, 0].set_ylabel("Y (m)", fontsize=12)
    axes[1, 0].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[1, 0].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im3, ax=axes[1, 0], label="Te (km)")
    
    # Panel 4: Residual Moho
    im4 = axes[1, 1].imshow(
        residual_moho, extent=extent, cmap="RdBu_r", origin="lower", aspect="equal"
    )
    axes[1, 1].set_title("Residual Moho (m)", fontsize=14, fontweight="bold")
    axes[1, 1].set_xlabel("X (m)", fontsize=12)
    axes[1, 1].set_ylabel("Y (m)", fontsize=12)
    axes[1, 1].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[1, 1].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im4, ax=axes[1, 1], label="Residual Moho (m)")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f" Saved: {output_path}")
    plt.show()
    
    return fig


def main():
    """Main function"""
    # Get output folder path
    if len(sys.argv) > 1:
        output_folder = sys.argv[1]
    else:
        output_folder = input("Enter path to output folder: ").strip().strip('"').strip("'")
    
    npz_path = os.path.join(output_folder, "inversion_results.npz")
    
    if not os.path.exists(npz_path):
        print(f" Error: {npz_path} not found!")
        return
    
    print(f"\nLoading data from: {npz_path}")
    data = load_results(npz_path)
    
    # Get Te map coordinates
    te_map = data["te_map"]
    te_x_centers = data["x_centers"]
    te_y_centers = data["y_centers"]
    shift_km = data["shift_km"]
    
    # Get full grid info
    X = data["X"]
    Y = data["Y"]
    x_m = X[0, :]
    y_m = Y[:, 0]
    dx = x_m[1] - x_m[0] if len(x_m) > 1 else 1000
    dy = y_m[1] - y_m[0] if len(y_m) > 1 else 1000
    
    # Calculate residual Moho
    print("\nCalculating predicted and residual Moho...")
    
    # Create plot
    output_path = os.path.join(output_folder, "input_output_comparison.png")
    plot_input_output_comparison(
        X, Y,
        data["topography"],
        data["moho_depth"],
        te_map,
        te_x_centers,
        te_y_centers,
        data["moho_undulation"],
        data["topo_anom"],
        output_path,
        shift_km,
        dx,
        dy
    )
    
    print(f"\n Complete! Results saved to: {output_path}")


if __name__ == "__main__":
    main()
