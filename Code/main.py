"""
Main script for Elastic Thickness Inversion
Updated version with single large Te map
"""

import numpy as np
import os
import sys
from datetime import datetime
import matplotlib.pyplot as plt

from data_loader import (
    read_surfer_grd,
    write_surfer_grd,
    check_grid_compatibility,
    prepare_data_for_inversion,
<<<<<<< HEAD
    write_surfer_grd,
=======
>>>>>>> e164d71a87412ccae8297068cbd15b9f29754aad
)
from elastic_thickness_inversion import ElasticThicknessInversion
from moving_window_analysis import MovingWindowAnalysis
from gravity_moho_inversion import GravityMohoInversion
from constants import (
    RHO_LOAD,
    RHO_MANTLE,
    RHO_INFILL,
    GRAVITY,
    YOUNGS_MODULUS,
    POISSONS_RATIO,
    DENSITY_CONTRAST,
    DEFAULT_REFERENCE_MOHO_DEPTH,
    GRAVITATIONAL_CONSTANT,
)


class TeeOutput:
    """Class to duplicate output to both console and file"""
    def __init__(self, file_path):
        self.terminal = sys.stdout
        self.log_file = open(file_path, 'w', encoding='utf-8')
    
    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()  # Ensure immediate write
    
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()
    
    def close(self):
        if self.log_file:
            self.log_file.close()
            sys.stdout = self.terminal


def create_output_folder():
    """Create timestamped output folder"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"Output_{timestamp}"
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), folder_name)
    os.makedirs(output_path, exist_ok=True)
    print(f"\nCreated output folder: {output_path}")
    return output_path


def create_filename_prefix(rho_load, ref_moho_km, window_size_km):
    """
    Create filename prefix: crustaldensity_referencedepth(in km)_windowsize
    
    Parameters:
    -----------
    rho_load : float
        Crustal density in kg/m³
    ref_moho_km : float
        Reference Moho depth in km
    window_size_km : float
        Window size in km
    
    Returns:
    --------
    prefix : str
        Filename prefix (e.g., "2.9_30_2160_")
    """
    # Convert density from kg/m³ to g/cm³ (divide by 1000)
    density_gcm3 = rho_load / 1000.0
    
    # Format: density_refmoho_km_window_km (no trailing underscore)
    prefix = f"{density_gcm3:.1f}_{ref_moho_km:.0f}_{window_size_km:.0f}_"
    return prefix


def plot_final_results(
    X, Y, topography, moho_depth,
    moho_undulation, topo_anom, output_folder, dx, dy, use_gravity=False, filename_prefix=""
):
    """Plot final results: Topography and Moho Depth (1x2 grid, 2 panels)"""
    from matplotlib.ticker import FuncFormatter
    
    # Set toolbar to None to remove buttons
    plt.rcParams['toolbar'] = 'None'

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))  # Changed to 1x2 layout for 2 panels
    # Flatten axes for easier indexing (though not needed for 1x2, keeping for consistency)
    axes = axes.flatten()

    # Keep axes in METERS to match target image style
    x_m = X[0, :]
    y_m = Y[:, 0]
    extent = [x_m.min(), x_m.max(), y_m.min(), y_m.max()]

    # Formatter to disable scientific notation
    def format_func(value, tick_number):
        return f"{int(value)}"

    # Panel 1: Topography (top-left)
    im1 = axes[0].imshow(
        topography, extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    axes[0].set_title("Topography (m)", fontsize=14, fontweight="bold")
    axes[0].set_xlabel("X (m)", fontsize=12)
    axes[0].set_ylabel("Y (m)", fontsize=12)
    axes[0].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[0].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im1, ax=axes[0], label="Topography (m)")

    # Panel 2: Moho depth (top-right)
    moho_title = "Moho depth (m) [from gravity]" if use_gravity else "Moho depth (m)"
    im2 = axes[1].imshow(
        moho_depth, extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    axes[1].set_title(moho_title, fontsize=14, fontweight="bold")
    axes[1].set_xlabel("X (m)", fontsize=12)
    axes[1].set_ylabel("Y (m)", fontsize=12)
    axes[1].xaxis.set_major_formatter(FuncFormatter(format_func))
    axes[1].yaxis.set_major_formatter(FuncFormatter(format_func))
    plt.colorbar(im2, ax=axes[1], label="Moho depth (m)")

    # ========== RESIDUAL MOHO PLOT (COMMENTED OUT - TEMPORARY) ==========
    # TODO: Uncomment this section if residual Moho plot is needed in the future
    # 
    # # Panel 3: Residual Moho (bottom-left)
    # # Calculate predicted Moho using flexure model and compare with observed/predicted Moho
    # inverter = ElasticThicknessInversion(
    #     dx=dx, dy=dy, rho_load=2900, rho_m=3500, rho_infill=0, g=3.72
    # )
    # 
    # if use_gravity:
    #     # If Moho from gravity: compare gravity-predicted vs flexure-predicted
    #     # Use a reasonable default Te for flexure prediction
    #     te_mean = 30000  # 30 km default
    #     moho_pred_flexure = inverter.predict_moho_flexure(topo_anom, te_mean)
    #     # Residual = gravity Moho - flexure Moho
    #     residual_moho = moho_undulation - moho_pred_flexure
    # else:
    #     # If Moho from observation: compare observed vs flexure-predicted
    #     # Use mean Te from moving window if available, otherwise default
    #     te_mean = 30000  # Default 30 km
    #     moho_pred = inverter.predict_moho_flexure(topo_anom, te_mean)
    #     residual_moho = moho_undulation - moho_pred
    # 
    # # Mask edge regions to reduce edge effect artifacts
    # # Calculate edge mask based on taper region (typically 10% of domain)
    # ny, nx = residual_moho.shape
    # edge_fraction = 0.1  # Same as taper_alpha
    # edge_x = int(nx * edge_fraction)
    # edge_y = int(ny * edge_fraction)
    # 
    # # Create mask for edge regions
    # edge_mask = np.ones_like(residual_moho, dtype=bool)
    # edge_mask[:edge_y, :] = False  # Top edge
    # edge_mask[-edge_y:, :] = False  # Bottom edge
    # edge_mask[:, :edge_x] = False  # Left edge
    # edge_mask[:, -edge_x:] = False  # Right edge
    # 
    # # Apply mask to residual (set edge regions to NaN for visualization)
    # residual_moho_masked = residual_moho.copy()
    # residual_moho_masked[~edge_mask] = np.nan
    # 
    # im3 = axes[2].imshow(
    #     residual_moho_masked, extent=extent, cmap="RdBu_r", origin="lower", aspect="equal"
    # )
    # axes[2].set_title("Residual Moho (m)", fontsize=14, fontweight="bold")
    # axes[2].set_xlabel("X (m)", fontsize=12)
    # axes[2].set_ylabel("Y (m)", fontsize=12)
    # axes[2].xaxis.set_major_formatter(FuncFormatter(format_func))
    # axes[2].yaxis.set_major_formatter(FuncFormatter(format_func))
    # plt.colorbar(im3, ax=axes[2], label="Residual Moho (m)")
    # 
    # # Hide the 4th subplot (bottom-right) - only needed for 2x2 layout
    # # axes[3].axis('off')
    #========================================================

    plt.tight_layout()
    
    # Save figure with prefix
    filename = f"{filename_prefix}input_data.png" if filename_prefix else "input_data.png"
    plt.savefig(
        os.path.join(output_folder, filename), dpi=300, bbox_inches="tight"
    )
    
    # Hide toolbar after figure is created
    try:
        if hasattr(fig.canvas, 'toolbar') and fig.canvas.toolbar is not None:
            fig.canvas.toolbar.pack_forget()
        elif hasattr(fig.canvas.manager, 'toolbar') and fig.canvas.manager.toolbar is not None:
            fig.canvas.manager.toolbar.pack_forget()
    except (AttributeError, KeyError, TypeError):
        pass
    
    plt.show(block=False)

    return fig


def plot_input_data_3d(
    X, Y, topography, moho_depth, topo_anom, moho_undulation, output_folder_3d, filename_prefix=""
):
    """Create combined 3D visualization with topography and Moho depth surfaces."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  - needed for 3D plots

    # Coordinates in meters
    x_m = X[0, :]
    y_m = Y[:, 0]
    Xg, Yg = np.meshgrid(x_m, y_m)

    # Create the combined 3D plot
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Plot Topography surface (using topography anomaly for better visualization)
    # Use gist_earth colormap for topography
    surf_topo = ax.plot_surface(
        Xg,
        Yg,
        topo_anom,
        cmap="gist_earth",
        edgecolor="none",
        alpha=0.8,
        linewidth=0,
        antialiased=True,
    )

    # Plot Moho Depth surface (using Moho undulation)
    # Use a colormap that goes from light pink/orange to purple/grey
    # We'll use 'coolwarm' or create a custom one, but 'RdBu_r' works well
    surf_moho = ax.plot_surface(
        Xg,
        Yg,
        moho_undulation,
        cmap="RdBu_r",
        edgecolor="none",
        alpha=0.6,
        linewidth=0,
        antialiased=True,
    )

    # Add reference planes
    # Dark blue plane at Z=0 (sea level)
    z_zero = np.zeros_like(Xg)
    ax.plot_surface(
        Xg,
        Yg,
        z_zero,
        color="darkblue",
        alpha=0.3,
        linewidth=0,
    )

    # Semi-transparent yellow plane at average Moho depth
    moho_mean = np.mean(moho_undulation)
    z_moho_ref = np.full_like(Xg, moho_mean)
    ax.plot_surface(
        Xg,
        Yg,
        z_moho_ref,
        color="yellow",
        alpha=0.2,
        linewidth=0,
    )

    # Set title
    ax.set_title(
        "3D Visualization: Topography and Moho Depth (Airy Isostasy Model)",
        fontsize=15,
        fontweight="bold",
        pad=20,
    )

    # Set axis labels
    ax.set_xlabel("X (m)", fontsize=12)
    ax.set_ylabel("Y (m)", fontsize=12)
    ax.set_zlabel("Elevation/Depth (m)", fontsize=12)

    # Set Z-axis limits to show both surfaces clearly
    z_min = min(moho_undulation.min(), -80000)
    z_max = max(topo_anom.max(), 20000)
    ax.set_zlim(z_min, z_max)

    # Add two separate colorbars
    # Left colorbar for Topography
    cbar_topo = fig.colorbar(
        surf_topo,
        ax=ax,
        shrink=0.5,
        aspect=20,
        pad=0.05,
        location="left",
    )
    cbar_topo.set_label("Topography (m)", fontsize=11, rotation=90, labelpad=15)

    # Right colorbar for Moho Depth
    cbar_moho = fig.colorbar(
        surf_moho,
        ax=ax,
        shrink=0.5,
        aspect=20,
        pad=0.05,
        location="right",
    )
    cbar_moho.set_label("Moho Depth (m)", fontsize=11, rotation=90, labelpad=15)

    plt.tight_layout()

    # Save figure with prefix
    filename = f"{filename_prefix}topography_moho_combined_3d.png" if filename_prefix else "topography_moho_combined_3d.png"
    save_path = os.path.join(output_folder_3d, filename)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f" Saved: {save_path}")
    plt.close(fig)


def main():
    """Main execution function"""

    print("\n" + "=" * 70)
    print(" ELASTIC THICKNESS INVERSION - BRAITENBERG METHOD")
    print("=" * 70)

    # DATA SELECTION
    print("\n" + "=" * 70)
    print(" INPUT FILE PATHS")
    print("=" * 70)
    
    # Get topography file path from user
    topography_file = input("\nEnter path to topography file (.grd): ").strip()
    # Remove quotes if user pasted path with quotes
    topography_file = topography_file.strip('"').strip("'")
    
    # Check if gravity data is available
    use_gravity = input("\nDo you have Bouguer gravity data? (y/n): ").strip().lower()
    
    if use_gravity == "y":
        # Get gravity file path
        gravity_file = input("Enter path to Bouguer gravity file (.grd): ").strip()
        gravity_file = gravity_file.strip('"').strip("'")
        moho_file = None  # Will predict Moho from gravity
        print("→ Will predict Moho from gravity data")
    else:
        # Get Moho depth file path from user
        moho_file = input("Enter path to Moho depth file (.grd): ").strip()
        moho_file = moho_file.strip('"').strip("'")
        gravity_file = None
        print("→ Will use provided Moho depth data")
    
    # No need for true Te value input
    true_Te = None

    # Get reference Moho depth from user (in km only)
    ref_moho_input = input(
        "Enter reference Moho depth (km) [press Enter to use mean]: "
    ).strip()
    if ref_moho_input:
        try:
            ref_moho_km = float(ref_moho_input)
            ref_moho = ref_moho_km * 1000.0  # Convert km to meters for internal calculations
            print(f"→ Using reference Moho depth: {ref_moho_km:.1f} km ({ref_moho:.1f} m)")
        except:
            ref_moho = None
            ref_moho_km = None
            print("→ Invalid input, will use mean Moho depth as reference")
    else:
        ref_moho = None
        ref_moho_km = None
        print("→ Will use mean Moho depth as reference")

    # Check file existence
    print(f"\nChecking files...")
    print(f"  Topography file: {topography_file}")
    
    if not os.path.exists(topography_file):
        print(f"\n ERROR: Topography file not found!")
        print(f"   {topography_file}")
        return
    
    if use_gravity == "y":
        print(f"  Gravity file: {gravity_file}")
        if not os.path.exists(gravity_file):
            print(f"\n ERROR: Gravity file not found!")
            print(f"   {gravity_file}")
            return
        print(" Both files found!")
    else:
        print(f"  Moho file: {moho_file}")
        if not os.path.exists(moho_file):
            print(f"\n ERROR: Moho file not found!")
            print(f"   {moho_file}")
            return
        print(" Both files found!")

    # Create output folder
    output_folder = create_output_folder()
    output_folder_3d = os.path.join(output_folder, "3D")
    os.makedirs(output_folder_3d, exist_ok=True)
    
    # Set up terminal output logging to file
    log_file_path = os.path.join(output_folder, "terminal_output.txt")
    tee = TeeOutput(log_file_path)
    sys.stdout = tee
    print(f"\n{'='*70}")
    print(f" TERMINAL OUTPUT LOGGING ENABLED")
    print(f"{'='*70}")
    print(f"All output is being saved to: {log_file_path}")
    print(f"{'='*70}\n")

    # LOAD DATA
    print("\n" + "=" * 70)
    print(" LOADING DATA")
    print("=" * 70)

    X_topo, Y_topo, topography, dx, dy, *_ = read_surfer_grd(topography_file)
    
    if use_gravity == "y":
        # Load gravity data and predict Moho
        X_grav, Y_grav, bouguer_gravity, *_ = read_surfer_grd(gravity_file)
        
        # Check grid compatibility
        compatible, message = check_grid_compatibility(X_topo, Y_topo, X_grav, Y_grav)
        if not compatible:
            print(f"\n⚠️  WARNING: {message}")
            print("Proceeding with caution...")
        else:
            print(f"\n {message}")
        
        # Predict Moho from gravity
        print("\n" + "=" * 70)
        print(" PREDICTING MOHO FROM GRAVITY")
        print("=" * 70)
        
        # Use user-provided reference Moho, or fallback default if not provided
        if ref_moho is None:
            ref_moho_for_gravity = DEFAULT_REFERENCE_MOHO_DEPTH
            print(f"  Using default reference Moho for gravity inversion: {ref_moho_for_gravity/1000:.1f} km")
        else:
            ref_moho_for_gravity = ref_moho
            print(f"  Using user-provided reference Moho: {ref_moho_for_gravity/1000:.1f} km")
        
        gravity_inverter = GravityMohoInversion(
            dx=dx, dy=dy, density_contrast=DENSITY_CONTRAST, reference_moho=ref_moho_for_gravity
        )
        moho_depth = gravity_inverter.predict_moho_from_gravity(bouguer_gravity)
        X_moho, Y_moho = X_grav, Y_grav
        
        print(f" Moho predicted from gravity data")
        print(f"  Moho depth range: {moho_depth.min():.1f} to {moho_depth.max():.1f} m")
    else:
        # Load observed Moho data
        X_moho, Y_moho, moho_depth, *_ = read_surfer_grd(moho_file)

        # Check grid compatibility
        compatible, message = check_grid_compatibility(X_topo, Y_topo, X_moho, Y_moho)
        if not compatible:
            print(f"\n⚠️  WARNING: {message}")
            print("Proceeding with caution...")
        else:
            print(f"\n {message}")

    # PREPARE DATA
    print("\n" + "=" * 70)
    print(" PREPARING DATA FOR INVERSION")
    print("=" * 70)

    # Apply taper and remove mean (or user-provided reference)
    if ref_moho is not None:
        # Use user-provided reference Moho
        topo_anom, moho_undulation, stats = prepare_data_for_inversion(
            topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1, ref_moho=ref_moho
        )
        if ref_moho_km is not None:
            print(f"  Using user-provided reference Moho: {ref_moho_km:.1f} km ({ref_moho:.1f} m)")
        else:
            print(f"  Using user-provided reference Moho: {ref_moho:.1f} m")
    else:
        # Use mean as reference (default)
        topo_anom, moho_undulation, stats = prepare_data_for_inversion(
            topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1
        )
        print(f"  Using mean Moho depth as reference: {stats['moho_mean']:.1f} m")
    
    # Create initial filename prefix (will be updated with actual window_size later)
    if ref_moho_km is not None:
        # ref_moho_km already set from user input
        pass
    elif ref_moho is not None:
        # Fallback: convert from meters if ref_moho_km wasn't set
        ref_moho_km = ref_moho / 1000.0
    else:
        ref_moho_km = stats['moho_mean'] / 1000.0
    # Use default window size for now (will be updated when window_size is known)
    window_size_km_default = 1000.0  # Default 1000 km
    filename_prefix = create_filename_prefix(RHO_LOAD, ref_moho_km, window_size_km_default)
    print(f"\n→ Initial filename prefix: {filename_prefix}")

    print("\nData Statistics:")
    print(f"  Topography mean: {stats['topo_mean']:.1f} m")
    print(f"  Topography std (original): {stats['topo_std']:.1f} m")
    print(f"  Topography std (anomaly): {stats['topo_anom_std']:.1f} m")
    print(f"  Moho mean depth: {stats['moho_mean']:.1f} m")
    print(f"  Moho std (original): {stats['moho_std']:.1f} m")
    print(f"  Moho undulation std: {stats['moho_undulation_std']:.1f} m")
    print(f"  Taper applied: {stats['tapered']}")

    # PHYSICAL PARAMETERS
    print("\n" + "=" * 70)
    print(" PHYSICAL PARAMETERS (MARS)")
    print("=" * 70)
    print(f"  Crustal density (ρ_crust): {RHO_LOAD} kg/m³")
    print(f"  Mantle density (ρ_mantle): {RHO_MANTLE} kg/m³")
    print(f"  Density contrast (Δρ): {DENSITY_CONTRAST} kg/m³")
    print(f"  Gravity (g): {GRAVITY} m/s²")
    print(f"  Young's modulus (E): {YOUNGS_MODULUS:.1e} Pa")
    print(f"  Poisson's ratio (ν): {POISSONS_RATIO}")

    # DIAGNOSTIC CHECK
    print("\n" + "=" * 70)
    print(" DIAGNOSTIC CHECK")
    print("=" * 70)

    inverter = ElasticThicknessInversion(
        dx=dx,
        dy=dy,
        rho_load=RHO_LOAD,
        rho_m=RHO_MANTLE,
        rho_infill=RHO_INFILL,
        g=GRAVITY,
    )

    diagnostics = inverter.diagnostic_check(topo_anom, moho_undulation)

    # Manual correlation check
    topo_flat = topo_anom.flatten()
    moho_flat = moho_undulation.flatten()
    correlation = np.corrcoef(topo_flat, moho_flat)[0, 1]

    print("\n=== SIGN CONVENTION CHECK ===")
    print(f"Topo-Moho correlation: {correlation:.3f}")

    # Find location of max topography
    max_topo_idx = np.unravel_index(np.abs(topo_anom).argmax(), topo_anom.shape)
    max_topo_val = topo_anom[max_topo_idx]
    moho_at_max_topo = moho_undulation[max_topo_idx]

    print(f"At max |topo| ({max_topo_val:.1f} m):")
    print(f"  Moho undulation: {moho_at_max_topo:.1f} m")

    if max_topo_val > 0 and moho_at_max_topo > 0:
        print("  → Positive topo, Positive Moho deflection")
        print("  → This means: Mountain → Moho goes UP (shallower)")
        print("  ⚠️  This is BACKWARDS! Moho should go DOWN under mountains")
        flip_sign = input("\nFlip Moho sign? (y/n): ").lower()
        if flip_sign == "y":
            print("→ Flipping Moho undulation sign...")
            moho_undulation = -moho_undulation
    elif max_topo_val > 0 and moho_at_max_topo < 0:
        print("  → Positive topo, Negative Moho deflection")
        print("  → This means: Mountain → Moho goes DOWN (deeper)")
        print("   Sign convention is CORRECT")

    print("============================\n")

    # Note: Final results plot will be created after moving window analysis

    # 3D visualization of input data
    plot_input_data_3d(
        X_topo,
        Y_topo,
        topography,
        moho_depth,
        topo_anom,
        moho_undulation,
        output_folder_3d,
        filename_prefix,
    )

    # GLOBAL INVERSION
    # perform_global = input("\nPerform global (full-region) inversion? (y/n): ").lower()
    perform_global = "n"
    if perform_global == "y":
        print("\n" + "=" * 70)
        print(" GLOBAL INVERSION")
        print("=" * 70)

        Te_range = (1000, 100000)
        print(f"Te search range: {Te_range[0] / 1000:.0f}-{Te_range[1] / 1000:.0f} km")

        global_result = inverter.invert_elastic_thickness(
            topo_anom, moho_undulation, Te_range=Te_range, method="bounded"
        )

        print(f"\n Global Te estimate: {global_result['Te_best'] / 1000:.2f} km")
        print(f"  RMS misfit: {global_result['rms_best']:.2f} m")

        if true_Te is not None:
            error = abs(global_result["Te_best"] - true_Te) / true_Te * 100
            print(f"  True Te: {true_Te / 1000:.1f} km")
            print(f"  Error: {error:.1f}%")

        # Sensitivity analysis
        perform_sensitivity = input(
            "\nPerform sensitivity analysis (Te vs RMS curve)? (y/n): "
        ).lower()

        if perform_sensitivity == "y":
            print("\n→ Running sensitivity analysis...")
            Te_values, rms_values = inverter.sensitivity_analysis(
                topo_anom, moho_undulation, Te_range=Te_range, n_points=100
            )

            fig_sens = plt.figure(figsize=(10, 6))
            plt.plot(
                Te_values / 1000,
                rms_values / 1000,
                "b-",
                linewidth=2,
                label="RMS Misfit",
            )
            plt.axvline(
                global_result["Te_best"] / 1000,
                color="r",
                linestyle="--",
                linewidth=2,
                label=f"Best Te = {global_result['Te_best'] / 1000:.1f} km",
            )

            if true_Te is not None:
                plt.axvline(
                    true_Te / 1000,
                    color="g",
                    linestyle="--",
                    linewidth=2,
                    label=f"True Te = {true_Te / 1000:.1f} km",
                )

            plt.xlabel("Elastic Thickness (km)", fontsize=12)
            plt.ylabel("RMS Misfit (km)", fontsize=12)
            plt.title(
                "Sensitivity Analysis: Te vs RMS Misfit", fontsize=14, fontweight="bold"
            )
            plt.grid(True, alpha=0.3)
            plt.legend(fontsize=10)
            plt.tight_layout()

            # Save with prefix if available
            sens_filename = f"{filename_prefix}sensitivity_analysis.png" if filename_prefix else "sensitivity_analysis.png"
            fig_sens.savefig(
                os.path.join(output_folder, sens_filename), dpi=300
            )
            plt.show(block=False)

            min_rms = global_result["rms_best"]
            threshold = min_rms * 1.1
            confidence_mask = rms_values <= threshold
            if np.any(confidence_mask):
                Te_conf = Te_values[confidence_mask]
                print(f"\n  Confidence interval (RMS ≤ {threshold:.1f} m):")
                print(
                    f"    Te range: {Te_conf.min() / 1000:.1f} - {Te_conf.max() / 1000:.1f} km"
                )
                print(f"    Spread: ±{(Te_conf.max() - Te_conf.min()) / 2000:.1f} km")

            rms_range = rms_values.max() - rms_values.min()
            if rms_range / min_rms > 2:
                print("\n   Solution appears well-constrained (clear minimum in RMS)")
            else:
                print("\n  ⚠️  Solution may be poorly constrained (shallow RMS minimum)")

    # MOVING WINDOW ANALYSIS
    # Moving window analysis is needed to estimate Te, regardless of Moho source
    # perform_mw = input("\nPerform moving window analysis? (y/n): ").lower()
    perform_mw = "y"
    
    if perform_mw == "y":
        if use_gravity == "y":
            print("\n→ Performing moving window analysis using gravity-predicted Moho")
            print("  (Will estimate Te by comparing gravity-predicted vs flexure-predicted Moho)")
        else:
            print("\n→ Performing moving window analysis using observed Moho")
            print("  (Will estimate Te by comparing observed vs flexure-predicted Moho)")
        print("\n" + "=" * 70)
        print(" MOVING WINDOW ANALYSIS")
        print("=" * 70)

        try:
            window_size = (
                float(input("Window size (km) [default 1000]: ") or "1000") * 1000
            )
            shift_min = float(input("Minimum shift (km) [default 20]: ") or "20") * 1000
            shift_max = float(input("Maximum shift (km) [default 80]: ") or "80") * 1000
            shift_step = float(input("Shift step (km) [default 20]: ") or "20") * 1000
            Te_min = float(input("Minimum Te (km) [default 5]: ") or "5") * 1000
            Te_max = float(input("Maximum Te (km) [default 80]: ") or "80") * 1000
        except:
            print("Using default parameters...")
            window_size = 1_000_000
            shift_min, shift_max, shift_step = 20_000, 80_000, 20_000
            Te_min, Te_max = 5_000, 80_000

        print("\n→ Analysis Parameters:")
        print(f"  Window size: {window_size / 1000:.0f} km")
        print(
            f"  Shift range: {shift_min / 1000:.0f}-{shift_max / 1000:.0f} km (step {shift_step / 1000:.0f} km)"
        )
        print(f"  Te search range: {Te_min / 1000:.0f}-{Te_max / 1000:.0f} km")
        
        # Update filename prefix with actual window_size
        if ref_moho_km is not None:
            # ref_moho_km already set from user input
            pass
        elif ref_moho is not None:
            # Fallback: convert from meters if ref_moho_km wasn't set
            ref_moho_km = ref_moho / 1000.0
        else:
            ref_moho_km = stats['moho_mean'] / 1000.0
        
        window_size_km = window_size / 1000.0
        filename_prefix = create_filename_prefix(RHO_LOAD, ref_moho_km, window_size_km)
        print(f"\n→ Updated filename prefix with window size: {filename_prefix}")

        mw_analyzer = MovingWindowAnalysis(
            dx=dx, dy=dy, rho_load=RHO_LOAD, rho_m=RHO_MANTLE, rho_infill=RHO_INFILL, g=GRAVITY
        )

        mw_results_dict = mw_analyzer.analyze_multiple_shifts(
            topo_anom,
            moho_undulation,
            window_size=window_size,
            shift_min=shift_min,
            shift_max=shift_max,
            shift_step=shift_step,
            Te_range=(Te_min, Te_max),
        )

        # CREATE TE MAPS FOR EACH SHIFT
        print("\n" + "=" * 70)
        print(" CREATING TE MAPS")
        print("=" * 70)

        from matplotlib.ticker import FuncFormatter

        def format_func(value, tick_number):
            return f"{int(value)}"

        def create_te_map_figure(
            result,
            shift_dist,
            X_topo,
            Y_topo,
            Te_min,
            Te_max,
            output_folder,
            output_folder_3d,
            filename_prefix,
        ):
            """Create and save a Te map figure for a given shift distance"""
            # 2D figure
            fig_te = plt.figure(figsize=(12, 10))
            ax_te = fig_te.add_subplot(111)

            cmap_te = plt.colormaps.get_cmap("jet")

            # Calculate extent
            x_centers_m = result["x_centers"] + X_topo[0, 0]
            y_centers_m = result["y_centers"] + Y_topo[0, 0]
            extent = [
                x_centers_m.min(),
                x_centers_m.max(),
                y_centers_m.min(),
                y_centers_m.max(),
            ]

            # Prepare Te data - fill NaN with minimum
            Te_data = result["Te_map"] / 1000
            Te_data_filled = Te_data.copy()
            Te_data_filled[np.isnan(Te_data_filled)] = Te_min / 1000

            # Plot
            im_te = ax_te.imshow(
                Te_data_filled,
                extent=extent,
                cmap=cmap_te,
                origin="lower",
                aspect="equal",
                vmin=Te_min / 1000,
                vmax=Te_max / 1000,
                interpolation="nearest",
            )

            shift_km = int(shift_dist / 1000)
            ax_te.set_title(
                f"Elastic Thickness (Te) - Window: {result['window_size'] / 1000:.0f} km, Shift: {shift_km} km",
                fontsize=14,
                fontweight="bold",
            )
            ax_te.set_xlabel("X (m)", fontsize=12)
            ax_te.set_ylabel("Y (m)", fontsize=12)
            ax_te.xaxis.set_major_formatter(FuncFormatter(format_func))
            ax_te.yaxis.set_major_formatter(FuncFormatter(format_func))

            cbar_te = plt.colorbar(
                im_te, ax=ax_te, label="Te (km)", fraction=0.046, pad=0.04
            )
            cbar_te.formatter.set_useOffset(False)

            plt.tight_layout()

            # Ensure output folder exists (create with full path)
            output_folder = os.path.abspath(output_folder)
            os.makedirs(output_folder, exist_ok=True)
            
            # Save figure with prefix
            base_filename = f"te_map_shift_{shift_km}km.png"
            filename = f"{filename_prefix}{base_filename}" if filename_prefix else base_filename
            te_map_path = os.path.join(output_folder, filename)
            # Normalize path to handle any path issues
            te_map_path = os.path.normpath(te_map_path)
            
            try:
                fig_te.savefig(te_map_path, dpi=300, bbox_inches="tight")
                print(f" Saved: {te_map_path}")
            except (FileNotFoundError, OSError) as e:
                print(f" Error saving {te_map_path}: {e}")
                # Ensure parent directory exists
                parent_dir = os.path.dirname(te_map_path)
                if not os.path.exists(parent_dir):
                    os.makedirs(parent_dir, exist_ok=True)
                fig_te.savefig(te_map_path, dpi=300, bbox_inches="tight")
                print(f" Saved (retry): {te_map_path}")

            plt.close(fig_te)  # Close to free memory

            # 3D Te map
            fig_te_3d = plt.figure(figsize=(12, 10))
            ax_te_3d = fig_te_3d.add_subplot(111, projection="3d")

            Xg, Yg = np.meshgrid(x_centers_m, y_centers_m)
            surf_te = ax_te_3d.plot_surface(
                Xg,
                Yg,
                Te_data_filled,
                cmap=cmap_te,
                linewidth=0,
                antialiased=True,
                vmin=Te_min / 1000,
                vmax=Te_max / 1000,
            )

            ax_te_3d.set_title(
                f"Elastic Thickness (Te) 3D - Window: {result['window_size'] / 1000:.0f} km, Shift: {shift_km} km",
                fontsize=14,
                fontweight="bold",
            )
            ax_te_3d.set_xlabel("X (m)", fontsize=12)
            ax_te_3d.set_ylabel("Y (m)", fontsize=12)
            ax_te_3d.set_zlabel("Te (km)", fontsize=12)

            fig_te_3d.colorbar(surf_te, shrink=0.5, aspect=10, pad=0.1, label="Te (km)")
            plt.tight_layout()

            # Ensure 3D output folder exists (create with full path)
            output_folder_3d = os.path.abspath(output_folder_3d)
            os.makedirs(output_folder_3d, exist_ok=True)
            
            # Save 3D figure with prefix
            base_filename_3d = f"te_map_shift_{shift_km}km_3d.png"
            filename_3d = f"{filename_prefix}{base_filename_3d}" if filename_prefix else base_filename_3d
            te_map_3d_path = os.path.join(output_folder_3d, filename_3d)
            # Normalize path to handle any path issues
            te_map_3d_path = os.path.normpath(te_map_3d_path)
            
            try:
                fig_te_3d.savefig(te_map_3d_path, dpi=300, bbox_inches="tight")
                print(f" Saved: {te_map_3d_path}")
            except (FileNotFoundError, OSError) as e:
                print(f" Error saving {te_map_3d_path}: {e}")
                # Try creating parent directory explicitly
                parent_dir = os.path.dirname(te_map_3d_path)
                os.makedirs(parent_dir, exist_ok=True)
                fig_te_3d.savefig(te_map_3d_path, dpi=300, bbox_inches="tight")
                print(f" Saved (retry): {te_map_3d_path}")

            plt.close(fig_te_3d)

        def create_rms_map_figure(
            result,
            shift_dist,
            X_topo,
            Y_topo,
            output_folder,
            output_folder_3d,
            filename_prefix,
        ):
            """Create and save an RMS misfit map figure for a given shift distance"""
            # 2D figure
            fig_rms = plt.figure(figsize=(12, 10))
            ax_rms = fig_rms.add_subplot(111)

            # Use a colormap that shows low RMS (good fit) in green and high RMS (poor fit) in red
            cmap_rms = plt.colormaps.get_cmap("RdYlGn_r")  # Reversed: green=low RMS, red=high RMS

            # Calculate extent
            x_centers_m = result["x_centers"] + X_topo[0, 0]
            y_centers_m = result["y_centers"] + Y_topo[0, 0]
            extent = [
                x_centers_m.min(),
                x_centers_m.max(),
                y_centers_m.min(),
                y_centers_m.max(),
            ]

            # Prepare RMS data - convert to km for display
            rms_data = np.abs(result["rms_map"]) / 1000  # Convert m to km, ensure positive
            rms_data_valid = rms_data[~np.isnan(rms_data)]
            
            if len(rms_data_valid) > 0:
                vmin_rms = np.nanpercentile(rms_data, 5)  # 5th percentile
                vmax_rms = np.nanpercentile(rms_data, 95)  # 95th percentile
            else:
                vmin_rms = 0
                vmax_rms = 1

            # Plot
            im_rms = ax_rms.imshow(
                rms_data,
                extent=extent,
                cmap=cmap_rms,
                origin="lower",
                aspect="equal",
                vmin=vmin_rms,
                vmax=vmax_rms,
                interpolation="nearest",
            )

            shift_km = int(shift_dist / 1000)
            ax_rms.set_title(
                f"RMS Misfit - Window: {result['window_size'] / 1000:.0f} km, Shift: {shift_km} km",
                fontsize=14,
                fontweight="bold",
            )
            ax_rms.set_xlabel("X (m)", fontsize=12)
            ax_rms.set_ylabel("Y (m)", fontsize=12)
            ax_rms.xaxis.set_major_formatter(FuncFormatter(format_func))
            ax_rms.yaxis.set_major_formatter(FuncFormatter(format_func))

            cbar_rms = plt.colorbar(
                im_rms, ax=ax_rms, label="RMS Misfit (km)", fraction=0.046, pad=0.04
            )
            cbar_rms.formatter.set_useOffset(False)

            # Add statistics text
            if len(rms_data_valid) > 0:
                mean_rms = np.nanmean(rms_data)
                min_rms = np.nanmin(rms_data)
                max_rms = np.nanmax(rms_data)
                stats_text = (
                    f"Mean: {mean_rms:.3f} km\n"
                    f"Min: {min_rms:.3f} km\n"
                    f"Max: {max_rms:.3f} km"
                )
                ax_rms.text(
                    0.02, 0.98, stats_text,
                    transform=ax_rms.transAxes,
                    fontsize=10,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )

            plt.tight_layout()

            # Ensure output folder exists
            output_folder = os.path.abspath(output_folder)
            os.makedirs(output_folder, exist_ok=True)
            
            # Save figure with prefix
            base_filename = f"rms_map_shift_{shift_km}km.png"
            filename = f"{filename_prefix}{base_filename}" if filename_prefix else base_filename
            rms_map_path = os.path.join(output_folder, filename)
            rms_map_path = os.path.normpath(rms_map_path)
            
            try:
                fig_rms.savefig(rms_map_path, dpi=300, bbox_inches="tight")
                print(f" Saved: {rms_map_path}")
            except (FileNotFoundError, OSError) as e:
                print(f" Error saving {rms_map_path}: {e}")
                parent_dir = os.path.dirname(rms_map_path)
                if not os.path.exists(parent_dir):
                    os.makedirs(parent_dir, exist_ok=True)
                fig_rms.savefig(rms_map_path, dpi=300, bbox_inches="tight")
                print(f" Saved (retry): {rms_map_path}")

            plt.close(fig_rms)  # Close to free memory

            # 3D RMS map
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  - needed for 3D plots
            fig_rms_3d = plt.figure(figsize=(12, 10))
            ax_rms_3d = fig_rms_3d.add_subplot(111, projection="3d")

            Xg, Yg = np.meshgrid(x_centers_m, y_centers_m)
            surf_rms = ax_rms_3d.plot_surface(
                Xg,
                Yg,
                rms_data,
                cmap=cmap_rms,
                linewidth=0,
                antialiased=True,
                vmin=vmin_rms,
                vmax=vmax_rms,
            )

            ax_rms_3d.set_title(
                f"RMS Misfit 3D - Window: {result['window_size'] / 1000:.0f} km, Shift: {shift_km} km",
                fontsize=14,
                fontweight="bold",
            )
            ax_rms_3d.set_xlabel("X (m)", fontsize=12)
            ax_rms_3d.set_ylabel("Y (m)", fontsize=12)
            ax_rms_3d.set_zlabel("RMS Misfit (km)", fontsize=12)

            fig_rms_3d.colorbar(surf_rms, shrink=0.5, aspect=10, pad=0.1, label="RMS Misfit (km)")
            plt.tight_layout()

            # Ensure 3D output folder exists
            output_folder_3d = os.path.abspath(output_folder_3d)
            os.makedirs(output_folder_3d, exist_ok=True)
            
            # Save 3D figure with prefix
            base_filename_3d = f"rms_map_shift_{shift_km}km_3d.png"
            filename_3d = f"{filename_prefix}{base_filename_3d}" if filename_prefix else base_filename_3d
            rms_map_3d_path = os.path.join(output_folder_3d, filename_3d)
            rms_map_3d_path = os.path.normpath(rms_map_3d_path)
            
            try:
                fig_rms_3d.savefig(rms_map_3d_path, dpi=300, bbox_inches="tight")
                print(f" Saved: {rms_map_3d_path}")
            except (FileNotFoundError, OSError) as e:
                print(f" Error saving {rms_map_3d_path}: {e}")
                parent_dir = os.path.dirname(rms_map_3d_path)
                os.makedirs(parent_dir, exist_ok=True)
                fig_rms_3d.savefig(rms_map_3d_path, dpi=300, bbox_inches="tight")
                print(f" Saved (retry): {rms_map_3d_path}")

            plt.close(fig_rms_3d)

        # Create Te map figures for each shift distance
        for shift_dist, result in mw_results_dict.items():
            create_te_map_figure(
                result,
                shift_dist,
                X_topo,
                Y_topo,
                Te_min,
                Te_max,
                output_folder,
                output_folder_3d,
                filename_prefix,
            )
        
        # CREATE RMS MAPS FOR EACH SHIFT
        print("\n" + "=" * 70)
        print(" CREATING RMS MAPS")
        print("=" * 70)
        
        # Create RMS map figures for each shift distance
        for shift_dist, result in mw_results_dict.items():
            create_rms_map_figure(
                result,
                shift_dist,
                X_topo,
                Y_topo,
                output_folder,
                output_folder_3d,
                filename_prefix,
            )

        # SAVE DATA
        print("\n" + "=" * 70)
        print(" SAVING DATA")
        print("=" * 70)

        save_dict = {
            "topography": topography,
            "topography_anomaly": topo_anom,
            "moho_depth": moho_depth,
            "moho_undulation": moho_undulation,
            "X": X_topo,
            "Y": Y_topo,
            "stats": stats,
        }

        if perform_global == "y":
            save_dict["global_Te"] = global_result["Te_best"]
            save_dict["global_rms"] = global_result["rms_best"]
            if "Te_values" in locals() and "rms_values" in locals():
                save_dict["sensitivity_Te_values"] = Te_values
                save_dict["sensitivity_rms_values"] = rms_values

        for shift_dist, res in mw_results_dict.items():
            km = int(shift_dist / 1000)
            save_dict[f"Te_map_shift_{km}km"] = res["Te_map"]
            save_dict[f"rms_map_shift_{km}km"] = res["rms_map"]
            # Save center coordinates for Te map interpolation
            save_dict[f"x_centers_shift_{km}km"] = res["x_centers"]
            save_dict[f"y_centers_shift_{km}km"] = res["y_centers"]

<<<<<<< HEAD
        # Additionally, export key outputs as Surfer DSAA grids for visualization in Surfer
        print("\n Exporting key grids to Surfer (.grd) format...")

        # Input grids
        write_surfer_grd(
            os.path.join(output_folder, f"{filename_prefix}topography.grd"),
            X_topo,
            Y_topo,
            topography,
        )
        write_surfer_grd(
            os.path.join(output_folder, f"{filename_prefix}moho_depth.grd"),
            X_moho,
            Y_moho,
            moho_depth,
        )

        # Anomaly / undulation grids
        write_surfer_grd(
            os.path.join(output_folder, f"{filename_prefix}topography_anomaly.grd"),
            X_topo,
            Y_topo,
            topo_anom,
        )
        write_surfer_grd(
            os.path.join(output_folder, f"{filename_prefix}moho_undulation.grd"),
            X_moho,
            Y_moho,
            moho_undulation,
        )

        # Moving-window Te and RMS maps for each shift (if present)
        for shift_dist, res in mw_results_dict.items():
            shift_km = int(shift_dist / 1000)

            # Build coordinate grids for window centers
            x_centers_m = res["x_centers"] + X_topo[0, 0]
            y_centers_m = res["y_centers"] + Y_topo[0, 0]
            Xc, Yc = np.meshgrid(x_centers_m, y_centers_m)

            # Te in km (optional, easier to read in Surfer)
            te_grid_km = res["Te_map"] / 1000.0
            write_surfer_grd(
                os.path.join(
                    output_folder, f"{filename_prefix}Te_map_shift_{shift_km}km.grd"
                ),
                Xc,
                Yc,
                te_grid_km,
            )

            # RMS misfit in km, absolute value
            rms_grid_km = np.abs(res["rms_map"]) / 1000.0
            write_surfer_grd(
                os.path.join(
                    output_folder, f"{filename_prefix}rms_map_shift_{shift_km}km.grd"
                ),
                Xc,
                Yc,
                rms_grid_km,
            )

=======
>>>>>>> e164d71a87412ccae8297068cbd15b9f29754aad
        # Save with prefix
        npz_filename = f"{filename_prefix}inversion_results.npz" if filename_prefix else "inversion_results.npz"
        np.savez(os.path.join(output_folder, npz_filename), **save_dict)

        # Save Surfer .grd outputs (same grid as input for topo/Moho; center grid for Te/rms)
        print(" Writing Surfer .grd files...")
        grd_prefix = filename_prefix if filename_prefix else ""
        write_surfer_grd(
            os.path.join(output_folder, f"{grd_prefix}topography.grd"),
            X_topo, Y_topo, topography,
        )
        write_surfer_grd(
            os.path.join(output_folder, f"{grd_prefix}topography_anomaly.grd"),
            X_topo, Y_topo, topo_anom,
        )
        write_surfer_grd(
            os.path.join(output_folder, f"{grd_prefix}moho_depth.grd"),
            X_topo, Y_topo, moho_depth,
        )
        write_surfer_grd(
            os.path.join(output_folder, f"{grd_prefix}moho_undulation.grd"),
            X_topo, Y_topo, moho_undulation,
        )
        if perform_mw == "y" and "mw_results_dict" in locals():
            for shift_dist, res in mw_results_dict.items():
                km = int(shift_dist / 1000)
                x_abs = res["x_centers"] + X_topo[0, 0]
                y_abs = res["y_centers"] + Y_topo[0, 0]
                write_surfer_grd(
                    os.path.join(output_folder, f"{grd_prefix}Te_map_shift_{km}km.grd"),
                    x_abs, y_abs, res["Te_map"],
                )
                write_surfer_grd(
                    os.path.join(output_folder, f"{grd_prefix}rms_map_shift_{km}km.grd"),
                    x_abs, y_abs, res["rms_map"],
                )
        print(" Surfer .grd files saved.")

    # PLOT FINAL RESULTS
    # Plot results whether or not moving window was performed
    print("\n" + "=" * 70)
    print(" PLOTTING FINAL RESULTS")
    print("=" * 70)
    
    plot_final_results(
        X_topo,
        Y_topo,
        topography,
        moho_depth,
        moho_undulation,
        topo_anom,
        output_folder,
        dx,
        dy,
        use_gravity=(use_gravity == "y"),
        filename_prefix=filename_prefix,
    )

    print(f"\n Results saved to: {output_folder}")
    print("  - inversion_results.npz")
    print("  - Surfer .grd: topography, topography_anomaly, moho_depth, moho_undulation" + (
          " (+ Te/rms maps per shift)" if perform_mw == "y" and "mw_results_dict" in locals() else ""
    ))
    print("  - input_data.png")
    print("  - 3D visualizations in subfolder: 3D")
    # List all Te and RMS map figures (if moving window was performed)
    if perform_mw == "y" and "mw_results_dict" in locals():
        for shift_dist in mw_results_dict.keys():
            shift_km = int(shift_dist / 1000)
            print(f"  - te_map_shift_{shift_km}km.png")
            print(f"  - rms_map_shift_{shift_km}km.png")
    if perform_global == "y" and "fig_sens" in locals():
        print("  - sensitivity_analysis.png")

    # COMPLETION
    print("\n" + "=" * 70)
    print(" ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nAll outputs saved in: {output_folder}")
    
    # Close terminal logging and restore stdout
    if 'tee' in locals():
        tee.close()
        print(f"\n Terminal output saved to: {log_file_path}")

    plt.show()


if __name__ == "__main__":
    main()
