"""
Main script for 3D Elastic Thickness Inversion
Based on Braitenberg Convolution Method, extended to 3D volumetric data
"""

import numpy as np
import os
import sys
from datetime import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Add parent directory to path to import from Code
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from data_loader_3d import (
    read_3d_volume_from_slices,
    prepare_3d_data_for_inversion,
    check_3d_grid_compatibility,
)
from elastic_thickness_inversion_3d import ElasticThicknessInversion3D
from moving_window_analysis_3d import MovingWindowAnalysis3D


def create_output_folder():
    """Create timestamped output folder"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"Output_3D_{timestamp}"
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), folder_name)
    os.makedirs(output_path, exist_ok=True)
    print(f"\nCreated output folder: {output_path}")
    return output_path


def plot_3d_input_data(X, Y, Z, topography_3d, moho_3d, topo_anom_3d, moho_undulation_3d, output_folder):
    """Plot 3D input data for verification"""
    from matplotlib.ticker import FuncFormatter

    nz, ny, nx = topography_3d.shape
    
    # Plot representative slices
    fig = plt.figure(figsize=(16, 12))
    
    # Middle depth slice
    z_mid = nz // 2
    
    # Original topography - middle slice
    ax1 = fig.add_subplot(2, 3, 1)
    x_m = X[z_mid, 0, :]
    y_m = Y[z_mid, :, 0]
    extent = [x_m.min(), x_m.max(), y_m.min(), y_m.max()]
    
    im1 = ax1.imshow(
        topography_3d[z_mid, :, :], extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    ax1.set_title(f"Topography - Depth Slice {z_mid}", fontsize=12, fontweight="bold")
    ax1.set_xlabel("X (m)", fontsize=11)
    ax1.set_ylabel("Y (m)", fontsize=11)
    plt.colorbar(im1, ax=ax1)
    
    # Topography anomaly
    ax2 = fig.add_subplot(2, 3, 2)
    im2 = ax2.imshow(
        topo_anom_3d[z_mid, :, :], extent=extent, cmap="RdBu_r", origin="lower", aspect="equal"
    )
    ax2.set_title(f"Topography Anomaly - Depth Slice {z_mid}", fontsize=12, fontweight="bold")
    ax2.set_xlabel("X (m)", fontsize=11)
    ax2.set_ylabel("Y (m)", fontsize=11)
    plt.colorbar(im2, ax=ax2)
    
    # Original Moho depth
    ax3 = fig.add_subplot(2, 3, 3)
    im3 = ax3.imshow(
        moho_3d[z_mid, :, :], extent=extent, cmap="jet", origin="lower", aspect="equal"
    )
    ax3.set_title(f"Moho Depth - Depth Slice {z_mid}", fontsize=12, fontweight="bold")
    ax3.set_xlabel("X (m)", fontsize=11)
    ax3.set_ylabel("Y (m)", fontsize=11)
    plt.colorbar(im3, ax=ax3)
    
    # Moho undulation
    ax4 = fig.add_subplot(2, 3, 4)
    im4 = ax4.imshow(
        moho_undulation_3d[z_mid, :, :], extent=extent, cmap="RdBu_r", origin="lower", aspect="equal"
    )
    ax4.set_title(f"Moho Undulation - Depth Slice {z_mid}", fontsize=12, fontweight="bold")
    ax4.set_xlabel("X (m)", fontsize=11)
    ax4.set_ylabel("Y (m)", fontsize=11)
    plt.colorbar(im4, ax=ax4)
    
    # 3D surface plot of topography
    ax5 = fig.add_subplot(2, 3, 5, projection='3d')
    x_sample = X[z_mid, 0, ::max(1, nx//50)]
    y_sample = Y[z_mid, ::max(1, ny//50), 0]
    X_sample, Y_sample = np.meshgrid(x_sample / 1000, y_sample / 1000)
    Z_sample = topography_3d[z_mid, ::max(1, ny//50), ::max(1, nx//50)]
    surf = ax5.plot_surface(X_sample, Y_sample, Z_sample, cmap='jet', alpha=0.8)
    ax5.set_xlabel("X (km)", fontsize=10)
    ax5.set_ylabel("Y (km)", fontsize=10)
    ax5.set_zlabel("Topography (m)", fontsize=10)
    ax5.set_title("3D Topography Surface", fontsize=12, fontweight="bold")
    
    # 3D surface plot of Moho
    ax6 = fig.add_subplot(2, 3, 6, projection='3d')
    Z_moho = moho_3d[z_mid, ::max(1, ny//50), ::max(1, nx//50)]
    surf2 = ax6.plot_surface(X_sample, Y_sample, Z_moho, cmap='jet', alpha=0.8)
    ax6.set_xlabel("X (km)", fontsize=10)
    ax6.set_ylabel("Y (km)", fontsize=10)
    ax6.set_zlabel("Moho Depth (m)", fontsize=10)
    ax6.set_title("3D Moho Surface", fontsize=12, fontweight="bold")
    
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_folder, "input_data_3d.png"), dpi=300, bbox_inches="tight"
    )
    plt.show(block=False)
    
    return fig


def main():
    """Main execution function for 3D inversion"""

    print("\n" + "=" * 70)
    print(" 3D ELASTIC THICKNESS INVERSION - BRAITENBERG METHOD")
    print("=" * 70)

    # ============== DATA SELECTION ==============
    print("\nFor 3D data, you need to provide multiple 2D slices or a 3D volume.")
    print("Currently, the code supports loading multiple 2D .grd files as depth slices.")
    
    data_type = input(
        "\nSelect data input method:\n  1. Multiple 2D slices (file pattern)\n  2. Single 2D file (will be treated as single slice)\nChoice: "
    )
    
    if data_type == "1":
        print("\n→ Using MULTIPLE 2D SLICES")
        print("Example pattern: 'path/to/topo_slice_*.grd' or list of files")
        file_pattern = input("Enter file pattern or path to first file: ").strip()
        
        # Try to find files
        import glob
        if '*' in file_pattern or '?' in file_pattern:
            topo_files = sorted(glob.glob(file_pattern))
        else:
            # Assume pattern like "topo_slice_0.grd", "topo_slice_1.grd", etc.
            base_path = os.path.dirname(file_pattern) if os.path.dirname(file_pattern) else '.'
            base_name = os.path.basename(file_pattern)
            # Try to find similar files
            topo_files = sorted(glob.glob(os.path.join(base_path, base_name.replace('0', '*'))))
            if not topo_files:
                topo_files = [file_pattern]
        
        if not topo_files:
            print(f"\n❌ ERROR: No files found matching pattern: {file_pattern}")
            return
        
        print(f"Found {len(topo_files)} topography files")
        
        # For Moho files, try similar pattern
        moho_pattern = input("Enter Moho file pattern (or press Enter to use same pattern): ").strip()
        if not moho_pattern:
            moho_pattern = file_pattern.replace('topo', 'moho').replace('Topo', 'Moho')
        
        if '*' in moho_pattern or '?' in moho_pattern:
            moho_files = sorted(glob.glob(moho_pattern))
        else:
            moho_files = [moho_pattern]
        
        if not moho_files or len(moho_files) != len(topo_files):
            print(f"⚠️  Warning: Found {len(moho_files)} Moho files, expected {len(topo_files)}")
            print("Using first Moho file for all slices...")
            if moho_files:
                moho_files = [moho_files[0]] * len(topo_files)
            else:
                print("❌ ERROR: No Moho files found!")
                return
        
        # Read 3D volumes
        try:
            X_topo, Y_topo, Z_topo, topography_3d, dx, dy, dz = read_3d_volume_from_slices(topo_files)
            X_moho, Y_moho, Z_moho, moho_3d, dx_m, dy_m, dz_m = read_3d_volume_from_slices(moho_files)
        except Exception as e:
            print(f"\n❌ ERROR loading files: {e}")
            return
        
    elif data_type == "2":
        print("\n→ Using SINGLE 2D FILE (treated as single depth slice)")
        from Code.data_loader import read_surfer_grd
        
        topography_file = input("Enter topography file path: ").strip()
        moho_file = input("Enter Moho file path: ").strip()
        
        if not os.path.exists(topography_file) or not os.path.exists(moho_file):
            print("\n❌ ERROR: Input files not found!")
            return
        
        # Read 2D data
        X_2d, Y_2d, topography_2d, dx, dy, nx, ny, xmin, xmax, ymin, ymax = read_surfer_grd(topography_file)
        X_m_2d, Y_m_2d, moho_2d, dx_m, dy_m, nx_m, ny_m, xmin_m, xmax_m, ymin_m, ymax_m = read_surfer_grd(moho_file)
        
        # Convert to 3D (single slice)
        topography_3d = topography_2d[np.newaxis, :, :]
        moho_3d = moho_2d[np.newaxis, :, :]
        
        # Create 3D coordinate arrays
        X_topo = np.zeros((1, ny, nx))
        Y_topo = np.zeros((1, ny, nx))
        Z_topo = np.zeros((1, ny, nx))
        X_topo[0, :, :] = X_2d
        Y_topo[0, :, :] = Y_2d
        Z_topo[0, :, :] = 0.0
        
        X_moho = np.zeros((1, ny_m, nx_m))
        Y_moho = np.zeros((1, ny_m, nx_m))
        Z_moho = np.zeros((1, ny_m, nx_m))
        X_moho[0, :, :] = X_m_2d
        Y_moho[0, :, :] = Y_m_2d
        Z_moho[0, :, :] = 0.0
        
        dz = 1000.0  # Default depth spacing
        dz_m = 1000.0
        
    else:
        print("Invalid option selected.")
        return

    # Check grid compatibility
    compatible, message = check_3d_grid_compatibility(X_topo, Y_topo, Z_topo, X_moho, Y_moho, Z_moho)
    if not compatible:
        print(f"\n⚠️  WARNING: {message}")
        print("Proceeding with caution...")
    else:
        print(f"\n✓ {message}")

    # Create output folder
    output_folder = create_output_folder()

    # ============== PREPARE DATA ==============
    print("\n" + "=" * 70)
    print(" PREPARING 3D DATA FOR INVERSION")
    print("=" * 70)

    topo_anom_3d, moho_undulation_3d, stats = prepare_3d_data_for_inversion(
        topography_3d, moho_3d, apply_taper_flag=True, taper_alpha=0.1
    )

    print("\n3D Data Statistics:")
    print(f"  Shape: {stats['shape']}")
    print(f"  Topography mean: {stats['topo_mean']:.1f} m")
    print(f"  Topography std (original): {stats['topo_std']:.1f} m")
    print(f"  Topography std (anomaly): {stats['topo_anom_std']:.1f} m")
    print(f"  Moho mean depth: {stats['moho_mean']:.1f} m")
    print(f"  Moho std (original): {stats['moho_std']:.1f} m")
    print(f"  Moho undulation std: {stats['moho_undulation_std']:.1f} m")
    print(f"  Taper applied: {stats['tapered']}")

    # ============== PHYSICAL PARAMETERS ==============
    print("\n" + "=" * 70)
    print(" PHYSICAL PARAMETERS (MARS)")
    print("=" * 70)
    print("  Crustal density (ρ_crust): 2900 kg/m³")
    print("  Mantle density (ρ_mantle): 3500 kg/m³")
    print("  Density contrast (Δρ): 600 kg/m³")
    print("  Gravity (g): 3.72 m/s²")
    print("  Young's modulus (E): 1.0×10¹¹ Pa")
    print("  Poisson's ratio (ν): 0.25")

    # Plot input data
    plot_3d_input_data(
        X_topo, Y_topo, Z_topo, topography_3d, moho_3d, topo_anom_3d, moho_undulation_3d, output_folder
    )

    # ============== GLOBAL INVERSION ==============
    perform_global = input("\nPerform global (full-region) 3D inversion? (y/n): ").lower()

    if perform_global == "y":
        print("\n" + "=" * 70)
        print(" GLOBAL 3D INVERSION")
        print("=" * 70)

        Te_range = (1000, 100000)
        print(f"Te search range: {Te_range[0] / 1000:.0f}-{Te_range[1] / 1000:.0f} km")

        inverter = ElasticThicknessInversion3D(
            dx=dx, dy=dy, dz=dz, rho_load=2900, rho_m=3500, rho_infill=2900, g=3.72
        )

        global_result = inverter.invert_elastic_thickness_3d(
            topo_anom_3d, moho_undulation_3d, Te_range=Te_range, method="bounded"
        )

        print(f"\n✓ Global Te estimate: {global_result['Te_best'] / 1000:.2f} km")
        print(f"  RMS misfit: {global_result['rms_best']:.2f} m")

    # ============== MOVING WINDOW ANALYSIS ==============
    perform_mw = input("\nPerform 3D moving window analysis? (y/n): ").lower()

    if perform_mw == "y":
        print("\n" + "=" * 70)
        print(" 3D MOVING WINDOW ANALYSIS")
        print("=" * 70)

        try:
            window_size = (
                float(input("Window size in X-Y plane (km) [default 1000]: ") or "1000") * 1000
            )
            window_depth_input = input("Window depth (km) [press Enter for full depth]: ").strip()
            window_depth = float(window_depth_input) * 1000 if window_depth_input else None
            
            shift_min = float(input("Minimum shift in X-Y plane (km) [default 20]: ") or "20") * 1000
            shift_max = float(input("Maximum shift in X-Y plane (km) [default 80]: ") or "80") * 1000
            shift_step = float(input("Shift step (km) [default 20]: ") or "20") * 1000
            Te_min = float(input("Minimum Te (km) [default 5]: ") or "5") * 1000
            Te_max = float(input("Maximum Te (km) [default 80]: ") or "80") * 1000
        except:
            print("Using default parameters...")
            window_size = 1_000_000
            window_depth = None
            shift_min, shift_max, shift_step = 20_000, 80_000, 20_000
            Te_min, Te_max = 5_000, 80_000

        print("\n→ Analysis Parameters:")
        print(f"  Window size (X-Y): {window_size / 1000:.0f} km")
        if window_depth:
            print(f"  Window depth: {window_depth / 1000:.0f} km")
        else:
            print(f"  Window depth: Full depth ({topography_3d.shape[0]} slices)")
        print(
            f"  Shift range: {shift_min / 1000:.0f}-{shift_max / 1000:.0f} km (step {shift_step / 1000:.0f} km)"
        )
        print(f"  Te search range: {Te_min / 1000:.0f}-{Te_max / 1000:.0f} km")

        mw_analyzer = MovingWindowAnalysis3D(
            dx=dx, dy=dy, dz=dz, rho_load=2900, rho_m=3500, rho_infill=2900, g=3.72
        )

        mw_results_dict = mw_analyzer.analyze_multiple_shifts_3d(
            topo_anom_3d,
            moho_undulation_3d,
            window_size=window_size,
            window_depth=window_depth,
            shift_min=shift_min,
            shift_max=shift_max,
            shift_step=shift_step,
            Te_range=(Te_min, Te_max),
        )

        # ============== CREATE TE MAPS FOR EACH SHIFT ==============
        print("\n" + "=" * 70)
        print(" CREATING 3D TE MAPS")
        print("=" * 70)

        from matplotlib.ticker import FuncFormatter

        def format_func(value, tick_number):
            return f"{int(value)}"

        def create_te_map_figure_3d(result, shift_dist, X_topo, Y_topo, Te_min, Te_max, output_folder):
            """Create and save a Te map figure for a given shift distance"""
            fig_te = plt.figure(figsize=(12, 10))
            ax_te = fig_te.add_subplot(111)

            cmap_te = plt.colormaps.get_cmap("jet")

            x_centers_m = result["x_centers"] + X_topo[0, 0, 0]
            y_centers_m = result["y_centers"] + Y_topo[0, 0, 0]
            extent = [
                x_centers_m.min(),
                x_centers_m.max(),
                y_centers_m.min(),
                y_centers_m.max(),
            ]

            Te_data = result["Te_map"] / 1000
            Te_data_filled = Te_data.copy()
            Te_data_filled[np.isnan(Te_data_filled)] = Te_min / 1000

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
                f"3D Elastic Thickness (Te) - Window: {result['window_size'] / 1000:.0f} km, Shift: {shift_km} km",
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

            te_map_path = os.path.join(output_folder, f"te_map_3d_shift_{shift_km}km.png")
            fig_te.savefig(te_map_path, dpi=300, bbox_inches="tight")
            print(f"✓ Saved: {te_map_path}")
            
            plt.close(fig_te)

        # Create a figure for each shift distance
        for shift_dist, result in mw_results_dict.items():
            create_te_map_figure_3d(
                result, shift_dist, X_topo, Y_topo, Te_min, Te_max, output_folder
            )

        # ============== SAVE DATA ==============
        print("\n" + "=" * 70)
        print(" SAVING DATA")
        print("=" * 70)

        save_dict = {
            "topography_3d": topography_3d,
            "topography_anomaly_3d": topo_anom_3d,
            "moho_depth_3d": moho_3d,
            "moho_undulation_3d": moho_undulation_3d,
            "X": X_topo,
            "Y": Y_topo,
            "Z": Z_topo,
            "stats": stats,
        }

        if perform_global == "y":
            save_dict["global_Te"] = global_result["Te_best"]
            save_dict["global_rms"] = global_result["rms_best"]

        for shift_dist, res in mw_results_dict.items():
            km = int(shift_dist / 1000)
            save_dict[f"Te_map_shift_{km}km"] = res["Te_map"]
            save_dict[f"rms_map_shift_{km}km"] = res["rms_map"]

        np.savez(os.path.join(output_folder, "inversion_results_3d.npz"), **save_dict)

        print(f"\n✓ Results saved to: {output_folder}")
        print("  - inversion_results_3d.npz")
        print("  - input_data_3d.png")
        for shift_dist in mw_results_dict.keys():
            shift_km = int(shift_dist / 1000)
            print(f"  - te_map_3d_shift_{shift_km}km.png")

    # ============== COMPLETION ==============
    print("\n" + "=" * 70)
    print(" 3D ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nAll outputs saved in: {output_folder}")

    plt.show()


if __name__ == "__main__":
    main()

