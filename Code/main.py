
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
    def __init__(self, file_path):
        self.terminal = sys.stdout
        self.log_file = open(file_path, 'w', encoding='utf-8')
    
    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush() 
    
    def flush(self):
        self.terminal.flush()
        self.log_file.flush()
    
    def close(self):
        if self.log_file:
            self.log_file.close()
            sys.stdout = self.terminal


def create_output_folder():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"Output_{timestamp}"
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), folder_name)
    os.makedirs(output_path, exist_ok=True)
    print(f"\nCreated output folder: {output_path}")
    return output_path


def create_filename_prefix(rho_load, ref_moho_km, window_size_km):
   
    density_gcm3 = rho_load / 1000.0
    
    prefix = f"{density_gcm3:.1f}_{ref_moho_km:.0f}_{window_size_km:.0f}_"
    return prefix


def plot_final_results(
    X, Y, topography, moho_depth,
    moho_undulation, topo_anom, output_folder, dx, dy, use_gravity=False, filename_prefix=""
):
    from matplotlib.ticker import FuncFormatter

    # Set toolbar to None to remove buttons
    plt.rcParams['toolbar'] = 'None'

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    axes = axes.flatten()

    # Use extent in KM so axis labels are readable (not millions of meters)
    x_m = X[0, :]
    y_m = Y[:, 0]
    extent_km = [x_m.min()/1000, x_m.max()/1000, y_m.min()/1000, y_m.max()/1000]

    def format_km(value, tick_number):
        return f"{int(value)}"

    # Data in km for display
    topo_km = topography / 1000.0
    moho_km = moho_depth / 1000.0
    tmin, tmax = np.nanmin(topo_km), np.nanmax(topo_km)
    if tmax - tmin < 0.001:
        tmin, tmax = tmin - 0.0005, tmax + 0.0005
    mmin, mmax = np.nanmin(moho_km), np.nanmax(moho_km)
    if mmax - mmin < 0.001:
        mmin, mmax = mmin - 0.0005, mmax + 0.0005

    # Panel 1: Topography (km)
    im1 = axes[0].imshow(
        topo_km, extent=extent_km, cmap="jet", origin="lower", aspect="equal",
        vmin=tmin, vmax=tmax
    )
    axes[0].set_title("Topography (km)", fontsize=14, fontweight="bold")
    axes[0].set_xlabel("X (km)", fontsize=12)
    axes[0].set_ylabel("Y (km)", fontsize=12)
    axes[0].xaxis.set_major_formatter(FuncFormatter(format_km))
    axes[0].yaxis.set_major_formatter(FuncFormatter(format_km))
    plt.colorbar(im1, ax=axes[0], label="Topography (km)")

    # Panel 2: Moho depth (km)
    moho_title = "Moho depth (km) [from gravity]" if use_gravity else "Moho depth (km)"
    im2 = axes[1].imshow(
        moho_km, extent=extent_km, cmap="jet", origin="lower", aspect="equal",
        vmin=mmin, vmax=mmax
    )
    axes[1].set_title(moho_title, fontsize=14, fontweight="bold")
    axes[1].set_xlabel("X (km)", fontsize=12)
    axes[1].set_ylabel("Y (km)", fontsize=12)
    axes[1].xaxis.set_major_formatter(FuncFormatter(format_km))
    axes[1].yaxis.set_major_formatter(FuncFormatter(format_km))
    plt.colorbar(im2, ax=axes[1], label="Moho depth (km)")

    plt.tight_layout()
    
    filename = f"{filename_prefix}input_data.png" if filename_prefix else "input_data.png"
    plt.savefig(
        os.path.join(output_folder, filename), dpi=300, bbox_inches="tight"
    )
    
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
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  - needed for 3D plots

    # Coordinates in km for display
    x_km = X[0, :] / 1000.0
    y_km = Y[:, 0] / 1000.0
    Xg, Yg = np.meshgrid(x_km, y_km)
    topo_anom_km = topo_anom / 1000.0
    moho_undulation_km = moho_undulation / 1000.0

    # Create the combined 3D plot
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Plot Topography surface (using topography anomaly for better visualization)
    # Use gist_earth colormap for topography
    surf_topo = ax.plot_surface(
        Xg,
        Yg,
        topo_anom_km,
        cmap="gist_earth",
        edgecolor="none",
        alpha=0.8,
        linewidth=0,
        antialiased=True,
    )

    surf_moho = ax.plot_surface(
        Xg,
        Yg,
        moho_undulation_km,
        cmap="RdBu_r",
        edgecolor="none",
        alpha=0.6,
        linewidth=0,
        antialiased=True,
    )
    z_zero = np.zeros_like(Xg)
    ax.plot_surface(
        Xg,
        Yg,
        z_zero,
        color="darkblue",
        alpha=0.3,
        linewidth=0,
    )

    moho_mean_km = np.mean(moho_undulation_km)
    z_moho_ref = np.full_like(Xg, moho_mean_km)
    ax.plot_surface(
        Xg,
        Yg,
        z_moho_ref,
        color="yellow",
        alpha=0.2,
        linewidth=0,
    )

    ax.set_title(
        "3D Visualization: Topography and Moho Depth (Airy Isostasy Model)",
        fontsize=15,
        fontweight="bold",
        pad=20,
    )

    ax.set_xlabel("X (km)", fontsize=12)
    ax.set_ylabel("Y (km)", fontsize=12)
    ax.set_zlabel("Elevation/Depth (km)", fontsize=12)

    z_min = min(moho_undulation_km.min(), -80)
    z_max = max(topo_anom_km.max(), 20)
    ax.set_zlim(z_min, z_max)

    cbar_topo = fig.colorbar(
        surf_topo,
        ax=ax,
        shrink=0.5,
        aspect=20,
        pad=0.05,
        location="left",
    )
    cbar_topo.set_label("Topography (km)", fontsize=11, rotation=90, labelpad=15)

    cbar_moho = fig.colorbar(
        surf_moho,
        ax=ax,
        shrink=0.5,
        aspect=20,
        pad=0.05,
        location="right",
    )
    cbar_moho.set_label("Moho Depth (km)", fontsize=11, rotation=90, labelpad=15)

    plt.tight_layout()

    filename = f"{filename_prefix}topography_moho_combined_3d.png" if filename_prefix else "topography_moho_combined_3d.png"
    save_path = os.path.join(output_folder_3d, filename)
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f" Saved: {save_path}")
    plt.close(fig)


def main():
    
    
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
        print("Moho input: will be predicted from gravity.")
    else:
        # Get Moho depth file path from user
        moho_file = input("Enter path to Moho depth file (.grd): ").strip()
        moho_file = moho_file.strip('"').strip("'")
        gravity_file = None
        print("Moho input: using provided Moho depth grid.")
    
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
            print(f"Reference Moho depth: {ref_moho_km:.1f} km.")
        except:
            ref_moho = None
            ref_moho_km = None
            print("Reference Moho depth: invalid input; using mean Moho depth.")
    else:
        ref_moho = None
        ref_moho_km = None
        print("Reference Moho depth: using mean Moho depth.")

    # Check file existence
    print("\nChecking input files.")
    print(f"  Topography: {topography_file}")
    
    if not os.path.exists(topography_file):
        print("\nError: topography file not found.")
        print(f"  {topography_file}")
        return
    
    if use_gravity == "y":
        print(f"  Gravity: {gravity_file}")
        if not os.path.exists(gravity_file):
            print("\nError: gravity file not found.")
            print(f"  {gravity_file}")
            return
        print("Input files look good.")
    else:
        print(f"  Moho: {moho_file}")
        if not os.path.exists(moho_file):
            print("\nError: Moho file not found.")
            print(f"  {moho_file}")
            return
        print("Input files look good.")

    # Create output folder
    output_folder = create_output_folder()
    output_folder_3d = os.path.join(output_folder, "3D")
    os.makedirs(output_folder_3d, exist_ok=True)
    
    # Set up terminal output logging to file
    log_file_path = os.path.join(output_folder, "terminal_output.txt")
    tee = TeeOutput(log_file_path)
    sys.stdout = tee
    print("\nLogging terminal output.")
    print(f"  Log file: {log_file_path}\n")

    # LOAD DATA
    print("Loading grids.")
     

    X_topo, Y_topo, topography, dx, dy, *_ = read_surfer_grd(topography_file)
    
    if use_gravity == "y":
        # Load gravity data and predict Moho
        X_grav, Y_grav, bouguer_gravity, *_ = read_surfer_grd(gravity_file)
        
        # Check grid compatibility
        compatible, message = check_grid_compatibility(X_topo, Y_topo, X_grav, Y_grav)
        if not compatible:
            print(f"Warning: {message}")
            print("Proceeding anyway.")
        else:
            print(message)
        
        # Predict Moho from gravity
        print("\nPredicting Moho from gravity.")
         
        
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
        
        print("Moho prediction complete.")
        print(f"  Moho depth range: {moho_depth.min()/1000:.1f} to {moho_depth.max()/1000:.1f} km")
    else:
        # Load observed Moho data
        X_moho, Y_moho, moho_depth, *_ = read_surfer_grd(moho_file)

        # Check grid compatibility
        compatible, message = check_grid_compatibility(X_topo, Y_topo, X_moho, Y_moho)
        if not compatible:
            print(f"Warning: {message}")
            print("Proceeding anyway.")
        else:
            print(message)

    # PREPARE DATA
    print("\nPreparing data for inversion.")
     
    # Apply taper and remove mean (or user-provided reference)
    if ref_moho is not None:
        # Use user-provided reference Moho
        topo_anom, moho_undulation, stats = prepare_data_for_inversion(
            topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1, ref_moho=ref_moho
        )
        if ref_moho_km is not None:
            print(f"  Using user-provided reference Moho: {ref_moho_km:.1f} km")
        else:
            print(f"  Using user-provided reference Moho: {ref_moho/1000:.1f} km")
    else:
        # Use mean as reference (default)
        topo_anom, moho_undulation, stats = prepare_data_for_inversion(
            topography, moho_depth, apply_taper_flag=True, taper_alpha=0.1
        )
        print(f"  Using mean Moho depth as reference: {stats['moho_mean']/1000:.1f} km")
    
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
    print(f"\nFilename prefix: {filename_prefix}")

    print("\nData Statistics:")
    print(f"  Topography mean: {stats['topo_mean']/1000:.1f} km")
    print(f"  Topography std (original): {stats['topo_std']/1000:.1f} km")
    print(f"  Topography std (anomaly): {stats['topo_anom_std']/1000:.1f} km")
    print(f"  Moho mean depth: {stats['moho_mean']/1000:.1f} km")
    print(f"  Moho std (original): {stats['moho_std']/1000:.1f} km")
    print(f"  Moho undulation std: {stats['moho_undulation_std']/1000:.1f} km")
    print(f"  Taper applied: {stats['tapered']}")

    # PHYSICAL PARAMETERS
    print("\nPhysical parameters (Mars).")
     
    print(f"  Crustal density (ρ_crust): {RHO_LOAD} kg/m³")
    print(f"  Mantle density (ρ_mantle): {RHO_MANTLE} kg/m³")
    print(f"  Density contrast (Δρ): {DENSITY_CONTRAST} kg/m³")
    print(f"  Gravity (g): {GRAVITY} m/s²")
    print(f"  Young's modulus (E): {YOUNGS_MODULUS:.1e} Pa")
    print(f"  Poisson's ratio (ν): {POISSONS_RATIO}")

    # DIAGNOSTIC CHECK
    print("\nQuick diagnostic check.")
     

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

    print(f"At max |topo| ({max_topo_val/1000:.1f} km):")
    print(f"  Moho undulation: {moho_at_max_topo/1000:.1f} km")

    if max_topo_val > 0 and moho_at_max_topo > 0:
        print("   Positive topo, Positive Moho deflection")
        print("   This means: Mountain  Moho goes UP (shallower)")
        print("   This is BACKWARDS! Moho should go DOWN under mountains")
        flip_sign = input("\nFlip Moho sign? (y/n): ").lower()
        if flip_sign == "y":
            print(" Flipping Moho undulation sign...")
            moho_undulation = -moho_undulation
    elif max_topo_val > 0 and moho_at_max_topo < 0:
        print("   Positive topo, Negative Moho deflection")
        print("   This means: Mountain  Moho goes DOWN (deeper)")
        print("   Sign convention is CORRECT")

    print("\n")

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
    # perform_global = input("\nRun a global (full-region) inversion? (y/n): ").lower()
    perform_global = "n"
    if perform_global == "y":
        print("\nRunning global inversion.")
         

        Te_range = (1000, 100000)
        print(f"Te search range: {Te_range[0] / 1000:.0f}-{Te_range[1] / 1000:.0f} km")

        global_result = inverter.invert_elastic_thickness(
            topo_anom, moho_undulation, Te_range=Te_range, method="bounded"
        )

        print(f"\nGlobal Te estimate: {global_result['Te_best'] / 1000:.2f} km")
        print(f"  RMS misfit: {global_result['rms_best']/1000:.2f} km")

        if true_Te is not None:
            error = abs(global_result["Te_best"] - true_Te) / true_Te * 100
            print(f"  True Te: {true_Te / 1000:.1f} km")
            print(f"  Error: {error:.1f}%")

        # Sensitivity analysis
        perform_sensitivity = input(
            "\nPerform sensitivity analysis (Te vs RMS curve)? (y/n): "
        ).lower()

        if perform_sensitivity == "y":
            print("\n Running sensitivity analysis...")
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
                print(f"\n  Confidence interval (RMS ≤ {threshold/1000:.1f} km):")
                print(
                    f"    Te range: {Te_conf.min() / 1000:.1f} - {Te_conf.max() / 1000:.1f} km"
                )
                print(f"    Spread: ±{(Te_conf.max() - Te_conf.min()) / 2000:.1f} km")

            rms_range = rms_values.max() - rms_values.min()
            if rms_range / min_rms > 2:
                print("\n   Solution appears well-constrained (clear minimum in RMS)")
            else:
                print("\n   Solution may be poorly constrained (shallow RMS minimum)")

    # WINDOW ANALYSIS: single window or moving window
    print("\nWindow analysis.")
     
    choice = input("Press 1 for single window, 2 for moving window: ").strip()
    if choice == "1":
        analysis_mode = "single"
    else:
        analysis_mode = "moving"
    perform_mw = "y" if analysis_mode == "moving" else "n"
    mw_results_dict = None
    single_window_result = None

    if analysis_mode == "single":
        # ---------- SINGLE WINDOW ANALYSIS ----------
        if use_gravity == "y":
            print("\n Single-window analysis using gravity-predicted Moho")
        else:
            print("\n Single-window analysis using observed Moho")
        print("\nSingle-window analysis.")
         

        x_m = X_topo[0, :]
        y_m = Y_topo[:, 0]
        xmin_m, xmax_m = x_m.min(), x_m.max()
        ymin_m, ymax_m = y_m.min(), y_m.max()
        xmin_km, xmax_km = xmin_m / 1000, xmax_m / 1000
        ymin_km, ymax_km = ymin_m / 1000, ymax_m / 1000
        nx, ny = topo_anom.shape[1], topo_anom.shape[0]

        # 1. Ask window size (km)
        try:
            window_size_km_sw = float(input("Window size (km) [default 1000]: ") or "1000")
        except Exception:
            window_size_km_sw = 1000.0
        window_size_m = window_size_km_sw * 1000
        half_m = window_size_m / 2.0

        # 2. Interactive plot: window follows cursor; left-click to select and run Te
        try:
            Te_min_sw = float(input("Minimum Te (km) [default 5]: ") or "5") * 1000
            Te_max_sw = float(input("Maximum Te (km) [default 80]: ") or "80") * 1000
        except Exception:
            Te_min_sw, Te_max_sw = 5_000, 80_000
        Te_range_sw = (Te_min_sw, Te_max_sw)

        from matplotlib.patches import Rectangle

        half_km = window_size_km_sw / 2.0
        # Clamp so full window stays inside grid
        x_center_km = (xmin_km + xmax_km) / 2.0
        y_center_km = (ymin_km + ymax_km) / 2.0
        selected = {}  # will store {'x': x_km, 'y': y_km} on click

        fig_sw, axes_sw = plt.subplots(1, 2, figsize=(14, 6))
        extent_km = [xmin_km, xmax_km, ymin_km, ymax_km]
        axes_sw[0].imshow(
            topography / 1000.0, extent=extent_km, cmap="jet", origin="lower", aspect="equal"
        )
        axes_sw[0].set_title("Topography (km)")
        axes_sw[0].set_xlabel("X (km)")
        axes_sw[0].set_ylabel("Y (km)")
        axes_sw[1].imshow(
            moho_depth / 1000.0, extent=extent_km, cmap="jet", origin="lower", aspect="equal"
        )
        axes_sw[1].set_title("Moho depth (km)")
        axes_sw[1].set_xlabel("X (km)")
        axes_sw[1].set_ylabel("Y (km)")
        fig_sw.suptitle(
            "Move cursor to position window  LEFT-CLICK to run Te inversion for that area",
            fontsize=11,
        )
        plt.tight_layout()

        # Rectangle (left, bottom, width, height) in km
        rects = []
        for ax in axes_sw:
            rect = Rectangle(
                (x_center_km - half_km, y_center_km - half_km),
                window_size_km_sw,
                window_size_km_sw,
                linewidth=2,
                edgecolor="lime",
                facecolor="none",
            )
            ax.add_patch(rect)
            rects.append(rect)

        def on_motion(event):
            if event.inaxes not in axes_sw or event.xdata is None or event.ydata is None:
                return
            xc = np.clip(event.xdata, xmin_km + half_km, xmax_km - half_km)
            yc = np.clip(event.ydata, ymin_km + half_km, ymax_km - half_km)
            for r in rects:
                r.set_xy((xc - half_km, yc - half_km))
            fig_sw.canvas.draw_idle()

        def on_click(event):
            if event.button != 1 or event.inaxes not in axes_sw:
                return
            if event.xdata is None or event.ydata is None:
                return
            selected["x"] = np.clip(event.xdata, xmin_km + half_km, xmax_km - half_km)
            selected["y"] = np.clip(event.ydata, ymin_km + half_km, ymax_km - half_km)
            plt.close(fig_sw)

        fig_sw.canvas.mpl_connect("motion_notify_event", on_motion)
        fig_sw.canvas.mpl_connect("button_press_event", on_click)
        print(
            f"\n  Grid extent: X = {xmin_km:.1f} to {xmax_km:.1f} km, "
            f"Y = {ymin_km:.1f} to {ymax_km:.1f} km"
        )
        print(f"  Window size: {window_size_km_sw:.0f} km. Move mouse, then LEFT-CLICK to run Te.")
        plt.show(block=True)

        if not selected:
            x_center_km = (xmin_km + xmax_km) / 2.0
            y_center_km = (ymin_km + ymax_km) / 2.0
            print("  No click recorded; using grid center.")
        else:
            x_center_km = selected["x"]
            y_center_km = selected["y"]
        print(f"  Selected window center: ({x_center_km:.1f}, {y_center_km:.1f}) km")

        x_center_m = x_center_km * 1000
        y_center_m = y_center_km * 1000
        x_min_m = x_center_m - half_m
        x_max_m = x_center_m + half_m
        y_min_m = y_center_m - half_m
        y_max_m = y_center_m + half_m
        x_min_m = max(xmin_m, min(x_min_m, xmax_m))
        x_max_m = max(xmin_m, min(x_max_m, xmax_m))
        y_min_m = max(ymin_m, min(y_min_m, ymax_m))
        y_max_m = max(ymin_m, min(y_max_m, ymax_m))

        i0 = int(np.argmin(np.abs(x_m - x_min_m)))
        i1 = int(np.argmin(np.abs(x_m - x_max_m)))
        j0 = int(np.argmin(np.abs(y_m - y_min_m)))
        j1 = int(np.argmin(np.abs(y_m - y_max_m)))
        if i0 > i1:
            i0, i1 = i1, i0
        if j0 > j1:
            j0, j1 = j1, j0
        i1 = min(i1 + 1, nx)
        j1 = min(j1 + 1, ny)
        window_bounds_str = f"center_{x_center_km:.0f}_{y_center_km:.0f}km_size_{window_size_km_sw:.0f}km"
        region_ny, region_nx = j1 - j0, i1 - i0
        region_size_km = window_size_km_sw  # selected region side in km

        # Mode for selected region:
        #   "1" → one Te for whole region (default, no prompt)
        #   "2" → Te map within region (not used here)
        single_region_mode = "1"
        single_window_result = None
        single_region_has_map = False

        if single_region_mode == "2":
            # ---------- Te map WITHIN selected region ----------
            print("\n Building Te map within selected region...")
            topo_region = topo_anom[j0:j1, i0:i1].copy()
            moho_region = moho_undulation[j0:j1, i0:i1].copy()
            # Default: analysis window = min(region_size/3, 400 km), shift = window/5
            analysis_window_km = min(region_size_km / 3.0, 400.0)
            try:
                aw = input(f"Analysis window size within region (km) [default {analysis_window_km:.0f}]: ").strip()
                if aw:
                    analysis_window_km = float(aw)
            except Exception:
                pass
            analysis_window_m = analysis_window_km * 1000
            # Single-window Te map: use default shift (no prompt — shift is for moving-window only)
            shift_km_default = max(analysis_window_km / 5.0, 10.0)
            shift_km_sw = shift_km_default
            shift_m_sw = shift_km_sw * 1000
            print(f"  Region size: {region_size_km:.0f} km; analysis window: {analysis_window_km:.0f} km; shift: {shift_km_sw:.0f} km (default)")
            mw_analyzer_sw = MovingWindowAnalysis(
                dx=dx, dy=dy, rho_load=RHO_LOAD, rho_m=RHO_MANTLE, rho_infill=RHO_INFILL, g=GRAVITY
            )
            try:
                result_region = mw_analyzer_sw.analyze(
                    topo_region, moho_region,
                    window_size=analysis_window_m,
                    shift_distance=shift_m_sw,
                    Te_range=Te_range_sw,
                )
            except Exception as e:
                print(f"  Error running Te map within region: {e}. Falling back to one Te for whole region.")
                result_region = None
            if result_region is not None:
                x_off = X_topo[j0, i0]
                y_off = Y_topo[j0, i0]
                x_abs_sw = result_region["x_centers"] + x_off
                y_abs_sw = result_region["y_centers"] + y_off
                Te_map_region = result_region["Te_map"]
                rms_map_region = result_region["rms_map"]
                te_mean_sw = np.nanmean(Te_map_region)
                rms_mean_sw = np.nanmean(rms_map_region)
                single_window_result = {
                    "Te_map": Te_map_region,
                    "rms_map": rms_map_region,
                    "x_centers": x_abs_sw,
                    "y_centers": y_abs_sw,
                    "Te_mean": te_mean_sw,
                    "rms_mean": rms_mean_sw,
                    "window_bounds": window_bounds_str,
                    "mode": "map",
                }
                print(f"\n Te map within region: mean Te = {te_mean_sw/1000:.2f} km, mean RMS = {rms_mean_sw/1000:.2f} km")
                single_region_has_map = True
            else:
                single_region_has_map = False
        else:
            single_region_has_map = False

        if not single_region_has_map:
            # ---------- One Te for whole region ----------
            topo_window = topo_anom[j0:j1, i0:i1].copy()
            moho_window = moho_undulation[j0:j1, i0:i1].copy()
            topo_window = topo_window - np.mean(topo_window)
            moho_window = moho_window - np.mean(moho_window)
            print(f"  Te search range: {Te_min_sw/1000:.0f}-{Te_max_sw/1000:.0f} km")
            print(f"  Window: {window_bounds_str} ({topo_window.shape[0]} x {topo_window.shape[1]} points)")
            print("  Running Te inversion...")

            try:
                single_result = inverter.invert_elastic_thickness(
                    topo_window, moho_window, Te_range=Te_range_sw, method="bounded"
                )
                te_best = float(single_result["Te_best"])
                rms_best = float(single_result["rms_best"])
                if not (np.isfinite(te_best) and np.isfinite(rms_best) and te_best > 0):
                    print("  Warning: Inversion returned invalid Te or RMS; result may show N/A.")
                    te_best = np.nan
                    rms_best = np.nan
            except Exception as e:
                print(f"  Warning: Inversion failed ({e}). Result figures will show N/A.")
                te_best = np.nan
                rms_best = np.nan
            at_upper_bound = False
            if np.isfinite(te_best) and te_best > 0:
                tol = 0.01 * (Te_max_sw - Te_min_sw)
                if te_best >= Te_max_sw - tol:
                    at_upper_bound = True
                    print("  Note: Te is at the upper limit of the search range; true Te may be higher.")
                    print("        Consider increasing 'Maximum Te (km)' and re-running.")
            single_window_result = {
                "Te_best": te_best,
                "rms_best": rms_best,
                "window_bounds": window_bounds_str,
                "at_upper_bound": at_upper_bound,
                "mode": "single",
            }
            print(f"\n Single-window result:")
            if np.isfinite(te_best) and te_best > 0:
                print(f"  Te = {te_best / 1000:.2f} km")
                print(f"  RMS misfit = {rms_best/1000:.2f} km")
            else:
                print("  Te = N/A (inversion failed or undefined)")
                print("  RMS misfit = N/A")

        if ref_moho_km is None and ref_moho is not None:
            ref_moho_km = ref_moho / 1000.0
        elif ref_moho_km is None:
            ref_moho_km = stats["moho_mean"] / 1000.0
        filename_prefix = create_filename_prefix(RHO_LOAD, ref_moho_km, window_size_km_sw)
        print(f"\n Filename prefix: {filename_prefix}")

        # SAVE DATA (single window)
        print("\nSaving outputs.")
         
        save_dict = {
            "topography": topography,
            "topography_anomaly": topo_anom,
            "moho_depth": moho_depth,
            "moho_undulation": moho_undulation,
            "X": X_topo,
            "Y": Y_topo,
            "stats": stats,
            "single_window_bounds": window_bounds_str,
        }
        if single_window_result.get("mode") == "map":
            save_dict["single_Te_map"] = single_window_result["Te_map"]
            save_dict["single_rms_map"] = single_window_result["rms_map"]
            save_dict["single_x_centers"] = single_window_result["x_centers"]
            save_dict["single_y_centers"] = single_window_result["y_centers"]
            save_dict["single_Te_mean"] = single_window_result["Te_mean"]
            save_dict["single_rms_mean"] = single_window_result["rms_mean"]
        else:
            save_dict["single_Te"] = single_window_result["Te_best"]
            save_dict["single_rms"] = single_window_result["rms_best"]
        if perform_global == "y":
            save_dict["global_Te"] = global_result["Te_best"]
            save_dict["global_rms"] = global_result["rms_best"]
        npz_filename = (
            f"{filename_prefix}inversion_results.npz"
            if filename_prefix
            else "inversion_results.npz"
        )
        np.savez(os.path.join(output_folder, npz_filename), **save_dict)

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
        if single_window_result.get("mode") == "map":
            write_surfer_grd(
                os.path.join(output_folder, f"{grd_prefix}Te_single_region.grd"),
                single_window_result["x_centers"], single_window_result["y_centers"],
                single_window_result["Te_map"],
            )
            write_surfer_grd(
                os.path.join(output_folder, f"{grd_prefix}rms_single_region.grd"),
                single_window_result["x_centers"], single_window_result["y_centers"],
                single_window_result["rms_map"],
            )
        else:
            Te_const = np.full_like(topo_anom, single_window_result["Te_best"], dtype=float)
            rms_const = np.full_like(topo_anom, single_window_result["rms_best"], dtype=float)
            write_surfer_grd(
                os.path.join(output_folder, f"{grd_prefix}Te_single.grd"),
                X_topo, Y_topo, Te_const,
            )
            write_surfer_grd(
                os.path.join(output_folder, f"{grd_prefix}rms_single.grd"),
                X_topo, Y_topo, rms_const,
            )
            # Also create simple Te/RMS map images for the single-window result (full grid, constant values)
            extent_km_full = [xmin_km, xmax_km, ymin_km, ymax_km]
            from matplotlib.ticker import FuncFormatter as _FmtKm
            def _fmt_km(v, t): return f"{int(v)}"

            # Te map (km)
            fig_te_full, ax_te_full = plt.subplots(1, 1, figsize=(10, 8))
            im_te_full = ax_te_full.imshow(
                Te_const / 1000.0,
                extent=extent_km_full,
                cmap="jet",
                origin="lower",
                aspect="equal",
            )
            ax_te_full.set_title("Elastic thickness (single-window result)")
            ax_te_full.set_xlabel("X (km)")
            ax_te_full.set_ylabel("Y (km)")
            ax_te_full.xaxis.set_major_formatter(_FmtKm(_fmt_km))
            ax_te_full.yaxis.set_major_formatter(_FmtKm(_fmt_km))
            plt.colorbar(im_te_full, ax=ax_te_full, label="Te (km)")
            te_full_path = os.path.join(output_folder, f"{grd_prefix}Te_single_map.png")
            fig_te_full.savefig(te_full_path, dpi=300, bbox_inches="tight")
            plt.close(fig_te_full)
            print(f" Te map (single window): {te_full_path}")

            # RMS map (km)
            fig_rms_full, ax_rms_full = plt.subplots(1, 1, figsize=(10, 8))
            rms_plot_full = np.abs(rms_const) / 1000.0
            im_rms_full = ax_rms_full.imshow(
                rms_plot_full,
                extent=extent_km_full,
                cmap="RdYlGn_r",
                origin="lower",
                aspect="equal",
            )
            ax_rms_full.set_title("RMS misfit (single-window result)")
            ax_rms_full.set_xlabel("X (km)")
            ax_rms_full.set_ylabel("Y (km)")
            ax_rms_full.xaxis.set_major_formatter(_FmtKm(_fmt_km))
            ax_rms_full.yaxis.set_major_formatter(_FmtKm(_fmt_km))
            plt.colorbar(im_rms_full, ax=ax_rms_full, label="RMS misfit (km)")
            rms_full_path = os.path.join(output_folder, f"{grd_prefix}rms_single_map.png")
            fig_rms_full.savefig(rms_full_path, dpi=300, bbox_inches="tight")
            plt.close(fig_rms_full)
            print(f" RMS map (single window): {rms_full_path}")
        print(" Surfer .grd files saved.")

        # Single-window FIGURES (PNG)
        from matplotlib.patches import Rectangle as RectPatch
        extent_km_sw = [xmin_km, xmax_km, ymin_km, ymax_km]
        half_km_sw = window_size_km_sw / 2.0
        topo_km_sw = topography / 1000.0
        moho_km_sw = moho_depth / 1000.0
        tmin, tmax = np.nanmin(topo_km_sw), np.nanmax(topo_km_sw)
        if tmax - tmin < 0.001:
            tmin, tmax = tmin - 0.0005, tmax + 0.0005
        mmin, mmax = np.nanmin(moho_km_sw), np.nanmax(moho_km_sw)
        if mmax - mmin < 0.001:
            mmin, mmax = mmin - 0.0005, mmax + 0.0005

        fig_loc, ax_loc = plt.subplots(1, 2, figsize=(14, 6))
        ax_loc[0].imshow(topo_km_sw, extent=extent_km_sw, cmap="jet", origin="lower", aspect="equal", vmin=tmin, vmax=tmax)
        ax_loc[0].set_title("Topography (km)")
        ax_loc[0].set_xlabel("X (km)")
        ax_loc[0].set_ylabel("Y (km)")
        rect0 = RectPatch((x_center_km - half_km_sw, y_center_km - half_km_sw), window_size_km_sw, window_size_km_sw, linewidth=2, edgecolor="lime", facecolor="none")
        ax_loc[0].add_patch(rect0)
        ax_loc[1].imshow(moho_km_sw, extent=extent_km_sw, cmap="jet", origin="lower", aspect="equal", vmin=mmin, vmax=mmax)
        ax_loc[1].set_title("Moho depth (km)")
        ax_loc[1].set_xlabel("X (km)")
        ax_loc[1].set_ylabel("Y (km)")
        rect1 = RectPatch((x_center_km - half_km_sw, y_center_km - half_km_sw), window_size_km_sw, window_size_km_sw, linewidth=2, edgecolor="lime", facecolor="none")
        ax_loc[1].add_patch(rect1)
        fig_loc.suptitle(f"Single-window location (center {x_center_km:.0f}, {y_center_km:.0f} km; size {window_size_km_sw:.0f} km)")
        plt.tight_layout()
        loc_path = os.path.join(output_folder, f"{grd_prefix}single_window_location.png")
        fig_loc.savefig(loc_path, dpi=300, bbox_inches="tight")
        plt.close(fig_loc)
        print(f"  Saved figure: {loc_path}")

        if single_window_result.get("mode") == "map":
            # Te map and RMS map figures for region
            from matplotlib.ticker import FuncFormatter
            def fmt_km(v, t): return f"{int(v)}"
            xc, yc = single_window_result["x_centers"], single_window_result["y_centers"]
            Te_map_r = single_window_result["Te_map"]
            rms_map_r = single_window_result["rms_map"]
            ext = [xc.min()/1000, xc.max()/1000, yc.min()/1000, yc.max()/1000]
            Te_plot = Te_map_r / 1000.0
            finite_te = Te_plot[np.isfinite(Te_plot)]
            if len(finite_te) > 0:
                Te_plot = np.where(np.isfinite(Te_plot), Te_plot, np.nanmin(finite_te))
            fig_te_sw, ax_te_sw = plt.subplots(1, 1, figsize=(10, 8))
            im_te = ax_te_sw.imshow(Te_plot, extent=ext, cmap="jet", origin="lower", aspect="equal")
            ax_te_sw.set_title(f"Te within selected region (mean = {single_window_result['Te_mean']/1000:.2f} km)")
            ax_te_sw.set_xlabel("X (km)")
            ax_te_sw.set_ylabel("Y (km)")
            ax_te_sw.xaxis.set_major_formatter(FuncFormatter(fmt_km))
            ax_te_sw.yaxis.set_major_formatter(FuncFormatter(fmt_km))
            plt.colorbar(im_te, ax=ax_te_sw, label="Te (km)")
            plt.tight_layout()
            fig_te_sw.savefig(os.path.join(output_folder, f"{grd_prefix}single_region_Te_map.png"), dpi=300, bbox_inches="tight")
            plt.close(fig_te_sw)
            print(f"  Saved figure: {grd_prefix}single_region_Te_map.png")
            fig_rms_sw, ax_rms_sw = plt.subplots(1, 1, figsize=(10, 8))
            rms_plot = np.abs(rms_map_r) / 1000.0
            rms_plot = np.where(np.isfinite(rms_plot), rms_plot, np.nan)
            im_rms = ax_rms_sw.imshow(rms_plot, extent=ext, cmap="jet", origin="lower", aspect="equal")
            ax_rms_sw.set_title(f"RMS misfit within region (mean = {single_window_result['rms_mean']/1000:.2f} km)")
            ax_rms_sw.set_xlabel("X (km)")
            ax_rms_sw.set_ylabel("Y (km)")
            ax_rms_sw.xaxis.set_major_formatter(FuncFormatter(fmt_km))
            ax_rms_sw.yaxis.set_major_formatter(FuncFormatter(fmt_km))
            plt.colorbar(im_rms, ax=ax_rms_sw, label="RMS (km)")
            plt.tight_layout()
            fig_rms_sw.savefig(os.path.join(output_folder, f"{grd_prefix}single_region_rms_map.png"), dpi=300, bbox_inches="tight")
            plt.close(fig_rms_sw)
            print(f"  Saved figure: {grd_prefix}single_region_rms_map.png")
        else:
            fig_res, ax_res = plt.subplots(1, 1, figsize=(6, 3.5))
            ax_res.set_xlim(0, 1)
            ax_res.set_ylim(0, 1)
            ax_res.axis("off")
            te_val = single_window_result["Te_best"]
            rms_m = single_window_result["rms_best"]
            at_upper = single_window_result.get("at_upper_bound", False)
            trans = ax_res.transAxes
            ax_res.text(0.5, 0.9, "Single-window result", fontsize=14, fontweight="bold", ha="center", transform=trans)
            if np.isfinite(te_val) and te_val > 0:
                te_km = te_val / 1000.0
                ax_res.text(0.5, 0.6, f"Te = {te_km:.2f} km", fontsize=12, ha="center", transform=trans)
                if at_upper:
                    ax_res.text(0.5, 0.48, "(at upper limit of search range; true Te may be higher)", fontsize=9, ha="center", style="italic", color="gray", transform=trans)
            else:
                ax_res.text(0.5, 0.6, "Te = N/A (inversion failed or undefined)", fontsize=12, ha="center", transform=trans)
            if np.isfinite(rms_m):
                ax_res.text(0.5, 0.32, f"RMS misfit = {rms_m/1000:.2f} km", fontsize=12, ha="center", transform=trans)
            else:
                ax_res.text(0.5, 0.32, "RMS misfit = N/A", fontsize=12, ha="center", transform=trans)
            ax_res.text(0.5, 0.12, f"Window: {window_bounds_str}", fontsize=10, ha="center", style="italic", transform=trans)
            fig_res.suptitle("Elastic thickness (Te) inversion")
            plt.tight_layout()
            res_path = os.path.join(output_folder, f"{grd_prefix}single_window_result.png")
            fig_res.savefig(res_path, dpi=300, bbox_inches="tight")
            plt.close(fig_res)
            print(f"  Saved figure: {res_path}")

    elif analysis_mode == "moving":
        # ---------- MOVING WINDOW ANALYSIS ----------
        if use_gravity == "y":
            print("\n Performing moving window analysis using gravity-predicted Moho")
            print("  (Will estimate Te by comparing gravity-predicted vs flexure-predicted Moho)")
        else:
            print("\n Performing moving window analysis using observed Moho")
            print("  (Will estimate Te by comparing observed vs flexure-predicted Moho)")
        print("\nMoving-window analysis.")
         

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

        print("\n Analysis Parameters:")
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
        print(f"\n Updated filename prefix with window size: {filename_prefix}")

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
        print("\nCreating Te maps.")
         

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

            # Calculate extent in km
            x_centers_m = result["x_centers"] + X_topo[0, 0]
            y_centers_m = result["y_centers"] + Y_topo[0, 0]
            extent_km = [
                x_centers_m.min() / 1000,
                x_centers_m.max() / 1000,
                y_centers_m.min() / 1000,
                y_centers_m.max() / 1000,
            ]

            # Prepare Te data - fill NaN with minimum (already in km)
            Te_data = result["Te_map"] / 1000
            Te_data_filled = Te_data.copy()
            Te_data_filled[np.isnan(Te_data_filled)] = Te_min / 1000

            # Plot
            im_te = ax_te.imshow(
                Te_data_filled,
                extent=extent_km,
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
            ax_te.set_xlabel("X (km)", fontsize=12)
            ax_te.set_ylabel("Y (km)", fontsize=12)
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

            # 3D Te map (coordinates in km)
            fig_te_3d = plt.figure(figsize=(12, 10))
            ax_te_3d = fig_te_3d.add_subplot(111, projection="3d")

            x_centers_km = x_centers_m / 1000.0
            y_centers_km = y_centers_m / 1000.0
            Xg, Yg = np.meshgrid(x_centers_km, y_centers_km)
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
            ax_te_3d.set_xlabel("X (km)", fontsize=12)
            ax_te_3d.set_ylabel("Y (km)", fontsize=12)
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

            # Calculate extent in km
            x_centers_m = result["x_centers"] + X_topo[0, 0]
            y_centers_m = result["y_centers"] + Y_topo[0, 0]
            extent_km = [
                x_centers_m.min() / 1000,
                x_centers_m.max() / 1000,
                y_centers_m.min() / 1000,
                y_centers_m.max() / 1000,
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

            # Plot (NaN values remain NaN, so white = no RMS estimate)
            im_rms = ax_rms.imshow(
                rms_data,
                extent=extent_km,
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
            ax_rms.set_xlabel("X (km)", fontsize=12)
            ax_rms.set_ylabel("Y (km)", fontsize=12)
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

            # 3D RMS map (coordinates in km)
            from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  - needed for 3D plots
            fig_rms_3d = plt.figure(figsize=(12, 10))
            ax_rms_3d = fig_rms_3d.add_subplot(111, projection="3d")

            x_centers_km = x_centers_m / 1000.0
            y_centers_km = y_centers_m / 1000.0
            Xg, Yg = np.meshgrid(x_centers_km, y_centers_km)
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
            ax_rms_3d.set_xlabel("X (km)", fontsize=12)
            ax_rms_3d.set_ylabel("Y (km)", fontsize=12)
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
        print("\nCreating RMS maps.")
         
        
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
        print("\nSaving outputs.")
         

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
    print("\nPlotting final results.")
     
    
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

    print(f"\nAnalysis complete. All outputs saved in: {output_folder}")

    # Close terminal logging and restore stdout
    if 'tee' in locals():
        tee.close()
        print(f"\n Terminal output saved to: {log_file_path}")

    plt.show()


if __name__ == "__main__":
    main()
