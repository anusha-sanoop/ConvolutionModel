"""
Main program for inverse modeling of elastic thickness
Implements the Braitenberg et al. (2002) convolution method

Usage:
    python main.py [input_file1] [input_file2] [output_dir]
    
     Or run without arguments for interactive mode.
"""

import sys
import os
from pathlib import Path
from datetime import datetime
from grd_reader import GRDReader, GridData

# Try to import netCDF4 for binary GRD support
try:
    import netCDF4
    HAS_NETCDF = True
except ImportError:
    HAS_NETCDF = False

from inverse_model import InverseModel, InverseParams
import numpy as np

# Import visualization functions
try:
    from visualization import (
        save_te_map,
        save_deflection_map,
        save_misfit_map,
        save_comparison_plot,
        save_te_vs_misfit_plot,
        save_shift_comparison_grid
    )
    HAS_VISUALIZATION = True
except ImportError:
    HAS_VISUALIZATION = False
    print("Warning: matplotlib not available, visualization disabled")


def get_user_input(prompt: str, default_value: str = "") -> str:
    """Get input from user interactively"""
    if default_value:
        full_prompt = f"{prompt} [{default_value}]: "
    else:
        full_prompt = f"{prompt}: "
    
    user_input = input(full_prompt).strip()
    
    # Remove quotes if present
    if user_input.startswith('"') and user_input.endswith('"'):
        user_input = user_input[1:-1]
    if user_input.startswith("'") and user_input.endswith("'"):
        user_input = user_input[1:-1]
    
    # Use default if empty
    if not user_input and default_value:
        return default_value
    
    return user_input


def display_accepted_input(input_file1: str, input_file2: str, output_dir: str,
                          grid1: GridData, grid2: GridData):
    """Display accepted input information in a clear format"""
    print("\n" + "=" * 70)
    print("                     ACCEPTED INPUT")
    print("=" * 70)
    
    print("\n[INPUT FILE 1]")
    print(f"  Path: {input_file1}")
    print("  Type: Input Data (Topography/Bouguer Gravity)")
    exists1 = GRDReader.file_exists(input_file1)
    print(f"  Status: {'✓ Found' if exists1 else '✗ Not Found'}")
    if exists1 and grid1.data is not None:
        print(f"  Dimensions: {grid1.nx} columns × {grid1.ny} rows")
        print("  Spatial extent:")
        print(f"    X: {grid1.xmin:.6f} to {grid1.xmax:.6f} (spacing: {grid1.dx:.6f})")
        print(f"    Y: {grid1.ymin:.6f} to {grid1.ymax:.6f} (spacing: {grid1.dy:.6f})")
        print(f"  Data range: {grid1.z_min:.3f} to {grid1.z_max:.3f}")
    
    print("\n[INPUT FILE 2]")
    print(f"  Path: {input_file2}")
    print("  Type: Observed Moho Depth/Flexure")
    exists2 = GRDReader.file_exists(input_file2)
    print(f"  Status: {'✓ Found' if exists2 else '✗ Not Found'}")
    if exists2 and grid2.data is not None:
        print(f"  Dimensions: {grid2.nx} columns × {grid2.ny} rows")
        print("  Spatial extent:")
        print(f"    X: {grid2.xmin:.6f} to {grid2.xmax:.6f} (spacing: {grid2.dx:.6f})")
        print(f"    Y: {grid2.ymin:.6f} to {grid2.ymax:.6f} (spacing: {grid2.dy:.6f})")
        print(f"  Data range: {grid2.z_min:.3f} to {grid2.z_max:.3f}")
    
    print("\n[OUTPUT DIRECTORY]")
    print(f"  Path: {output_dir}")
    
    print("\n" + "=" * 70)
    print()


def main():
    """Main program entry point"""
    print("=" * 40)
    print("Inverse Modeling of Elastic Thickness")
    print("Convolution Method Implementation")
    print("Based on Braitenberg et al. (2002)")
    print("=" * 40)
    print()
    
    input_file1 = ""
    input_file2 = ""
    output_dir = "output/"
    grid1 = GridData()
    grid2 = GridData()
    
    # Parse command line arguments or get from user
    if len(sys.argv) >= 3:
        # Command line mode
        input_file1 = sys.argv[1]
        input_file2 = sys.argv[2]
        # Default output directory
        output_dir = sys.argv[3] if len(sys.argv) >= 4 else r"D:\Marine\Project\Synthetic Data\Convolution\Implementation\Outputs" + os.sep
        
        print("Command line mode detected.")
    else:
        # Interactive mode
        print("Interactive input mode.")
        print("Please provide the following information:\n")
        
        # Get first input file
        input_file1 = get_user_input("Enter path to Input File 1 (Topography/Bouguer Gravity GRD file)")
        if not input_file1:
            print("\n✗ ERROR: Input File 1 path cannot be empty!")
            return 1
        
        # Get second input file
        input_file2 = get_user_input("Enter path to Input File 2 (Observed Moho Depth/Flexure GRD file)")
        if not input_file2:
            print("\n✗ ERROR: Input File 2 path cannot be empty!")
            return 1
        
        # Use default output directory
        output_dir = r"D:\Marine\Project\Synthetic Data\Convolution\Implementation\Outputs" + os.sep
    
    # Check if files exist and read them
    print("\nValidating input files...")
    
    if not GRDReader.file_exists(input_file1):
        print(f"\n✗ ERROR: Input File 1 not found: {input_file1}")
        return 1
    
    if not GRDReader.file_exists(input_file2):
        print(f"\n✗ ERROR: Input File 2 not found: {input_file2}")
        return 1
    
    # Read input data
    print("Reading input data files...")
    try:
        grid1 = GRDReader.read_grd(input_file1)
        print(f"  File 1 read: {grid1.nx} x {grid1.ny}, data shape: {grid1.data.shape if grid1.data is not None else 'None'}")
    except Exception as e:
        print(f"\n✗ ERROR: Failed to read input file 1: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1
    
    try:
        grid2 = GRDReader.read_grd(input_file2)
        print(f"  File 2 read: {grid2.nx} x {grid2.ny}, data shape: {grid2.data.shape if grid2.data is not None else 'None'}")
    except Exception as e:
        print(f"\n✗ ERROR: Failed to read input file 2: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Display accepted input
    display_accepted_input(input_file1, input_file2, output_dir, grid1, grid2)
    
    # Check grid compatibility
    if grid1.nx != grid2.nx or grid1.ny != grid2.ny:
        print("\n✗ ERROR: Grid dimensions don't match!")
        print(f"  Input File 1: {grid1.nx} × {grid1.ny}")
        print(f"  Input File 2: {grid2.nx} × {grid2.ny}")
        print("\nBoth grids must have the same dimensions.")
        return 1
    
    # Assign to meaningful variable names
    topography = grid1
    observed_moho = grid2
    
    # Set up inverse modeling parameters
    # Use Mars parameters (can be changed to get_default_eastern_alps_params() for Earth)
    use_sliding_window = True  # Set to False for global Te estimation
    run_multiple_shifts = True  # Set to True to run analysis for multiple shift values
    
    # Define shift values to test (in km, will be converted to meters)
    shift_values_km = [20.0, 30.0, 40.0]  # Different shift values to compare
    
    if use_sliding_window:
        # Base parameters
        base_params = InverseModel.get_mars_params(
            window_size=1000000.0,    # 1000 km in meters
            shift_min=20000.0,        # 20 km in meters (initial value)
            shift_max=80000.0,        # 80 km in meters (initial value)
            shift_step=20000.0,       # 20 km in meters
            Te_min=5.0,
            Te_max=50.0,
            Te_step=1.0
        )
    else:
        base_params = InverseModel.get_mars_params(
            window_size=None,
            shift_min=None,
            shift_max=None,
            shift_step=None,
            Te_min=5.0,
            Te_max=50.0,
            Te_step=1.0
        )
        base_params.use_sliding_window = False
    
    print("Inverse modeling parameters:")
    print(f"  Te range: {base_params.Te_min} - {base_params.Te_max} km")
    print(f"  Te step: {base_params.Te_step} km")
    print(f"  Crustal density: {base_params.rho_c} kg/m³")
    print(f"  Mantle density: {base_params.rho_m} kg/m³")
    print(f"  Gravity: {base_params.g} m/s²")
    print(f"  Young's modulus: {base_params.E / 1e9:.1f} GPa")
    print(f"  Poisson's ratio: {base_params.nu}")
    if base_params.use_sliding_window:
        print(f"  Sliding window: Enabled")
        print(f"    Window size: {base_params.window_size / 1000.0:.1f} km")
        if run_multiple_shifts:
            print(f"    Multiple shifts: Will test {len(shift_values_km)} shift values: {shift_values_km} km")
        else:
            print(f"    Shift range: {base_params.shift_min / 1000.0:.1f} - {base_params.shift_max / 1000.0:.1f} km")
            print(f"    Shift step: {base_params.shift_step / 1000.0:.1f} km")
    else:
        print(f"  Sliding window: Disabled (global Te estimation)")
    print()
    
    # Perform inverse modeling
    results_list = []
    te_maps_list = []
    rms_maps_list = []
    
    try:
        modeler = InverseModel()
        
        if base_params.use_sliding_window and run_multiple_shifts:
            # Run analysis for multiple shift values
            print(f"Running analysis for {len(shift_values_km)} different shift values...")
            for shift_km in shift_values_km:
                print(f"\n{'='*70}")
                print(f"Processing shift = {shift_km:.0f} km")
                print(f"{'='*70}")
                
                # Create params for this shift value
                # shift_km is the spacing (step size) between window centers
                # We want to analyze windows spaced by this distance across the entire grid
                # Set shift_step to control spacing, leave min/max as None to cover full grid
                params = InverseModel.get_mars_params(
                    window_size=base_params.window_size,
                    shift_min=None,  # Will start from window edge
                    shift_max=None,  # Will go to grid edge
                    shift_step=shift_km * 1000.0,  # Step = spacing between window centers (in meters)
                    Te_min=base_params.Te_min,
                    Te_max=base_params.Te_max,
                    Te_step=base_params.Te_step
                )
                
                # Ensure use_sliding_window is enabled
                params.use_sliding_window = True
                
                # Run analysis
                result = modeler.estimate_te_spatial_variable(topography, observed_moho, params)
                results_list.append(result)
                te_maps_list.append(result.Te_map)
                
                # Calculate RMS map (misfit map already has squared differences)
                rms_map = GridData()
                rms_map.nx = result.misfit_map.nx
                rms_map.ny = result.misfit_map.ny
                rms_map.dx = result.misfit_map.dx
                rms_map.dy = result.misfit_map.dy
                rms_map.xmin = result.misfit_map.xmin
                rms_map.xmax = result.misfit_map.xmax
                rms_map.ymin = result.misfit_map.ymin
                rms_map.ymax = result.misfit_map.ymax
                rms_map.data = np.sqrt(result.misfit_map.data)  # Convert from squared to RMS
                rms_map.update_range()
                rms_maps_list.append(rms_map)
            
            # Use the middle shift result as the main result for file output
            main_result = results_list[len(results_list) // 2]
        else:
            # Single analysis
            if base_params.use_sliding_window:
                print("Using sliding window approach for spatially variable Te estimation...")
                main_result = modeler.estimate_te_spatial_variable(topography, observed_moho, base_params)
            else:
                print("Using global Te estimation...")
                main_result = modeler.estimate_te(topography, observed_moho, base_params)
            
            results_list = [main_result]
            shift_values_km = [base_params.shift_step / 1000.0 if base_params.shift_step else 20.0]
    
    except Exception as e:
        print(f"\n✗ ERROR during inverse modeling: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1
    
    # Create timestamped output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir_timestamped = os.path.join(output_dir, f"Output_{timestamp}")
    output_path = Path(output_dir_timestamped)
    
    try:
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"\nOutput directory: {output_dir_timestamped}")
    except Exception as e:
        print(f"\n✗ ERROR: Cannot create output directory {output_dir_timestamped}: {str(e)}")
        return 1
    
    # Save a copy of input parameters to the output folder
    try:
        params_file = output_path / "parameters.txt"
        with open(params_file, 'w') as f:
            f.write(f"Inverse Modeling Parameters - {timestamp}\n")
            f.write("="*50 + "\n\n")
            f.write(f"Input File 1: {input_file1}\n")
            f.write(f"Input File 2: {input_file2}\n")
            f.write(f"\nParameters:\n")
            f.write(f"  Te range: {base_params.Te_min} - {base_params.Te_max} km\n")
            f.write(f"  Te step: {base_params.Te_step} km\n")
            f.write(f"  Crustal density: {base_params.rho_c} kg/m³\n")
            f.write(f"  Mantle density: {base_params.rho_m} kg/m³\n")
            f.write(f"  Gravity: {base_params.g} m/s²\n")
            f.write(f"  Young's modulus: {base_params.E / 1e9:.1f} GPa\n")
            f.write(f"  Poisson's ratio: {base_params.nu}\n")
            if base_params.use_sliding_window:
                f.write(f"  Window size: {base_params.window_size / 1000.0:.1f} km\n")
                if run_multiple_shifts:
                    f.write(f"  Shift values: {shift_values_km} km\n")
    except Exception as e:
        print(f"Warning: Could not save parameters file: {e}")
    
    # Save results
    print("\nSaving results...")
    
    # Save main result files
    result = main_result
    te_output = str(output_path / "Te_estimated.grd")
    deflection_output = str(output_path / "deflection_modeled.grd")
    misfit_output = str(output_path / "misfit.grd")
    
    try:
        GRDReader.write_grd(te_output, result.Te_map)
        print(f"  Written: {te_output}")
        
        GRDReader.write_grd(deflection_output, result.modeled_deflection)
        print(f"  Written: {deflection_output}")
        
        GRDReader.write_grd(misfit_output, result.misfit_map)
        print(f"  Written: {misfit_output}")
        
        # Save results for each shift value if multiple shifts were used
        if len(results_list) > 1:
            for i, (res, shift_km) in enumerate(zip(results_list, shift_values_km)):
                shift_str = f"{shift_km:.0f}km"
                GRDReader.write_grd(str(output_path / f"Te_estimated_shift_{shift_str}.grd"), res.Te_map)
                GRDReader.write_grd(str(output_path / f"misfit_shift_{shift_str}.grd"), res.misfit_map)
    except Exception as e:
        print(f"\n✗ ERROR saving output files: {str(e)}")
        return 1
    
    # Save visualization images
    if HAS_VISUALIZATION:
        print("\nGenerating and saving visualization images...")
        try:
            # Save Te map
            te_image = str(output_path / "Te_estimated.png")
            save_te_map(result.Te_map, te_image, "Estimated Elastic Thickness (Te)")
            print(f"  Written: {te_image}")
            
            # Save deflection map
            deflection_image = str(output_path / "deflection_modeled.png")
            save_deflection_map(result.modeled_deflection, deflection_image, "Modeled Moho Deflection")
            print(f"  Written: {deflection_image}")
            
            # Save misfit map
            misfit_image = str(output_path / "misfit.png")
            save_misfit_map(result.misfit_map, misfit_image, "Misfit Map")
            print(f"  Written: {misfit_image}")
            
            # Save comparison plot
            comparison_image = str(output_path / "observed_vs_modeled.png")
            save_comparison_plot(observed_moho, result.modeled_deflection, comparison_image, 
                               "Observed vs Modeled Moho Depth Comparison")
            print(f"  Written: {comparison_image}")
            
            # Save Te vs Misfit plot (if available)
            if result.Te_values and result.misfits:
                te_misfit_image = str(output_path / "Te_vs_misfit.png")
                save_te_vs_misfit_plot(result.Te_values, result.misfits, te_misfit_image,
                                      "Elastic Thickness vs Misfit")
                print(f"  Written: {te_misfit_image}")
            
            # Save shift comparison grid if multiple shifts were analyzed
            if len(te_maps_list) > 1 and len(rms_maps_list) > 1:
                comparison_grid_image = str(output_path / "shift_comparison_grid.png")
                save_shift_comparison_grid(te_maps_list, rms_maps_list, shift_values_km,
                                         comparison_grid_image, "Te and RMS Maps Comparison")
                print(f"  Written: {comparison_grid_image}")
        
        except Exception as e:
            print(f"\n⚠ WARNING: Failed to generate some visualizations: {str(e)}")
            print("  Data files saved successfully, but some images may be missing.")
    else:
        print("\n⚠ Visualization disabled (matplotlib not available)")
    
    # Print summary
    print("\n" + "=" * 70)
    print("                         RESULTS SUMMARY")
    print("=" * 70)
    if len(results_list) > 1:
        print(f"Analysis completed for {len(results_list)} shift values:")
        for shift_km, res in zip(shift_values_km, results_list):
            print(f"  Shift {shift_km:.0f} km: Best Te = {res.best_Te:.2f} km, Misfit = {res.best_misfit:.3f}")
    else:
        print(f"Best Te: {result.best_Te:.2f} km")
        print(f"Best misfit (RMS): {result.best_misfit:.3f}")
        print(f"Best correlation: {result.best_correlation:.3f}")
    print("=" * 70)
    print("\nProcessing complete!")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

