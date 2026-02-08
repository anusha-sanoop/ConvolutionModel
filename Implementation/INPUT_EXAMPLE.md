# Input File Example

## Your Specified Input Files

1. **Input File 1**: `MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd`
   - Type: Bouguer Gravity data
   - Location: `Test data_DRK\Test data_DRK\Large and Small mount_New\`
   - This file contains the gravity anomaly data

2. **Input File 2**: `TwoMountainsFar_2Mohoflexure.grd`
   - Type: Observed Moho Flexure data
   - Location: `Test data_DRK\Test data_DRK\Large and Small mount_New\`
   - This file contains the observed Moho depth/flexure data

## Usage Options

### Option 1: Interactive Mode

Run the program without arguments:
```bash
./inverse_model
```

Then enter the paths when prompted:
```
Enter path to Input File 1 (Topography/Bouguer Gravity GRD file): D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd

Enter path to Input File 2 (Observed Moho Depth/Flexure GRD file): D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\TwoMountainsFar_2Mohoflexure.grd

Enter output directory [output/]: output/
```

### Option 2: Command Line Mode

Run with command line arguments:
```bash
./inverse_model.exe "D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd" "D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\TwoMountainsFar_2Mohoflexure.grd" "output/"
```

### Option 3: Using Batch File (Windows)

Edit `run_example.bat` with your paths and run:
```bash
run_example.bat
```

## Expected Output Display

The program will display the accepted input in a formatted way:

```
======================================================================
                     ACCEPTED INPUT
======================================================================

[INPUT FILE 1]
  Path: D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\MarsModel3TwoMountainsFar_3Gra_Bouguer_Coord_m.grd
  Type: Input Data (Topography/Bouguer Gravity)
  Status: ✓ Found
  Dimensions: [nx] columns × [ny] rows
  Spatial extent:
    X: [xmin] to [xmax] (spacing: [dx])
    Y: [ymin] to [ymax] (spacing: [dy])
  Data range: [z_min] to [z_max]

[INPUT FILE 2]
  Path: D:\Marine\Project\Synthetic Data\Convolution\Test data_DRK\Test data_DRK\Large and Small mount_New\TwoMountainsFar_2Mohoflexure.grd
  Type: Observed Moho Depth/Flexure
  Status: ✓ Found
  Dimensions: [nx] columns × [ny] rows
  Spatial extent:
    X: [xmin] to [xmax] (spacing: [dx])
    Y: [ymin] to [ymax] (spacing: [dy])
  Data range: [z_min] to [z_max]

[OUTPUT DIRECTORY]
  Path: output/

======================================================================
```

## Notes

- You can copy-paste the full paths directly (quotes will be automatically removed if present)
- The program will validate that both files exist before proceeding
- Both grids must have the same dimensions (nx × ny)
- Paths with spaces are handled automatically

