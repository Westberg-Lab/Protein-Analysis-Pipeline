# Directory Structure Changes

## Overview

This document outlines the changes made to the protein prediction pipeline to improve the organization of output files. The main goal was to simplify the directory structure by:

1. Saving CSV files in the PSE_FILES directory instead of a separate csv directory
2. Creating a matching subdirectory structure in the plots directory to mirror the PSE_FILES directory
3. Organizing all analysis outputs by analysis run

## Problem Statement

Previously, the pipeline had the following issues:

- CSV files were saved in a separate csv directory, while PSE files were saved in the PSE_FILES directory
- Plot images were saved in the plots directory without any subdirectory structure
- This made it difficult to associate related files from the same analysis run

## Changes Made

### 1. In `extract_motif_plddt.py`:

- Removed references to the csv directory
- Updated the code to save CSV files in the PSE_FILES directory
- Changed the CSV file path to save to the PSE_FILES directory:
  ```python
  csv_file = pse_dir / f'motif_plddt_{args.motif}.csv'
  ```

### 2. In `plot_plddt.py`:

- Updated the `find_plddt_csv_files` function to look in the PSE_FILES directory and its subdirectories
- Modified the output file path to save plots in a subdirectory of the plots directory based on the analysis run name
- Added code to extract the analysis run name from command-line arguments or input path
- Fixed a bug in the data emptiness check
- Added a note to the plot indicating that BOLTZ values represent complex pLDDT (whole protein)

### 3. In `plot_rmsd.py`:

- Updated the `find_rmsd_csv_files` function to look in the PSE_FILES directory and its subdirectories
- Added an `--analysis-run` command-line argument
- Modified the output file path to save plots in a subdirectory of the plots directory based on the analysis run name
- Added code to extract the analysis run name from command-line arguments or input path

## Key Code Changes

### Finding CSV Files in Subdirectories

The `find_plddt_csv_files` and `find_rmsd_csv_files` functions were updated to look for CSV files in the PSE_FILES directory and its subdirectories:

```python
# Look for motif_plddt_*.csv files in the PSE directory
csv_files = list(pse_dir.glob(f"motif_plddt_{motif_id}*.csv"))
if not csv_files:
    csv_files = list(pse_dir.glob(f"motif_plddt_{motif_id}.csv"))

# Also look in subdirectories
for subdir in pse_dir.iterdir():
    if subdir.is_dir():
        subdir_files = list(subdir.glob(f"motif_plddt_{motif_id}*.csv"))
        if subdir_files:
            csv_files.extend(subdir_files)
        else:
            csv_file = subdir / f"motif_plddt_{motif_id}.csv"
            if csv_file.exists():
                csv_files.append(csv_file)
```

### Creating Subdirectories in the Plots Directory

The code was updated to create subdirectories in the plots directory based on the analysis run name:

```python
# Get analysis run name
analysis_run_name = None
if args.analysis_run:
    analysis_run_name = args.analysis_run
elif args.input:
    # Try to extract analysis run name from input path
    input_path = Path(args.input)
    if input_path.is_dir():
        # The directory name might be the analysis run name
        analysis_run_name = input_path.name
    elif input_path.parent.name != "csv":
        # The parent directory might be the analysis run name
        analysis_run_name = input_path.parent.name

# Create analysis run subdirectory in plots directory if we have an analysis run
if analysis_run_name:
    plot_dir = plots_dir / analysis_run_name
    plot_dir.mkdir(exist_ok=True, parents=True)
else:
    plot_dir = plots_dir
```

### Saving Plot Files in Subdirectories

The code was updated to save plot files in the appropriate subdirectory:

```python
# Set output file path
if args.motif:
    output_file = plot_dir / f'motif_plddt_heatmap_{args.motif}.png'
elif analysis_run_name:
    output_file = plot_dir / f'plddt_heatmap.png'
else:
    output_file = plot_dir / 'plddt_heatmap.png'
```

## Benefits

1. **Simplified Directory Structure**: All analysis outputs (PSE files, CSV files, plot images) are now organized by analysis run
2. **Consistent Organization**: The plots directory now mirrors the structure of the PSE_FILES directory
3. **Reduced Complexity**: Eliminated the need for a separate csv directory
4. **Improved Discoverability**: Related files are now grouped together, making it easier to find and manage analysis outputs

## Example Usage

```bash
# Extract motif pLDDT values for a specific motif
python src/extract_motif_plddt.py --motif=binding_pocket_hM3D

# Plot pLDDT values for a specific motif and analysis run
python src/plot_plddt.py --motif=binding_pocket_hM3D --analysis-run=binding_pocket_hM3D_analysis

# Plot RMSD values for a specific motif and analysis run
python src/plot_rmsd.py --motif=binding_pocket_hM3D --analysis-run=binding_pocket_hM3D_analysis
```

## Future Improvements

1. **Automatic Analysis Run Detection**: Enhance the scripts to automatically detect the analysis run name from the motif ID or other context
2. **Configuration Option**: Add a configuration option to specify the default analysis run name
3. **Command-Line Help**: Update the command-line help to document the new directory structure and command-line arguments
