#!/usr/bin/env python3
"""
Unified script to generate heatmap visualizations of pLDDT values.

This script can:
1. Read pLDDT values from CSV files (motif_plddt_*.csv)
2. Extract pLDDT values directly from JSON files in OUTPUT/CHAI and OUTPUT/BOLTZ
3. Create heatmaps based on this data

Usage:
    python plot_plddt.py [--input INPUT_CSV] [--output OUTPUT_PNG]
                        [--chai-output CHAI_DIR] [--boltz-output BOLTZ_DIR]
                        [--motif MOTIF_ID] [--quiet]

Options:
    --input INPUT_CSV      Input CSV file or directory (default: auto-detect)
    --output OUTPUT_PNG    Output PNG file for the heatmap (default: plots/plddt_heatmap.png)
    --chai-output CHAI_DIR CHAI output directory (default: from config)
    --boltz-output BOLTZ_DIR BOLTZ output directory (default: from config)
    --motif MOTIF_ID       Specific motif to plot (default: all motifs)
    --quiet                Suppress detailed output
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
import argparse
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate heatmaps of pLDDT values.')
    parser.add_argument('--input', type=str, default=None,
                        help='Input CSV file or directory (default: auto-detect)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG file for the heatmap (default: plots/plddt_heatmap.png)')
    parser.add_argument('--motif', type=str, default=None,
                        help='Specific motif to plot (default: all motifs)')
    parser.add_argument('--analysis-run', type=str, default=None,
                        help='Analysis run ID to use (default: first enabled analysis run)')
    
    # Add common arguments
    parser = config_loader.add_common_args(parser, exclude=['motif'])
    
    return parser.parse_args()

def find_plddt_csv_files(pse_dir, motif_id=None, quiet=False):
    """Find pLDDT CSV files in the specified directory.
    
    Args:
        pse_dir: PSE files directory to search in
        motif_id: Optional motif ID to filter by
        quiet: Whether to suppress output
        
    Returns:
        List of paths to CSV files
    """
    if not pse_dir.exists() or not pse_dir.is_dir():
        print(f"{pse_dir} directory not found.")
        return []
    
    csv_files = []
    
    # If motif_id is specified, look for motif-specific pLDDT files
    if motif_id:
        if not quiet:
            print(f"Looking for motif-specific pLDDT files for motif {motif_id}...")
        
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
    else:
        # Look for all pLDDT files in the PSE directory
        csv_files = list(pse_dir.glob("plddt_values_*.csv"))
        if not csv_files:
            # If no plddt_values_*.csv files found, look for motif_plddt_*.csv files
            csv_files = list(pse_dir.glob("motif_plddt_*.csv"))
        
        # Also look in subdirectories
        for subdir in pse_dir.iterdir():
            if subdir.is_dir():
                subdir_files = list(subdir.glob("plddt_values_*.csv"))
                if subdir_files:
                    csv_files.extend(subdir_files)
                else:
                    subdir_files = list(subdir.glob("motif_plddt_*.csv"))
                    if subdir_files:
                        csv_files.extend(subdir_files)
    
    if not quiet:
        print(f"Found {len(csv_files)} pLDDT CSV files")
    
    return csv_files

def find_chai_json_files(root_dir, quiet=False):
    """Find all outs.json files in the CHAI output directory."""
    json_files = []
    if not quiet:
        print(f"Searching for outs.json files in {root_dir}...")
    
    if not root_dir.exists():
        print(f"Directory {root_dir} does not exist.")
        return json_files
    
    for path in root_dir.rglob('outs.json'):
        json_files.append(path)
    
    if not quiet:
        print(f"Found {len(json_files)} outs.json files in CHAI output")
    return json_files

def find_boltz_json_files(root_dir, quiet=False):
    """Find all confidence_*_model_0.json files in the BOLTZ output directory."""
    json_files = []
    if not quiet:
        print(f"Searching for confidence_*_model_0.json files in {root_dir}...")
    
    if not root_dir.exists():
        print(f"Directory {root_dir} does not exist.")
        return json_files
    
    # Look for confidence_*_model_0.json files in the predictions subdirectories
    for predictions_dir in root_dir.glob('*/boltz_results_*/predictions'):
        for subdir in predictions_dir.iterdir():
            if subdir.is_dir():
                confidence_file = subdir / f"confidence_{subdir.name}_model_0.json"
                if confidence_file.exists():
                    json_files.append(confidence_file)
    
    if not quiet:
        print(f"Found {len(json_files)} confidence_*_model_0.json files in BOLTZ output")
    return json_files

def get_chai_identifier_from_path(file_path, chai_dir):
    """Extract identifier from CHAI file path."""
    # Get the path relative to OUTPUT/CHAI
    try:
        rel_path = file_path.relative_to(chai_dir)
        # Extract the parts of the path (main folder and subfolder)
        parts = list(rel_path.parts)
        
        # Check if this is a with_MSA path
        is_msa = False
        if len(parts) >= 1 and parts[0].endswith('_with_MSA'):
            is_msa = True
        
        # Get the subfolder (ligand name)
        subfolder = ""
        if len(parts) >= 2:
            subfolder = parts[1]
        
        return (is_msa, subfolder)
    except ValueError:
        # If the file is not relative to chai_dir
        return (False, "")

def get_boltz_identifier_from_path(file_path, boltz_dir):
    """Extract identifier from BOLTZ file path."""
    # Get the path relative to OUTPUT/BOLTZ
    try:
        rel_path = file_path.relative_to(boltz_dir)
        parts = list(rel_path.parts)
        
        # The structure is typically:
        # OUTPUT/BOLTZ/[folder_name]/boltz_results_[base_name]/predictions/[name]/confidence_[name]_model_0.json
        # or for MSA:
        # OUTPUT/BOLTZ/[folder_name]_with_MSA/boltz_results_[base_name]/predictions/[name]/confidence_[name]_model_0.json
        
        # Check if this is a with_MSA path
        is_msa = False
        subfolder = ""
        
        if len(parts) >= 1:
            folder_name = parts[0]  # e.g., KORDshort or KORDshort_with_MSA
            
            # Check if this is an MSA path
            if folder_name.endswith('_with_MSA'):
                is_msa = True
        
        if len(parts) >= 2:
            results_dir = parts[1]  # e.g., boltz_results_KORDshort_Arodyn1-6
            
            # Extract the base name from boltz_results_[base_name]
            if results_dir.startswith('boltz_results_'):
                # Remove boltz_results_ prefix
                subfolder = results_dir[len('boltz_results_'):]
        
        return (is_msa, subfolder)
    except ValueError:
        # If the file is not relative to boltz_dir
        return (False, "")

def extract_chai_plddt(file_path):
    """Extract complex_plddt from a CHAI outs.json file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            # Extract the first value from the complex_plddt array
            complex_plddt = data.get("cand_0", {}).get("complex_plddt", [0])[0]
            return complex_plddt
    except Exception as e:
        print(f"Error extracting complex_plddt from CHAI file {file_path}: {e}")
        return None

def extract_boltz_plddt(file_path):
    """Extract complex_plddt from a BOLTZ confidence_*_model_0.json file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            # Extract the complex_plddt value
            complex_plddt = data.get("complex_plddt", 0)
            return complex_plddt
    except Exception as e:
        print(f"Error extracting complex_plddt from BOLTZ file {file_path}: {e}")
        return None

def organize_data_from_json(chai_files, boltz_files, chai_dir, boltz_dir, prediction_runs, quiet=False):
    """Extract and organize complex_plddt values from both CHAI and BOLTZ files.
    
    Organizes data with ligands on the y-axis and methods on the x-axis,
    similar to the RMSD heatmap.
    """
    # Initialize a nested dictionary with ligands as the outer key and methods as the inner key
    data = defaultdict(dict)
    
    # Process each prediction run
    for run in prediction_runs:
        run_id = run.get("id", "")
        methods = run.get("methods", {})
        
        # Process CHAI files if enabled for this run
        if methods.get("use_chai", True):
            if not quiet:
                print(f"Processing CHAI files for prediction run: {run_id}")
            
            for file_path in chai_files:
                is_msa, subfolder = get_chai_identifier_from_path(file_path, chai_dir)
                
                # Skip files that don't match the run's MSA configuration
                if is_msa != methods.get("use_msa", False):
                    continue
                    
                plddt = extract_chai_plddt(file_path)
                
                if plddt is not None and subfolder:
                    # Use simple, consistent method names
                    method = 'chai_with_MSA' if is_msa else 'chai'
                    
                    # Store the pLDDT value
                    data[subfolder][method] = plddt
        
        # Process BOLTZ files if enabled for this run
        if methods.get("use_boltz", True):
            if not quiet:
                print(f"Processing BOLTZ files for prediction run: {run_id}")
            
            for file_path in boltz_files:
                is_msa, subfolder = get_boltz_identifier_from_path(file_path, boltz_dir)
                
                # Skip files that don't match the run's MSA configuration
                if is_msa != methods.get("use_msa", False):
                    continue
                    
                plddt = extract_boltz_plddt(file_path)
                
                if plddt is not None and subfolder:
                    # Use simple, consistent method names
                    method = 'boltz_with_MSA' if is_msa else 'boltz'
                    
                    # Store the pLDDT value
                    data[subfolder][method] = plddt
    
    return data

def read_plddt_values(csv_files, quiet=False):
    """Read pLDDT values from CSV files.
    
    Args:
        csv_files: List of CSV files to read
        quiet: Whether to suppress output
        
    Returns:
        DataFrame with pLDDT values
    """
    all_data = []
    
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            if not quiet:
                print(f"Read {len(df)} pLDDT values from {csv_file}")
            all_data.append(df)
        except Exception as e:
            print(f"Error reading CSV file {csv_file}: {e}")
    
    if not all_data:
        print("No data found in CSV files")
        return None
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df

def create_plddt_heatmap(data, output_file, config, full_config=None, motif_id=None, quiet=False):
    """Create a heatmap visualization of pLDDT values.
    
    Args:
        data: DataFrame or dictionary with pLDDT values
        output_file: Path to save the heatmap
        config: Configuration dictionary
        full_config: Full configuration dictionary (for motif definitions)
        motif_id: Optional motif ID to filter by
        quiet: Whether to suppress output
        
    Returns:
        True if successful, False otherwise
    """
    # Convert data to DataFrame if it's a dictionary
    if isinstance(data, dict):
        df = pd.DataFrame(data).T.fillna(0)
    else:
        df = data
    
    # Check if DataFrame is valid
    if df is None or len(df) == 0 or df.shape[1] == 0:
        print("No data to visualize.")
        return False
    
    # Filter by motif if specified
    if motif_id and 'motif' in df.columns:
        df = df[df['motif'] == motif_id]
        if len(df) == 0:
            print(f"No data found for motif '{motif_id}'")
            return False
    
    # If DataFrame has 'ligand', 'method', and 'plddt' columns, pivot it
    if all(col in df.columns for col in ['ligand', 'method', 'plddt']):
        # Use pivot_table instead of pivot to handle duplicate entries by taking the mean
        pivot_df = df.pivot_table(index='ligand', columns='method', values='plddt', aggfunc='mean')
    else:
        # Assume it's already in the right format
        pivot_df = df
    
    # Get motif description if available
    motif_desc = motif_id
    if motif_id and full_config:
        motif_def = config_loader.get_motif_definition(full_config, motif_id)
        if motif_def:
            motif_desc = motif_def.get("description", motif_id)
    
    # Determine available methods based on configuration and data
    available_methods = []
    all_methods = ['chai', 'chai_with_MSA', 'boltz', 'boltz_with_MSA']
    
    # Filter methods based on configuration
    if not config.get("methods", {}).get("use_chai", True):
        all_methods = [m for m in all_methods if not m.startswith('chai')]
    if not config.get("methods", {}).get("use_boltz", True):
        all_methods = [m for m in all_methods if not m.startswith('boltz')]
    # Always include MSA methods regardless of the use_msa flag
    
    # Filter methods based on available data
    for method in all_methods:
        if method in pivot_df.columns:
            available_methods.append(method)
    
    if not available_methods:
        print("No methods available for visualization.")
        return False
    
    # Reindex with available methods
    pivot_df = pivot_df.reindex(columns=available_methods)
    
    # Sort the ligands (rows) alphabetically
    pivot_df = pivot_df.sort_index()
    
    # Create the heatmap
    plt.figure(figsize=(10, max(8, len(pivot_df) * 0.4)))
    
    # Use the RdYlGn colormap (Red-Yellow-Green)
    # For pLDDT, higher values are better, so we use the standard orientation
    # (green for high values, red for low values)
    
    # Create the heatmap
    ax = sns.heatmap(pivot_df, annot=True, cmap="RdYlGn", fmt=".3f", linewidths=.5)
    
    # Set title and labels
    if motif_id:
        plt.title(f"Motif-Specific pLDDT Values: {motif_desc}", fontsize=14)
    else:
        plt.title("Complex pLDDT Values by Method and Ligand", fontsize=14)
    
    plt.xlabel("Method", fontsize=12)
    plt.ylabel("Ligand", fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    
    # Add a note about BOLTZ using complex pLDDT if motif-specific
    if motif_id:
        # Adjust bottom margin to make room for the note
        plt.subplots_adjust(bottom=0.15)
        plt.figtext(0.1, 0.01, "Note: BOLTZ = complex pLDDT (whole protein)", 
                   ha="center", fontsize=10, style="italic")
    else:
        # Use tight layout for non-motif plots
        plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    if not quiet:
        print(f"Heatmap saved to {output_file}")
    
    return True

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    full_config = config_loader.load_config()
    
    # Get a merged configuration (combines global settings with a prediction run)
    config = config_loader.get_merged_config(full_config)
    
    # Update with command-line arguments
    config = config_loader.update_config_from_args(config, args)
    
    # Get directories from config
    # Handle both old and new configuration structures
    if "directories" in config:
        # Old configuration structure
        plots_dir = Path(config["directories"]["plots"])
        chai_dir = Path(config["directories"]["chai_output"])
        boltz_dir = Path(config["directories"]["boltz_output"])
        pse_dir = Path(config["directories"]["pse_files"])
    else:
        # New configuration structure
        plots_dir = Path(config.get("plots", "plots"))
        chai_dir = Path(config.get("chai_output", "OUTPUT/CHAI"))
        boltz_dir = Path(config.get("boltz_output", "OUTPUT/BOLTZ"))
        pse_dir = Path(config.get("pse_files", "PSE_FILES"))
    
    # Create output directories
    plots_dir.mkdir(exist_ok=True)
    pse_dir.mkdir(exist_ok=True, parents=True)
    
    # Determine input files and data source
    data = None
    
    if args.input:
        # Use specified input file or directory
        input_path = Path(args.input)
        if input_path.is_file():
            # Single file
            csv_files = [input_path]
            data = read_plddt_values(csv_files, args.quiet)
        elif input_path.is_dir():
            # Directory - find CSV files
            if args.motif:
                # Look for motif-specific pLDDT files
                csv_files = list(input_path.glob(f"motif_plddt_{args.motif}*.csv"))
                if not csv_files:
                    csv_files = list(input_path.glob("motif_plddt_*.csv"))
            else:
                # Look for all pLDDT files
                csv_files = list(input_path.glob("plddt_values_*.csv"))
                if not csv_files:
                    csv_files = list(input_path.glob("motif_plddt_*.csv"))
            
            if csv_files:
                data = read_plddt_values(csv_files, args.quiet)
        else:
            print(f"Input path {input_path} does not exist.")
            return
    else:
        # Auto-detect input files
        if args.motif:
            # Look for motif-specific pLDDT files in the PSE directory
            csv_files = find_plddt_csv_files(pse_dir, args.motif, args.quiet)
            if csv_files:
                data = read_plddt_values(csv_files, args.quiet)
        else:
            # Look for pLDDT files in the PSE directory
            csv_files = find_plddt_csv_files(pse_dir, None, args.quiet)
            if csv_files:
                data = read_plddt_values(csv_files, args.quiet)
            else:
                # If no CSV files found, extract pLDDT values from JSON files
                if not args.quiet:
                    print("No pLDDT CSV files found. Extracting pLDDT values from JSON files...")
                
                # Get the analysis run to use
                analysis_runs = full_config.get("analysis_runs", [])
                analysis_run = None
                
                if args.analysis_run:
                    # Find the specified analysis run
                    for run in analysis_runs:
                        if run.get("id") == args.analysis_run:
                            analysis_run = run
                            break
                    
                    if not analysis_run:
                        print(f"Analysis run '{args.analysis_run}' not found.")
                        return
                else:
                    # Get the first enabled analysis run
                    for run in analysis_runs:
                        if run.get("enabled", True):
                            analysis_run = run
                            break
                
                if not analysis_run:
                    print("No enabled analysis runs found.")
                    return
                
                # Get the source predictions for this analysis run
                source_prediction_ids = analysis_run.get("source_predictions", [])
                if not source_prediction_ids:
                    print(f"No source predictions specified for analysis run '{analysis_run.get('id')}'.")
                    return
                
                # Get the prediction runs
                prediction_runs = []
                for pred_id in source_prediction_ids:
                    pred_run = config_loader.get_prediction_run_by_id(full_config, pred_id)
                    if pred_run:
                        prediction_runs.append(pred_run)
                
                if not prediction_runs:
                    print("No valid prediction runs found.")
                    return
                
                # Find all JSON files with pLDDT values
                chai_files = find_chai_json_files(chai_dir, args.quiet)
                boltz_files = find_boltz_json_files(boltz_dir, args.quiet)
                
                total_files = len(chai_files) + len(boltz_files)
                if not args.quiet:
                    print(f"Found a total of {total_files} JSON files with pLDDT values")
                
                if total_files == 0:
                    print("No files found. Please check if the files exist and try again.")
                    return
                
                # Extract and organize the data
                data = organize_data_from_json(chai_files, boltz_files, chai_dir, boltz_dir, prediction_runs, args.quiet)
    
    if data is None or (isinstance(data, pd.DataFrame) and data.empty) or (isinstance(data, dict) and not data):
        print("No data found to visualize.")
        return
    
    # Determine output file
    if args.output:
        output_file = Path(args.output)
    else:
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
        elif analysis_run:  # Use the analysis run determined earlier in the code
            analysis_run_name = analysis_run.get("id")
        elif args.motif:
            # If no analysis run is specified, try to use the motif ID as part of the analysis run name
            # This is a common naming convention: motif_id + "_analysis"
            analysis_run_name = f"{args.motif}_analysis"
        
        # Create analysis run subdirectory in plots directory if we have an analysis run
        if analysis_run_name:
            plot_dir = plots_dir / analysis_run_name
            plot_dir.mkdir(exist_ok=True, parents=True)
        else:
            plot_dir = plots_dir
        
        # Set output file path
        if args.motif:
            output_file = plot_dir / f'motif_plddt_heatmap_{args.motif}.png'
        elif analysis_run_name:
            output_file = plot_dir / f'plddt_heatmap.png'
        else:
            output_file = plot_dir / 'plddt_heatmap.png'
    
    # Create heatmap
    success = create_plddt_heatmap(data, output_file, config, full_config, args.motif, args.quiet)
    
    if success:
        if not args.quiet:
            print(f"Heatmap created successfully: {output_file}")
        
        # Save the data to CSV for further analysis if it's not already from a CSV
        if isinstance(data, dict):
            df = pd.DataFrame(data).T
            if args.motif:
                # Save CSV file in the PSE directory with the analysis run name as a subdirectory
                if analysis_run_name:
                    pse_run_dir = pse_dir / analysis_run_name
                    pse_run_dir.mkdir(exist_ok=True, parents=True)
                    csv_file = pse_run_dir / f'motif_plddt_{args.motif}.csv'
                else:
                    csv_file = pse_dir / f'motif_plddt_{args.motif}.csv'
            elif args.analysis_run:
                # Save CSV file in the PSE directory with the analysis run name as a subdirectory
                pse_run_dir = pse_dir / args.analysis_run
                pse_run_dir.mkdir(exist_ok=True, parents=True)
                csv_file = pse_run_dir / f'plddt_values.csv'
            else:
                csv_file = pse_dir / 'plddt_values.csv'
            
            df.to_csv(csv_file)
            if not args.quiet:
                print(f"Data saved to {csv_file}")
    else:
        print(f"Failed to create heatmap")

if __name__ == "__main__":
    main()
