#!/usr/bin/env python3
"""
Script to recursively search through the subfolders within OUTPUT/CHAI and OUTPUT/BOLTZ
to find all the JSON files containing pLDDT values, extract the complex_plddt values,
and create a heatmap visualization.

For CHAI: Searches for outs.json files
For BOLTZ: Searches for confidence_[name]_model_0.json files

Usage:
    python plot_plddt_heatmap.py [--chai-output CHAI_DIR] [--boltz-output BOLTZ_DIR]
                                [--output OUTPUT_PNG] [--quiet]

Options:
    --chai-output CHAI_DIR   CHAI output directory (default: from config)
    --boltz-output BOLTZ_DIR BOLTZ output directory (default: from config)
    --output OUTPUT_PNG      Output PNG file for the heatmap (default: plots/plddt_heatmap.png)
    --quiet                  Suppress detailed output
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
    parser = argparse.ArgumentParser(description='Generate a heatmap of pLDDT values.')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG file for the heatmap (default: plots/plddt_heatmap.png)')
    
    # Add common arguments (which includes --quiet)
    parser = config_loader.add_common_args(parser)
    
    return parser.parse_args()

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
        
        # Check if this is a with_MSA path
        is_msa = False
        subfolder = ""
        
        if len(parts) >= 2:
            results_dir = parts[1]  # e.g., boltz_results_KORDshort_Arodyn1-6_with_MSA
            if '_with_MSA' in results_dir:
                is_msa = True
            
            # Extract the base name from boltz_results_[base_name]
            if results_dir.startswith('boltz_results_'):
                # Remove boltz_results_ prefix
                subfolder = results_dir[len('boltz_results_'):]
                # Remove _with_MSA suffix if present
                if subfolder.endswith('_with_MSA'):
                    subfolder = subfolder[:-9]  # Remove _with_MSA
        
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

def organize_data(chai_files, boltz_files, chai_dir, boltz_dir, config):
    """Extract and organize complex_plddt values from both CHAI and BOLTZ files.
    
    Organizes data with ligands on the y-axis and methods on the x-axis,
    similar to the RMSD heatmap.
    """
    # Initialize a nested dictionary with ligands as the outer key and methods as the inner key
    data = defaultdict(dict)
    
    # Process CHAI files if enabled in config
    if config["methods"]["use_chai"]:
        for file_path in chai_files:
            is_msa, subfolder = get_chai_identifier_from_path(file_path, chai_dir)
            
            # Skip MSA files if MSA is disabled
            if is_msa and not config["methods"]["use_msa"]:
                continue
                
            plddt = extract_chai_plddt(file_path)
            
            if plddt is not None and subfolder:
                # Determine the method based on whether it's MSA or not
                method = 'chai_with_MSA' if is_msa else 'chai'
                # Store the pLDDT value
                data[subfolder][method] = plddt
    
    # Process BOLTZ files if enabled in config
    if config["methods"]["use_boltz"]:
        for file_path in boltz_files:
            is_msa, subfolder = get_boltz_identifier_from_path(file_path, boltz_dir)
            
            # Skip MSA files if MSA is disabled
            if is_msa and not config["methods"]["use_msa"]:
                continue
                
            plddt = extract_boltz_plddt(file_path)
            
            if plddt is not None and subfolder:
                # Determine the method based on whether it's MSA or not
                method = 'boltz_with_MSA' if is_msa else 'boltz'
                # Store the pLDDT value
                data[subfolder][method] = plddt
    
    return data

def create_heatmap(data_dict, output_file, config, quiet=False):
    """Create a heatmap visualization of the complex_plddt values."""
    # Convert the nested dictionary to a pandas DataFrame
    df = pd.DataFrame(data_dict).T.fillna(0)
    
    # Determine available methods based on configuration and data
    available_methods = []
    all_methods = ['chai', 'chai_with_MSA', 'boltz', 'boltz_with_MSA']
    
    # Filter methods based on configuration
    if not config["methods"]["use_chai"]:
        all_methods = [m for m in all_methods if not m.startswith('chai')]
    if not config["methods"]["use_boltz"]:
        all_methods = [m for m in all_methods if not m.startswith('boltz')]
    if not config["methods"]["use_msa"]:
        all_methods = [m for m in all_methods if not '_with_MSA' in m]
    
    # Filter methods based on available data
    for method in all_methods:
        if method in df.columns:
            available_methods.append(method)
    
    if not available_methods:
        print("No methods available for visualization.")
        return None
    
    # Reindex with available methods
    df = df.reindex(columns=available_methods)
    
    # Sort the ligands (rows) alphabetically
    df = df.sort_index()
    
    # Create a figure with appropriate size
    plt.figure(figsize=(10, max(8, len(df) * 0.4)))
    
    # Use the RdYlGn colormap (Red-Yellow-Green)
    # For pLDDT, higher values are better, so we use the standard orientation
    # (green for high values, red for low values)
    
    # Create the heatmap
    ax = sns.heatmap(df, annot=True, cmap="RdYlGn", fmt=".3f", linewidths=.5)
    
    # Set title and labels
    plt.title("Complex pLDDT Values by Method and Ligand", fontsize=14)
    plt.xlabel("Method", fontsize=12)
    plt.ylabel("Ligand", fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    
    # Adjust layout
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300)
    if not quiet:
        print(f"Heatmap saved to {output_file}")
    
    return df

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    config = config_loader.load_config()
    config = config_loader.update_config_from_args(config, args)
    
    # Create plots directory if it doesn't exist
    plots_dir = Path(config["directories"]["plots"])
    if not plots_dir.exists():
        plots_dir.mkdir()
        if not args.quiet:
            print(f"Created directory: {plots_dir}")
    
    # Create csv directory if it doesn't exist
    csv_dir = Path(config["directories"]["csv"])
    if not csv_dir.exists():
        csv_dir.mkdir()
        if not args.quiet:
            print(f"Created directory: {csv_dir}")
    
    # Define the output file
    if args.output:
        output_file = Path(args.output)
    else:
        output_file = plots_dir / 'plddt_heatmap.png'
    
    # Define the root directories
    chai_dir = Path(config["directories"]["chai_output"])
    boltz_dir = Path(config["directories"]["boltz_output"])
    
    # Find all JSON files with pLDDT values
    chai_files = []
    boltz_files = []
    
    if config["methods"]["use_chai"]:
        chai_files = find_chai_json_files(chai_dir, args.quiet)
    
    if config["methods"]["use_boltz"]:
        boltz_files = find_boltz_json_files(boltz_dir, args.quiet)
    
    total_files = len(chai_files) + len(boltz_files)
    if not args.quiet:
        print(f"Found a total of {total_files} JSON files with pLDDT values")
    
    if total_files == 0:
        print("No files found. Please check if the files exist and try again.")
        return
    
    # Extract and organize the data
    data = organize_data(chai_files, boltz_files, chai_dir, boltz_dir, config)
    
    if not data:
        print("No data found to visualize.")
        return
    
    # Create the heatmap
    df = create_heatmap(data, output_file, config, args.quiet)
    
    if df is not None:
        # Save the data to CSV for further analysis
        csv_file = csv_dir / 'plddt_values.csv'
        df.to_csv(csv_file)
        if not args.quiet:
            print(f"Data saved to {csv_file}")

if __name__ == "__main__":
    main()
