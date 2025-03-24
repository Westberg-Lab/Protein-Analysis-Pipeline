#!/usr/bin/env python3
"""
Script to generate a heatmap visualization of RMSD values.

This script reads RMSD values from a CSV file and creates a heatmap with:
- x-axis: methods (chai, chai_with_MSA, boltz, boltz_with_MSA)
- y-axis: ligand names (each .pse file in PSE_FILES folder)
- Cell values: RMSD values

Usage:
    python plot_rmsd_heatmap.py [--input INPUT_CSV] [--output OUTPUT_PNG]
                               [--reference REFERENCE] [--quiet]

Options:
    --input INPUT_CSV    Input CSV file with RMSD values (default: auto-detect in PSE_FILES subdirectories)
    --output OUTPUT_PNG  Output PNG file for the heatmap (default: plots/rmsd_heatmap_[reference].png)
    --reference REF      Reference name to filter by (default: use all references)
    --quiet              Suppress detailed output
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate a heatmap of RMSD values.')
    parser.add_argument('--input', type=str, default=None,
                        help='Input CSV file with RMSD values (default: auto-detect in PSE_FILES subdirectories)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG file for the heatmap (default: plots/rmsd_heatmap_[reference].png)')
    parser.add_argument('--reference', type=str, default=None,
                        help='Reference name to filter by (default: use all references)')
    
    # Add common arguments (which includes --quiet)
    parser = config_loader.add_common_args(parser)
    
    return parser.parse_args()

def get_templates(config):
    """Get all template files from configuration."""
    # Handle both old and new configuration structures
    if "templates" in config:
        # Old configuration structure
        if "files" in config["templates"]:
            return [Path(template) for template in config["templates"]["files"]]
        elif "default_template" in config["templates"]:
            return [Path(config["templates"]["default_template"])]
    else:
        # New configuration structure
        if "files" in config:
            return [Path(template) for template in config["files"]]
        elif "default_template" in config:
            return [Path(config["default_template"])]
    return []

def find_rmsd_csv_files(pse_dir, template_name=None):
    """Find all rmsd_values.csv files in PSE_FILES subdirectories.
    
    If template_name is provided, only look for CSV files in that template's subdirectory.
    """
    if not pse_dir.exists() or not pse_dir.is_dir():
        print(f"{pse_dir} directory not found.")
        return []
    
    csv_files = []
    
    if template_name:
        # Look for CSV files in the specific template subdirectory
        template_dir = pse_dir / template_name
        if template_dir.exists() and template_dir.is_dir():
            csv_file = template_dir / 'rmsd_values.csv'
            if csv_file.exists():
                csv_files.append(csv_file)
    else:
        # Look for CSV files in all subdirectories
        for subdir in pse_dir.iterdir():
            if subdir.is_dir():
                csv_file = subdir / 'rmsd_values.csv'
                if csv_file.exists():
                    csv_files.append(csv_file)
    
    return csv_files

def read_rmsd_values(csv_file, quiet=False):
    """Read RMSD values from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        if not quiet:
            print(f"Read {len(df)} RMSD values from {csv_file}")
        return df
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None

def create_heatmap(df, output_file, config, quiet=False):
    """Create a heatmap visualization of RMSD values."""
    # Check if DataFrame is valid
    if df is None or len(df) == 0:
        print("No data to visualize.")
        return False
    
    # Check if required columns exist
    required_columns = ['ligand', 'method', 'rmsd']
    if not all(col in df.columns for col in required_columns):
        print(f"CSV file must contain columns: {', '.join(required_columns)}")
        return False
    
    # Get reference name from the data if available
    reference_name = "Unknown"
    if 'reference' in df.columns and not df['reference'].empty:
        reference_name = df['reference'].iloc[0]
    
    # Pivot the data for the heatmap
    # This transforms the data from long format to wide format
    # where rows are ligands, columns are methods, and values are RMSD values
    try:
        pivot_df = df.pivot(index='ligand', columns='method', values='rmsd')
        
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
        
        # Use a natural color spectrum: 'RdYlGn_r' (Red-Yellow-Green reversed)
        # This is intuitive for RMSD values: green for good alignments (low RMSD),
        # yellow for medium, red for poor alignments (high RMSD)
        
        # Apply normalization to better visualize small differences
        from matplotlib.colors import Normalize
        
        # Get visualization parameters
        if "visualization" in config:
            vmin = config["visualization"]["rmsd_vmin"]
            vmax = config["visualization"]["rmsd_vmax"]
        else:
            vmin = config.get("rmsd_vmin", 0.2)
            vmax = config.get("rmsd_vmax", 6.2)
        
        # Create the heatmap with seaborn using the natural colormap and normalization
        ax = sns.heatmap(pivot_df, annot=True, cmap='RdYlGn_r', fmt='.4f', 
                     norm=Normalize(vmin=vmin, vmax=vmax),
                     cbar_kws={'label': 'RMSD (Å)'})
        
        # Set labels and title with reference protein name from the data
        plt.title(f'RMSD Values by Method and Ligand (Reference: {reference_name})', fontsize=14)
        plt.xlabel('Method', fontsize=12)
        plt.ylabel('Ligand', fontsize=12)
        
        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        
        # Adjust layout to make room for labels
        plt.tight_layout()
        
        # Save the figure
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        if not quiet:
            print(f"Heatmap saved to {output_file}")
        return True
    
    except Exception as e:
        print(f"Error creating heatmap: {e}")
        return False

def process_template(template_name, config, args, plots_dir, csv_dir):
    """Process RMSD data for a specific template."""
    if not args.quiet:
        print(f"\nProcessing template: {template_name}")
    
    # Determine input file(s)
    if args.input:
        # Use specified input file
        input_files = [Path(args.input)]
        if not input_files[0].exists():
            print(f"Input file {input_files[0]} does not exist.")
            return False
    else:
        # Auto-detect input files for this template
        # Handle both old and new configuration structures
        if "directories" in config:
            # Old configuration structure
            pse_dir = Path(config["directories"]["pse_files"])
        else:
            # New configuration structure
            pse_dir = Path(config.get("pse_files", "PSE_FILES"))
        
        input_files = find_rmsd_csv_files(pse_dir, template_name)
        if not input_files:
            print(f"No RMSD CSV files found for template {template_name} in {pse_dir} subdirectories.")
            return False
        if not args.quiet:
            print(f"Found {len(input_files)} RMSD CSV files for template {template_name}: {', '.join(str(f) for f in input_files)}")
    
    # Process each input file
    for input_file in input_files:
        # Read RMSD values
        df = read_rmsd_values(input_file, args.quiet)
        if df is None:
            continue
        
        # Filter by reference if specified
        if args.reference and 'reference' in df.columns:
            df = df[df['reference'] == args.reference]
            if len(df) == 0:
                print(f"No data found for reference '{args.reference}' in {input_file}")
                continue
        
        # Determine output file
        if args.output:
            output_file = Path(args.output)
        else:
            # Get reference name for the filename
            reference_name = template_name
            if 'reference' in df.columns and not df['reference'].empty:
                reference_name = df['reference'].iloc[0]
            
            # Save in plots directory
            output_file = plots_dir / f'rmsd_heatmap_{reference_name}.png'
        
        # Create heatmap
        success = create_heatmap(df, output_file, config, args.quiet)
        
        if success:
            # Save the data to CSV for further analysis
            csv_file = csv_dir / f'rmsd_values_{reference_name}.csv'
            df.to_csv(csv_file)
            if not args.quiet:
                print(f"Data saved to {csv_file}")
                print(f"Heatmap created successfully: {output_file}")
            return True
        else:
            print(f"Failed to create heatmap for {input_file}")
    
    return False

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
        csv_dir = Path(config["directories"]["csv"])
    else:
        # New configuration structure
        plots_dir = Path(config.get("plots", "plots"))
        csv_dir = Path(config.get("csv", "csv"))
    
    # Create plots directory if it doesn't exist
    if not plots_dir.exists():
        plots_dir.mkdir()
        if not args.quiet:
            print(f"Created directory: {plots_dir}")
    
    # Create csv directory if it doesn't exist
    if not csv_dir.exists():
        csv_dir.mkdir()
        if not args.quiet:
            print(f"Created directory: {csv_dir}")
    
    # Get all templates from configuration
    # First try to get templates from the merged config
    templates = get_templates(config)
    
    # If no templates found, try to get them from the global config
    if not templates and "global" in full_config and "templates" in full_config["global"]:
        templates = get_templates(full_config["global"])
    if not templates:
        print("No templates found in configuration.")
        return
    
    if not args.quiet:
        print(f"Found {len(templates)} templates: {', '.join(str(t) for t in templates)}")
    
    # Process each template
    templates_processed = 0
    for template in templates:
        template_name = template.stem
        if process_template(template_name, config, args, plots_dir, csv_dir):
            templates_processed += 1
    
    if templates_processed > 0:
        print(f"Successfully processed {templates_processed} out of {len(templates)} templates.")
    else:
        print("No templates were successfully processed.")

if __name__ == "__main__":
    main()
