#!/usr/bin/env python3
"""
Unified script to generate heatmap visualizations of RMSD values.

This script can read RMSD values from CSV files and create heatmaps for:
1. Whole protein RMSD values (from rmsd_values.csv files)
2. Motif-specific RMSD values (from motif_rmsd_*.csv files)

Usage:
    python plot_rmsd.py [--input INPUT_CSV] [--output OUTPUT_PNG]
                       [--reference REFERENCE] [--motif MOTIF_ID]
                       [--vmin VMIN] [--vmax VMAX] [--quiet]

Options:
    --input INPUT_CSV    Input CSV file or directory (default: auto-detect)
    --output OUTPUT_PNG  Output PNG file for the heatmap (default: plots/rmsd_heatmap_[reference].png)
    --reference REF      Reference name to filter by (default: use all references)
    --motif MOTIF_ID     Specific motif to plot (default: all motifs)
    --vmin VMIN          Minimum value for colormap (default: from config)
    --vmax VMAX          Maximum value for colormap (default: from config)
    --quiet              Suppress detailed output
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import config_loader
from matplotlib.colors import Normalize

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate heatmaps of RMSD values.')
    parser.add_argument('--input', type=str, default=None,
                        help='Input CSV file or directory (default: auto-detect)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG file for the heatmap (default: plots/rmsd_heatmap_[reference].png)')
    parser.add_argument('--reference', type=str, default=None,
                        help='Reference name to filter by (default: use all references)')
    parser.add_argument('--motif', type=str, default=None,
                        help='Specific motif to plot (default: all motifs)')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Minimum value for colormap (default: from config)')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Maximum value for colormap (default: from config)')
    parser.add_argument('--pse-files', type=str,
                        help='PSE files directory (default: from config)')
    parser.add_argument('--analysis-run', type=str, default=None,
                        help='Analysis run ID to use (default: auto-detect from input)')
    
    # Add common arguments
    parser = config_loader.add_common_args(parser, exclude=['motif', 'pse-files'])
    
    return parser.parse_args()

def find_rmsd_csv_files(pse_dir, template_name=None, motif_id=None, quiet=False):
    """Find RMSD CSV files in the specified directory.
    
    Args:
        pse_dir: Directory to search in
        template_name: Optional template name to filter by
        motif_id: Optional motif ID to filter by
        quiet: Whether to suppress output
        
    Returns:
        List of paths to CSV files
    """
    if not pse_dir.exists() or not pse_dir.is_dir():
        print(f"{pse_dir} directory not found.")
        return []
    
    csv_files = []
    
    # If motif_id is specified, look for motif-specific RMSD files
    if motif_id:
        if not quiet:
            print(f"Looking for motif-specific RMSD files for motif {motif_id}...")
        
        # Look for motif_rmsd_*.csv files in the PSE directory and its subdirectories
        if template_name:
            # If template is specified, look in that specific directory
            csv_file = pse_dir / f"motif_rmsd_{motif_id}_{template_name}.csv"
            if csv_file.exists():
                csv_files.append(csv_file)
            
            # Also look in subdirectories
            for subdir in pse_dir.iterdir():
                if subdir.is_dir():
                    csv_file = subdir / f"motif_rmsd_{motif_id}_{template_name}.csv"
                    if csv_file.exists():
                        csv_files.append(csv_file)
        else:
            # Otherwise, look for all motif_rmsd_*.csv files in the PSE directory
            for csv_file in pse_dir.glob(f"motif_rmsd_{motif_id}*.csv"):
                csv_files.append(csv_file)
            
            # Also look in subdirectories
            for subdir in pse_dir.iterdir():
                if subdir.is_dir():
                    for csv_file in subdir.glob(f"motif_rmsd_{motif_id}*.csv"):
                        csv_files.append(csv_file)
            
            # If no specific motif files found, look for general motif files
            if not csv_files:
                for csv_file in pse_dir.glob(f"motif_rmsd_{motif_id}.csv"):
                    csv_files.append(csv_file)
                
                # Also look in subdirectories
                for subdir in pse_dir.iterdir():
                    if subdir.is_dir():
                        csv_file = subdir / f"motif_rmsd_{motif_id}.csv"
                        if csv_file.exists():
                            csv_files.append(csv_file)
    else:
        # Look for whole protein RMSD files (rmsd_values.csv)
        if template_name:
            # If template is specified, look in that specific directory
            template_dir = pse_dir / template_name
            if template_dir.exists() and template_dir.is_dir():
                csv_file = template_dir / 'rmsd_values.csv'
                if csv_file.exists():
                    csv_files.append(csv_file)
        else:
            # Otherwise, look in all subdirectories
            for subdir in pse_dir.iterdir():
                if subdir.is_dir():
                    # Check for rmsd_values.csv directly in the subdirectory
                    csv_file = subdir / 'rmsd_values.csv'
                    if csv_file.exists():
                        csv_files.append(csv_file)
                    
                    # Also check for rmsd_values.csv in subdirectories of the subdirectory
                    for sub_subdir in subdir.iterdir():
                        if sub_subdir.is_dir():
                            csv_file = sub_subdir / 'rmsd_values.csv'
                            if csv_file.exists():
                                csv_files.append(csv_file)
    
    if not quiet:
        print(f"Found {len(csv_files)} RMSD CSV files")
    
    return csv_files

def read_rmsd_values(csv_files, quiet=False):
    """Read RMSD values from CSV files.
    
    Args:
        csv_files: List of CSV files to read
        quiet: Whether to suppress output
        
    Returns:
        DataFrame with RMSD values
    """
    all_data = []
    
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            if not quiet:
                print(f"Read {len(df)} RMSD values from {csv_file}")
            all_data.append(df)
        except Exception as e:
            print(f"Error reading CSV file {csv_file}: {e}")
    
    if not all_data:
        print("No data found in CSV files")
        return None
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    return combined_df

def create_rmsd_heatmap(df, output_file, config, full_config=None, motif_id=None, vmin=None, vmax=None, quiet=False):
    """Create a heatmap visualization of RMSD values.
    
    Args:
        df: DataFrame with RMSD values
        output_file: Path to save the heatmap
        config: Configuration dictionary
        full_config: Full configuration dictionary (for motif definitions)
        motif_id: Optional motif ID to filter by
        vmin: Minimum value for colormap
        vmax: Maximum value for colormap
        quiet: Whether to suppress output
        
    Returns:
        True if successful, False otherwise
    """
    # Check if DataFrame is valid
    if df is None or len(df) == 0:
        print("No data to visualize.")
        return False
    
    # Check if required columns exist
    required_columns = ['ligand', 'method', 'rmsd']
    if not all(col in df.columns for col in required_columns):
        print(f"CSV file must contain columns: {', '.join(required_columns)}")
        return False
    
    # Filter by motif if specified
    if motif_id and 'motif' in df.columns:
        df = df[df['motif'] == motif_id]
        if len(df) == 0:
            print(f"No data found for motif '{motif_id}'")
            return False
    
    # Get reference name from the data if available
    reference_name = "Unknown"
    if 'reference' in df.columns and not df['reference'].empty:
        reference_name = df['reference'].iloc[0]
    
    # Get motif description if available
    motif_desc = motif_id
    if motif_id and full_config:
        motif_def = config_loader.get_motif_definition(full_config, motif_id)
        if motif_def:
            motif_desc = motif_def.get("description", motif_id)
    
    # Check for duplicate entries
    if df.duplicated(subset=['ligand', 'method']).any():
        if not quiet:
            print("Warning: Found duplicate entries. Aggregating by taking the mean of RMSD values.")
        # Aggregate duplicate entries by taking the mean
        df = df.groupby(['ligand', 'method'], as_index=False)['rmsd'].mean()
    
    # Pivot the data for the heatmap
    pivot_df = df.pivot(index='ligand', columns='method', values='rmsd')
    
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
    
    # Get visualization parameters
    if vmin is None:
        vmin = config.get("visualization", {}).get("rmsd_vmin", 0.2)
    if vmax is None:
        vmax = config.get("visualization", {}).get("rmsd_vmax", 6.2)
    
    # Create the heatmap with seaborn
    ax = sns.heatmap(pivot_df, annot=True, cmap='RdYlGn_r', fmt='.4f', 
                 norm=Normalize(vmin=vmin, vmax=vmax),
                 cbar_kws={'label': 'RMSD (Å)'})
    
    # Set labels and title
    if motif_id:
        plt.title(f'Motif-Specific RMSD Values: {motif_desc}', fontsize=14)
    else:
        plt.title(f'RMSD Values by Method and Ligand (Reference: {reference_name})', fontsize=14)
    
    plt.xlabel('Method', fontsize=12)
    plt.ylabel('Ligand', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout
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
        csv_dir = Path(config["directories"]["csv"])
        plots_dir = Path(config["directories"]["plots"])
        
        # Use pse_files from command line if provided
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config["directories"]["pse_files"])
    else:
        # New configuration structure
        csv_dir = Path(config.get("csv", "csv"))
        plots_dir = Path(config.get("plots", "plots"))
        
        # Use pse_files from command line if provided
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config.get("pse_files", "PSE_FILES"))
    
    # Create output directories
    plots_dir.mkdir(exist_ok=True)
    csv_dir.mkdir(exist_ok=True)
    
    # Determine input files
    if args.input:
        # Use specified input file or directory
        input_path = Path(args.input)
        if input_path.is_file():
            # Single file
            csv_files = [input_path]
        elif input_path.is_dir():
            # Directory - find CSV files
            if args.motif:
                # Look for motif-specific RMSD files
                csv_files = list(input_path.glob(f"motif_rmsd_{args.motif}*.csv"))
                if not csv_files:
                    csv_files = list(input_path.glob("motif_rmsd_*.csv"))
            else:
                # Look for whole protein RMSD files
                csv_files = []
                for subdir in input_path.iterdir():
                    if subdir.is_dir():
                        csv_file = subdir / 'rmsd_values.csv'
                        if csv_file.exists():
                            csv_files.append(csv_file)
        else:
            print(f"Input path {input_path} does not exist.")
            return
    else:
        # Auto-detect input files
        if args.motif:
            # Look for motif-specific RMSD files in the PSE directory
            csv_files = find_rmsd_csv_files(pse_dir, None, args.motif, args.quiet)
        else:
            # Look for whole protein RMSD files in the pse_dir
            csv_files = find_rmsd_csv_files(pse_dir, args.reference, None, args.quiet)
    
    if not csv_files:
        print("No RMSD CSV files found. Please check if the files exist and try again.")
        return
    
    # Read RMSD values
    df = read_rmsd_values(csv_files, args.quiet)
    if df is None:
        return
    
    # Filter by reference if specified
    if args.reference and 'reference' in df.columns:
        df = df[df['reference'] == args.reference]
        if len(df) == 0:
            print(f"No data found for reference '{args.reference}'")
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
        
        # Create analysis run subdirectory in plots directory if we have an analysis run
        if analysis_run_name:
            plot_dir = plots_dir / analysis_run_name
            plot_dir.mkdir(exist_ok=True, parents=True)
        else:
            plot_dir = plots_dir
        
        # Get reference name for the filename
        if args.motif:
            output_file = plot_dir / f'motif_rmsd_heatmap_{args.motif}.png'
        else:
            reference_name = args.reference
            if not reference_name and 'reference' in df.columns and not df['reference'].empty:
                reference_name = df['reference'].iloc[0]
            output_file = plot_dir / f'rmsd_heatmap_{reference_name}.png'
    
    # Create heatmap
    success = create_rmsd_heatmap(df, output_file, config, full_config, args.motif, args.vmin, args.vmax, args.quiet)
    
    if success:
        if not args.quiet:
            print(f"Heatmap created successfully: {output_file}")
    else:
        print(f"Failed to create heatmap")

if __name__ == "__main__":
    main()
