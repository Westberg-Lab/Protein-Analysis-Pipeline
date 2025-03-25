#!/usr/bin/env python3
"""
Script to generate a heatmap visualization of motif-specific RMSD values.

This script reads motif-specific RMSD values from CSV files and creates heatmaps
with separate visualizations for each motif.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate heatmaps of motif-specific RMSD values.')
    parser.add_argument('--motif', type=str, default=None,
                        help='Specific motif to plot (default: all motifs)')
    parser.add_argument('--reference', type=str, default=None,
                        help='Reference name to filter by (default: all references)')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Minimum value for colormap (default: from config)')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Maximum value for colormap (default: from config)')
    
    # Add common arguments, excluding the motif argument
    parser = config_loader.add_common_args(parser, exclude=['motif'])
    
    return parser.parse_args()

def create_motif_heatmap(df, motif_id, output_file, config, full_config, vmin=None, vmax=None, quiet=False):
    """Create a heatmap visualization of motif-specific RMSD values."""
    # Filter for the specific motif
    motif_df = df[df['motif'] == motif_id]
    
    # Check if DataFrame is valid
    if motif_df.empty:
        print(f"No data for motif {motif_id}")
        return False
    
    # Get motif definition
    motif_def = config_loader.get_motif_definition(full_config, motif_id)
    motif_desc = motif_id
    if motif_def:
        motif_desc = motif_def.get("description", motif_id)
    
    # Pivot the data for the heatmap
    pivot_df = motif_df.pivot(index='ligand', columns='method', values='rmsd')
    
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
        print(f"No methods available for visualization of motif {motif_id}")
        return False
    
    # Reindex with available methods
    pivot_df = pivot_df.reindex(columns=available_methods)
    
    # Sort the ligands (rows) alphabetically
    pivot_df = pivot_df.sort_index()
    
    # Create the heatmap
    plt.figure(figsize=(10, max(8, len(pivot_df) * 0.4)))
    
    # Use a natural color spectrum: 'RdYlGn_r' (Red-Yellow-Green reversed)
    from matplotlib.colors import Normalize
    
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
    plt.title(f'Motif-Specific RMSD Values: {motif_desc}', fontsize=14)
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
        print(f"Heatmap for motif {motif_id} saved to {output_file}")
    
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
    else:
        # New configuration structure
        csv_dir = Path(config.get("csv", "csv"))
        plots_dir = Path(config.get("plots", "plots"))
    
    if not csv_dir.exists():
        print(f"Error: CSV directory {csv_dir} does not exist")
        return
    
    # Create output directory
    plots_dir.mkdir(exist_ok=True)
    
    # Find all motif RMSD CSV files
    if args.motif:
        csv_files = list(csv_dir.glob(f'motif_rmsd_{args.motif}.csv'))
    else:
        csv_files = list(csv_dir.glob('motif_rmsd_*.csv'))
    
    if not csv_files:
        print(f"Error: No motif RMSD CSV files found in {csv_dir}")
        return
    
    if not args.quiet:
        print(f"Found {len(csv_files)} motif RMSD CSV files")
    
    # Read and combine all CSV files
    all_data = []
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            all_data.append(df)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
    
    if not all_data:
        print("No data found in CSV files")
        return
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Filter by reference if specified
    if args.reference:
        combined_df = combined_df[combined_df['reference'] == args.reference]
        if combined_df.empty:
            print(f"No data found for reference {args.reference}")
            return
    
    # Get unique motifs
    motifs = combined_df['motif'].unique()
    if args.motif:
        if args.motif in motifs:
            motifs = [args.motif]
        else:
            print(f"Motif {args.motif} not found in data")
            return
    
    if not args.quiet:
        print(f"Processing {len(motifs)} motifs: {', '.join(motifs)}")
    
    # Create heatmap for each motif
    for motif_id in motifs:
        output_file = plots_dir / f'motif_rmsd_heatmap_{motif_id}.png'
        create_motif_heatmap(combined_df, motif_id, output_file, config, full_config, args.vmin, args.vmax, args.quiet)

if __name__ == "__main__":
    main()
