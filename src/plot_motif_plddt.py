#!/usr/bin/env python3
"""
Script to generate a heatmap visualization of motif-specific pLDDT values.

This script reads motif-specific pLDDT values from CSV files and creates heatmaps
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
    parser = argparse.ArgumentParser(description='Generate heatmaps of motif-specific pLDDT values.')
    parser.add_argument('--motif', type=str, default=None,
                        help='Specific motif to plot (default: all motifs)')
    
    # Add common arguments, excluding the motif argument
    parser = config_loader.add_common_args(parser, exclude=['motif'])
    
    return parser.parse_args()

def create_motif_heatmap(df, motif_id, output_file, config, quiet=False):
    """Create a heatmap visualization of motif-specific pLDDT values."""
    # Filter for the specific motif
    motif_df = df[df['motif'] == motif_id]
    
    # Check if DataFrame is valid
    if motif_df.empty:
        print(f"No data for motif {motif_id}")
        return False
    
    # Get motif definition
    motif_def = config_loader.get_motif_definition(config, motif_id)
    motif_desc = motif_id
    if motif_def:
        motif_desc = motif_def.get("description", motif_id)
    
    # Pivot the data for the heatmap
    pivot_df = motif_df.pivot(index='ligand', columns='method', values='plddt')
    
    # Determine available methods based on configuration and data
    available_methods = []
    all_methods = ['chai', 'chai_with_MSA', 'boltz', 'boltz_with_MSA']
    
    # Filter methods based on configuration
    if not config.get("methods", {}).get("use_chai", True):
        all_methods = [m for m in all_methods if not m.startswith('chai')]
    if not config.get("methods", {}).get("use_boltz", True):
        all_methods = [m for m in all_methods if not m.startswith('boltz')]
    if not config.get("methods", {}).get("use_msa", True):
        all_methods = [m for m in all_methods if not '_with_MSA' in m]
    
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
    
    # Use the RdYlGn colormap (Red-Yellow-Green)
    # For pLDDT, higher values are better, so we use the standard orientation
    # (green for high values, red for low values)
    
    # Create the heatmap
    ax = sns.heatmap(pivot_df, annot=True, cmap="RdYlGn", fmt=".3f", linewidths=.5)
    
    # Set title and labels
    plt.title(f"Motif-Specific pLDDT Values: {motif_desc}", fontsize=14)
    plt.xlabel("Method", fontsize=12)
    plt.ylabel("Ligand", fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha="right")
    
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
    config = config_loader.load_config()
    config = config_loader.update_config_from_args(config, args)
    
    # Determine input directory
    csv_dir = Path(config["directories"]["csv"])
    if not csv_dir.exists():
        print(f"Error: CSV directory {csv_dir} does not exist")
        return
    
    # Determine output directory
    plots_dir = Path(config["directories"]["plots"])
    plots_dir.mkdir(exist_ok=True)
    
    # Find all motif pLDDT CSV files
    if args.motif:
        csv_files = list(csv_dir.glob(f'motif_plddt_{args.motif}.csv'))
    else:
        csv_files = list(csv_dir.glob('motif_plddt_*.csv'))
    
    if not csv_files:
        print(f"Error: No motif pLDDT CSV files found in {csv_dir}")
        return
    
    if not args.quiet:
        print(f"Found {len(csv_files)} motif pLDDT CSV files")
    
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
        output_file = plots_dir / f'motif_plddt_heatmap_{motif_id}.png'
        create_motif_heatmap(combined_df, motif_id, output_file, config, args.quiet)

if __name__ == "__main__":
    main()
