"""
Script to generate a heatmap visualization of RMSD values.

This script reads RMSD values from a CSV file and creates a heatmap with:
- x-axis: methods (chai, chai_with_MSA, boltz, boltz_with_MSA)
- y-axis: ligand names (each .pse file in PSE_FILES folder)
- Cell values: RMSD values

Usage:
    python plot_rmsd_heatmap.py [--input INPUT_CSV] [--output OUTPUT_PNG]

Options:
    --input INPUT_CSV    Input CSV file with RMSD values (default: rmsd_values.csv)
    --output OUTPUT_PNG  Output PNG file for the heatmap (default: rmsd_heatmap.png)
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate a heatmap of RMSD values.')
    parser.add_argument('--input', type=str, default=None,
                        help='Input CSV file with RMSD values (default: auto-detect in PSE_FILES subdirectories)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output PNG file for the heatmap (default: same directory as input with name rmsd_heatmap.png)')
    parser.add_argument('--reference', type=str, default=None,
                        help='Reference name to filter by (default: use all references)')
    return parser.parse_args()

def find_rmsd_csv_files():
    """Find all rmsd_values.csv files in PSE_FILES subdirectories."""
    pse_dir = Path('PSE_FILES')
    if not pse_dir.exists() or not pse_dir.is_dir():
        print(f"PSE_FILES directory not found.")
        return []
    
    csv_files = []
    for subdir in pse_dir.iterdir():
        if subdir.is_dir():
            csv_file = subdir / 'rmsd_values.csv'
            if csv_file.exists():
                csv_files.append(csv_file)
    
    return csv_files

def read_rmsd_values(csv_file):
    """Read RMSD values from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        print(f"Read {len(df)} RMSD values from {csv_file}")
        return df
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return None

def create_heatmap(df, output_file):
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
        
        # Sort the methods in the desired order
        method_order = ['chai', 'chai_with_MSA', 'boltz', 'boltz_with_MSA']
        pivot_df = pivot_df.reindex(columns=method_order)
        
        # Sort the ligands (rows) alphabetically
        pivot_df = pivot_df.sort_index()
        
        # Create the heatmap
        plt.figure(figsize=(10, max(8, len(pivot_df) * 0.4)))
        
        # Use a natural color spectrum: 'RdYlGn_r' (Red-Yellow-Green reversed)
        # This is intuitive for RMSD values: green for good alignments (low RMSD),
        # yellow for medium, red for poor alignments (high RMSD)
        
        # Apply logarithmic normalization to better visualize small differences
        from matplotlib.colors import Normalize
        
        vmin = 0.2
        vmax = 6.2
        
        # Create the heatmap with seaborn using the natural colormap and log normalization
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
        
        print(f"Heatmap saved to {output_file}")
        return True
    
    except Exception as e:
        print(f"Error creating heatmap: {e}")
        return False

def main():
    """Main function."""
    args = parse_arguments()
    
    # Create plots directory if it doesn't exist
    plots_dir = Path('plots')
    if not plots_dir.exists():
        plots_dir.mkdir()
        print(f"Created directory: {plots_dir}")
    
    # Determine input file(s)
    if args.input:
        # Use specified input file
        input_files = [Path(args.input)]
        if not input_files[0].exists():
            print(f"Input file {input_files[0]} does not exist.")
            return
    else:
        # Auto-detect input files
        input_files = find_rmsd_csv_files()
        if not input_files:
            print("No RMSD CSV files found in PSE_FILES subdirectories.")
            return
        print(f"Found {len(input_files)} RMSD CSV files: {', '.join(str(f) for f in input_files)}")
    
    # Process each input file
    for input_file in input_files:
        # Read RMSD values
        df = read_rmsd_values(input_file)
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
            reference_name = "Unknown"
            if 'reference' in df.columns and not df['reference'].empty:
                reference_name = df['reference'].iloc[0]
            
            # Save in plots directory
            output_file = plots_dir / f'rmsd_heatmap_{reference_name}.png'
        
        # Create heatmap
        success = create_heatmap(df, output_file)
        
        if success:
            print(f"Heatmap created successfully: {output_file}")
        else:
            print(f"Failed to create heatmap for {input_file}")

if __name__ == "__main__":
    main()
