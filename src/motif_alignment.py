#!/usr/bin/env python3
"""
Script to perform motif-specific alignment and RMSD calculation.

This script:
1. Takes protein structures and motif definitions
2. Extracts the specified motif regions
3. Performs alignment on just those regions
4. Calculates motif-specific RMSD values
"""

import argparse
import pymol
from pymol import cmd
import pandas as pd
from pathlib import Path
import json
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Perform motif-specific alignment and RMSD calculation.')
    parser.add_argument('--motif', type=str, required=True,
                        help='Motif ID to use for alignment')
    
    # Add common arguments, excluding the motif argument
    parser = config_loader.add_common_args(parser, exclude=['motif'])
    
    return parser.parse_args()

def find_pse_files_for_molecule(pse_dir, molecule_name):
    """Find PSE files that contain the specified molecule."""
    pse_files = []
    
    # Look in all template directories
    for template_dir in pse_dir.iterdir():
        if template_dir.is_dir():
            # Check each PSE file
            for pse_file in template_dir.glob('*.pse'):
                # If the molecule name is in the PSE filename, include it
                if molecule_name in pse_file.stem:
                    pse_files.append(pse_file)
    
    return pse_files

def process_pse_file(pse_file, motif_def, molecule, output_dir, quiet=False):
    """Process a PyMOL session file for motif-specific alignment."""
    # Reset PyMOL
    cmd.reinitialize()
    
    # Load the PyMOL session
    cmd.load(str(pse_file))
    
    if not quiet:
        print(f"Loaded PyMOL session: {pse_file}")
    
    # Get all objects in the session
    objects = cmd.get_names('objects')
    
    # Identify the template object (usually starts with the template name)
    template_obj = None
    for obj in objects:
        if not (obj.startswith('chai_') or obj.startswith('boltz_')):
            template_obj = obj
            break
    
    if not template_obj:
        print(f"Warning: Could not identify template object in {pse_file}")
        return None
    
    # Create a selection for the motif residues
    residue_list = ','.join(str(r) for r in motif_def.get("residues", []))
    chain = motif_def.get("chain", "A")
    motif_selection = f"chain {chain} and resi {residue_list}"
    
    # Create a dictionary to store RMSD values
    rmsd_values = []
    
    # Process each prediction object
    for obj in objects:
        if obj == template_obj:
            continue
        
        # Extract method and ligand name from object name
        if obj.startswith('chai_msa_'):
            method = 'chai_with_MSA'
            ligand = obj[len('chai_msa_'):]
        elif obj.startswith('chai_'):
            method = 'chai'
            ligand = obj[len('chai_'):]
        elif obj.startswith('boltz_msa_'):
            method = 'boltz_with_MSA'
            ligand = obj[len('boltz_msa_'):]
        elif obj.startswith('boltz_'):
            method = 'boltz'
            ligand = obj[len('boltz_'):]
        else:
            continue
        
        # Perform motif-specific alignment
        try:
            # Create temporary selections for the motif regions
            cmd.select("template_motif", f"{template_obj} and {motif_selection}")
            cmd.select("pred_motif", f"{obj} and {motif_selection}")
            
            # Check if selections are valid
            if cmd.count_atoms("template_motif") == 0 or cmd.count_atoms("pred_motif") == 0:
                print(f"Warning: Empty motif selection for {obj}")
                continue
            
            # Align the motif regions
            alignment_result = cmd.align("pred_motif", "template_motif")
            rmsd = alignment_result[0]  # First element is RMSD
            
            # Store the RMSD value
            rmsd_values.append({
                'ligand': ligand,
                'method': method,
                'rmsd': rmsd,
                'motif': motif_def.get("id"),
                'molecule': molecule,
                'reference': Path(pse_file).parent.name
            })
            
            if not quiet:
                print(f"  Aligned {obj} to {template_obj} using motif {motif_def.get('id')}, RMSD: {rmsd:.4f}")
        
        except pymol.CmdException as e:
            print(f"  Error aligning {obj}: {e}")
    
    # Clean up selections
    cmd.delete("template_motif")
    cmd.delete("pred_motif")
    
    return rmsd_values

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    config = config_loader.load_config()
    config = config_loader.update_config_from_args(config, args)
    
    # Get motif definition
    motif_def = config_loader.get_motif_definition(config, args.motif)
    if not motif_def:
        print(f"Error: Motif '{args.motif}' not found in configuration")
        return
    
    # Get molecules for this motif
    molecules = motif_def.get("molecules", [])
    if not molecules:
        print(f"Error: No molecules specified for motif '{args.motif}'")
        return
    
    # Determine PSE files directory
    pse_dir = Path(config["directories"]["pse_files"])
    if not pse_dir.exists():
        print(f"Error: PSE files directory {pse_dir} does not exist")
        return
    
    # Determine output directory
    csv_dir = Path(config["directories"]["csv"])
    csv_dir.mkdir(exist_ok=True)
    
    # Process each molecule
    all_rmsd_values = []
    for molecule in molecules:
        if not args.quiet:
            print(f"Processing motif {args.motif} for molecule {molecule}")
        
        # Find PSE files for this molecule
        pse_files = find_pse_files_for_molecule(pse_dir, molecule)
        
        if not pse_files:
            print(f"Warning: No PSE files found for molecule {molecule}")
            continue
        
        if not args.quiet:
            print(f"Found {len(pse_files)} PSE files for molecule {molecule}")
        
        # Process each PSE file
        for pse_file in pse_files:
            rmsd_values = process_pse_file(pse_file, motif_def, molecule, csv_dir, args.quiet)
            if rmsd_values:
                all_rmsd_values.extend(rmsd_values)
    
    if not all_rmsd_values:
        print(f"No RMSD values generated for motif {args.motif}")
        return
    
    # Write RMSD values to CSV
    csv_file = csv_dir / f'motif_rmsd_{args.motif}.csv'
    df = pd.DataFrame(all_rmsd_values)
    df.to_csv(csv_file, index=False)
    
    if not args.quiet:
        print(f"Wrote {len(all_rmsd_values)} RMSD values to {csv_file}")

if __name__ == "__main__":
    main()
