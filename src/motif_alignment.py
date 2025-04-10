#!/usr/bin/env python3
"""
Script to calculate motif-specific RMSD values using whole-protein alignment.

This script:
1. Takes protein structures and motif definitions
2. Aligns the whole proteins
3. Extracts the specified motif regions
4. Calculates RMSD values for the motif regions using the whole-protein alignment
"""

import argparse
import pymol
from pymol import cmd
import pandas as pd
from pathlib import Path
import json
import config_loader

def get_template_for_motif(motif_def, templates_dir):
    """Get the template file for a motif if specified."""
    if "template" in motif_def:
        template_path = Path(motif_def["template"])
        # If it's a relative path, prepend templates_dir
        if not template_path.is_absolute():
            template_path = templates_dir / template_path
        return template_path
    return None

def find_pse_files_for_molecule(pse_dir, molecule_name):
    """Find PSE files that contain the specified molecule."""
    pse_files = []
    
    # Look directly in the pse_dir first
    for pse_file in pse_dir.glob('*.pse'):
        # If the molecule name is in the PSE filename, include it
        if molecule_name in pse_file.stem:
            pse_files.append(pse_file)
    
    # Also look in all template directories
    for template_dir in pse_dir.iterdir():
        if template_dir.is_dir():
            # Check each PSE file
            for pse_file in template_dir.glob('*.pse'):
                # If the molecule name is in the PSE filename, include it
                if molecule_name in pse_file.stem:
                    pse_files.append(pse_file)
    
    return pse_files

def process_pse_file(pse_file, motif_def, molecule, output_dir, quiet=False, motif_pse_dir=None):
    """Process a PyMOL session file to calculate motif-specific RMSD."""
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
    
    # Get chain from motif definition
    chain = motif_def.get("chain", "A")
    
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
        
        try:
            # Check if this is a whole protein motif
            if motif_def.get("whole_protein", False):
                # For whole protein motifs, align the entire chains and use the RMSD directly
                whole_chain = f"chain {chain}"
                cmd.select("template_whole", f"{template_obj} and {whole_chain}")
                cmd.select("pred_whole", f"{obj} and {whole_chain}")
                
                # Check if whole protein selections are valid
                if cmd.count_atoms("template_whole") == 0 or cmd.count_atoms("pred_whole") == 0:
                    print(f"Warning: Empty whole protein selection for {obj}")
                    continue
                
                # Align the whole proteins and get RMSD directly
                alignment_result = cmd.align("pred_whole", "template_whole")
                rmsd = alignment_result[0]  # First element is RMSD
                
                if not quiet:
                    print(f"  Aligned whole protein {obj} to {template_obj}, RMSD: {rmsd:.4f}")
            else:
                # For specific residue motifs, we need both residues and template_residues
                if "residues" not in motif_def or "template_residues" not in motif_def:
                    print(f"Error: Both 'residues' and 'template_residues' must be specified for motif {motif_def.get('id')}")
                    continue
                
                # Create residue lists for selection
                residue_list = ','.join(str(r) for r in motif_def.get("residues", []))
                template_residue_list = ','.join(str(r) for r in motif_def.get("template_residues", []))
                
                # Create selections for the specific residues
                cmd.select("template_motif", f"{template_obj} and chain {chain} and resi {template_residue_list}")
                cmd.select("pred_motif", f"{obj} and chain {chain} and resi {residue_list}")
                
                # Check if motif selections are valid
                if cmd.count_atoms("template_motif") == 0 or cmd.count_atoms("pred_motif") == 0:
                    print(f"Warning: Empty motif selection for {obj}")
                    continue
                
                # Align the specific residues and get RMSD directly
                alignment_result = cmd.align("pred_motif", "template_motif")
                rmsd = alignment_result[3]  # 4th element is the RMSD before refinement
                
                if not quiet:
                    print(f"  Aligned motif residues {obj} to {template_obj}, RMSD: {rmsd:.4f}")
            
            # Store the RMSD value
            rmsd_values.append({
                'ligand': ligand,
                'method': method,
                'rmsd': rmsd,
                'motif': motif_def.get("id"),
                'molecule': molecule,
                'reference': Path(pse_file).parent.name
            })
        
        except pymol.CmdException as e:
            print(f"  Error calculating RMSD for {obj}: {e}")
    
    # Save motif-aligned PSE file if requested
    if motif_pse_dir is not None:
        # Create a visually enhanced representation of the motif regions
        
        # Hide everything first
        cmd.hide("everything")
        
        # Show the full proteins as transparent cartoons for context
        cmd.show("cartoon", "all")
        cmd.set("cartoon_transparency", 0.7, "all")
        
        # Get the chain
        chain = motif_def.get("chain", "A")
        
        # Process each object (template and predictions)
        all_objects = [template_obj] + [obj for obj in objects if obj != template_obj]
        
        for obj in all_objects:
            # Create a selection for this object's motif
            motif_sel_name = f"{obj}_motif"
            
            # Check if this is a whole protein motif
            if motif_def.get("whole_protein", False):
                # Select the entire chain
                cmd.select(motif_sel_name, f"{obj} and chain {chain}")
            else:
                # Use the specified residue list (for the prediction objects) or template residue list (for the template)
                if obj == template_obj:
                    residue_list = ','.join(str(r) for r in motif_def.get("template_residues", []))
                else:
                    residue_list = ','.join(str(r) for r in motif_def.get("residues", []))
                
                cmd.select(motif_sel_name, f"{obj} and chain {chain} and resi {residue_list}")
            
            # Show as sticks only (no spheres)
            cmd.show("sticks", motif_sel_name)
            
            # No need to set colors - the motif will inherit the color from the parent object
            
            if not quiet:
                print(f"  Showing motif for {obj} as sticks with original color")
        
        # Save the aligned session - use a simpler filename without the motif ID
        motif_pse_file = motif_pse_dir / f"{pse_file.stem}_motif.pse"
        cmd.save(str(motif_pse_file))
        
        if not quiet:
            print(f"  Saved session with highlighted motif regions to {motif_pse_file}")
        
        # Clean up the selections
        for obj in all_objects:
            cmd.delete(f"{obj}_motif")
    
    # Clean up selections
    cmd.delete("template_motif")
    cmd.delete("pred_motif")
    cmd.delete("template_whole")
    cmd.delete("pred_whole")
    cmd.delete("all_motifs")
    
    return rmsd_values

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Calculate motif-specific RMSD values using whole-protein alignment.')
    parser.add_argument('--motif', type=str, required=True,
                        help='Motif ID to use for RMSD calculation after whole-protein alignment')
    parser.add_argument('--save-motif-pse', action='store_true', default=True,
                        help='Save PSE files with highlighted motif regions (default: True)')
    parser.add_argument('--no-save-motif-pse', dest='save_motif_pse', action='store_false',
                        help='Do not save PSE files with highlighted motif regions')
    parser.add_argument('--pse-files', type=str,
                        help='PSE files directory (default: from config)')
    
    # Add common arguments, excluding the motif argument and pse-files
    parser = config_loader.add_common_args(parser, exclude=['motif', 'pse-files'])
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    full_config = config_loader.load_config()
    
    # Get a merged configuration (combines global settings with a prediction run)
    config = config_loader.get_merged_config(full_config)
    
    # Update with command-line arguments
    config = config_loader.update_config_from_args(config, args)
    
    # Get motif definition
    motif_def = config_loader.get_motif_definition(full_config, args.motif)
    if not motif_def:
        print(f"Error: Motif '{args.motif}' not found in configuration")
        return
    
    # Validate motif definition
    if not motif_def.get("whole_protein", False):
        if "residues" not in motif_def:
            print(f"Error: 'residues' must be specified for motif {args.motif}")
            return
        if "template_residues" not in motif_def:
            print(f"Error: 'template_residues' must be specified for motif {args.motif}")
            return
    
    # Get molecules for this motif
    molecules = motif_def.get("molecules", [])
    if not molecules:
        print(f"Error: No molecules specified for motif '{args.motif}'")
        return
    
    # Get directories from config
    # Handle both old and new configuration structures
    if "directories" in config:
        # Old configuration structure
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config["directories"]["pse_files"])
        csv_dir = Path(config["directories"]["csv"])
        templates_dir = Path(config["directories"].get("templates", "templates"))
    else:
        # New configuration structure
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config.get("pse_files", "PSE_FILES"))
        csv_dir = Path(config.get("csv", "csv"))
        templates_dir = Path(config.get("templates_dir", "templates"))
    
    if not pse_dir.exists():
        print(f"Error: PSE files directory {pse_dir} does not exist")
        return
    
    # Create output directories
    csv_dir.mkdir(exist_ok=True)
    
    # Create directory for motif-aligned PSE files if saving is enabled
    if args.save_motif_pse:
        # Save directly in the pse_dir, not in a motifs subdirectory
        motif_pse_dir = pse_dir
        if not args.quiet:
            print(f"Using directory for PSE files with highlighted motif regions: {motif_pse_dir}")
    else:
        motif_pse_dir = None
    
    # Check if this motif has a specific template
    template_path = config_loader.get_motif_template_path(full_config, args.motif, templates_dir)
    if template_path and template_path.exists():
        if not args.quiet:
            print(f"Found motif-specific template: {template_path}")
    
    # Process each molecule
    all_rmsd_values = []
    for molecule in molecules:
        if not args.quiet:
            print(f"Processing motif {args.motif} for molecule {molecule}")
        
        # Find PSE files for this molecule
        pse_files = find_pse_files_for_molecule(pse_dir, molecule)
        
        if not pse_files:
            # If no PSE files found and we have a motif-specific template, we might need to create them
            if template_path and template_path.exists():
                print(f"No PSE files found for molecule {molecule}, but motif has a specific template.")
                print(f"You may need to run combine_cif_files.py first with the --template option.")
                print(f"Example: python src/combine_cif_files.py --template {template_path}")
            else:
                print(f"Warning: No PSE files found for molecule {molecule}")
            continue
        
        if not args.quiet:
            print(f"Found {len(pse_files)} PSE files for molecule {molecule}")
        
        # Process each PSE file
        for pse_file in pse_files:
            rmsd_values = process_pse_file(pse_file, motif_def, molecule, csv_dir, args.quiet, motif_pse_dir)
            if rmsd_values:
                all_rmsd_values.extend(rmsd_values)
    
    if not all_rmsd_values:
        print(f"No RMSD values generated for motif {args.motif}")
        return
    
    # Write RMSD values to CSV in the PSE directory
    csv_file = pse_dir / f'motif_rmsd_{args.motif}.csv'
    df = pd.DataFrame(all_rmsd_values)
    df.to_csv(csv_file, index=False)
    
    if not args.quiet:
        print(f"Wrote {len(all_rmsd_values)} RMSD values to {csv_file}")

if __name__ == "__main__":
    main()
