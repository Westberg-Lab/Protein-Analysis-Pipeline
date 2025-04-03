#!/usr/bin/env python3
"""
Script to extract pLDDT values for specific motifs from prediction outputs.

This script:
1. Finds all JSON files with pLDDT values
2. Extracts pLDDT values for the specified motif residues
3. Calculates average pLDDT values for the motif regions
"""

import os
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Extract pLDDT values for specific motifs.')
    parser.add_argument('--motif', type=str, required=True,
                        help='Motif ID to extract pLDDT values for')
    parser.add_argument('--pse-files', type=str,
                        help='PSE files directory (default: from config)')
    
    # Add common arguments, excluding the motif argument and pse-files
    parser = config_loader.add_common_args(parser, exclude=['motif', 'pse-files'])
    
    return parser.parse_args()

def find_chai_json_files(root_dir, molecule, quiet=False):
    """Find all outs.json files in the CHAI output directory for a specific molecule."""
    json_files = []
    if not quiet:
        print(f"Searching for outs.json files for molecule {molecule} in {root_dir}...")
    
    if not root_dir.exists():
        print(f"Directory {root_dir} does not exist.")
        return json_files
    
    # Look for directories that contain the molecule name
    for path in root_dir.glob(f"*{molecule}*"):
        if path.is_dir():
            # Look for outs.json files in subdirectories
            for json_file in path.rglob('outs.json'):
                json_files.append(json_file)
    
    if not quiet:
        print(f"Found {len(json_files)} outs.json files for molecule {molecule}")
    return json_files

def find_boltz_json_files(root_dir, molecule, quiet=False):
    """Find all confidence_*_model_0.json files in the BOLTZ output directory for a specific molecule."""
    json_files = []
    if not quiet:
        print(f"Searching for confidence_*_model_0.json files for molecule {molecule} in {root_dir}...")
    
    if not root_dir.exists():
        print(f"Directory {root_dir} does not exist.")
        return json_files
    
    # Look for directories that contain the molecule name
    for path in root_dir.glob(f"*{molecule}*"):
        if path.is_dir():
            # Look for boltz_results directories
            for results_dir in path.glob('boltz_results_*'):
                if results_dir.is_dir():
                    # Look for confidence files in predictions subdirectories
                    predictions_dir = results_dir / 'predictions'
                    if predictions_dir.exists() and predictions_dir.is_dir():
                        for subdir in predictions_dir.iterdir():
                            if subdir.is_dir():
                                confidence_file = subdir / f"confidence_{subdir.name}_model_0.json"
                                if confidence_file.exists():
                                    json_files.append(confidence_file)
    
    if not quiet:
        print(f"Found {len(json_files)} confidence_*_model_0.json files for molecule {molecule}")
    return json_files

def extract_chai_plddt_for_motif(file_path, motif_def):
    """Extract pLDDT values for motif residues from a CHAI outs.json file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            
            # Debug: Print keys in the JSON
            print(f"CHAI JSON keys: {list(data.keys())}")
            
            # Find the first candidate key (usually starts with "cand_")
            cand_key = None
            for key in data.keys():
                if key.startswith("cand_"):
                    cand_key = key
                    break
            
            if cand_key:
                print(f"CHAI: Found candidate key: {cand_key}")
                print(f"CHAI {cand_key} keys: {list(data[cand_key].keys())}")
                
                # Get the per_atom_plddt nested array
                per_atom_plddt = data.get(cand_key, {}).get("per_atom_plddt", [])
                
                # Flatten the nested array if it exists
                if per_atom_plddt and isinstance(per_atom_plddt, list) and len(per_atom_plddt) > 0:
                    # Check if it's a nested array
                    if isinstance(per_atom_plddt[0], list):
                        plddt_array = per_atom_plddt[0]  # Take the first array in the nested structure
                    else:
                        plddt_array = per_atom_plddt
                else:
                    plddt_array = []
            else:
                print(f"CHAI: No candidate key found in {file_path}")
                plddt_array = []
            
            # Debug: Print pLDDT array length
            print(f"CHAI pLDDT array length: {len(plddt_array)}")
            
            if not plddt_array:
                print(f"CHAI: No pLDDT array found in {file_path}")
                return None
            
            # Debug: Try to extract all pLDDT values
            if plddt_array:
                avg_all_plddt = sum(plddt_array) / len(plddt_array)
                print(f"CHAI: Average pLDDT for all residues: {avg_all_plddt}")
            
            # Extract pLDDT values for the motif residues
            motif_residues = motif_def.get("residues", [])
            print(f"CHAI: Motif residues: {motif_residues}")
            print(f"CHAI: Motif has whole_protein flag: {motif_def.get('whole_protein', False)}")
            
            # If this is a whole protein motif, use all residues
            if motif_def.get("whole_protein", False):
                print(f"CHAI: Using all residues for whole protein motif")
                motif_plddt_values = plddt_array
            else:
                motif_plddt_values = []
                
                for residue in motif_residues:
                    # Adjust for 0-based indexing if needed
                    idx = residue - 1  # Assuming residue numbers start from 1
                    if 0 <= idx < len(plddt_array):
                        motif_plddt_values.append(plddt_array[idx])
                    else:
                        print(f"CHAI: Residue {residue} (idx {idx}) is out of bounds (array length: {len(plddt_array)})")
            
            if not motif_plddt_values:
                print(f"CHAI: No valid pLDDT values found for motif residues")
                return None
            
            print(f"CHAI: Found {len(motif_plddt_values)} valid pLDDT values for motif residues")
            
            # Calculate average pLDDT for the motif
            avg_plddt = sum(motif_plddt_values) / len(motif_plddt_values)
            
            # Extract method and ligand information from the file path
            path_parts = file_path.parts
            
            # Determine if this is an MSA run
            is_msa = any("_with_MSA" in part for part in path_parts)
            method = "chai_with_MSA" if is_msa else "chai"
            
            # Extract ligand name from the parent directory
            ligand = file_path.parent.name
            
            return {
                "method": method,
                "ligand": ligand,
                "plddt": avg_plddt,
                "motif": motif_def.get("id"),
                "molecule": file_path.parent.parent.name.replace("_with_MSA", "")
            }
    
    except Exception as e:
        print(f"Error extracting pLDDT from CHAI file {file_path}: {e}")
        return None

def extract_boltz_plddt_for_motif(file_path, motif_def):
    """Extract pLDDT values for motif residues from a BOLTZ confidence_*_model_0.json file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            
            # Debug: Print keys in the JSON
            print(f"BOLTZ JSON keys: {list(data.keys())}")
            
            # Get the complex_plddt value (single value for the whole complex)
            complex_plddt = data.get("complex_plddt", 0)
            
            print(f"BOLTZ complex_plddt: {complex_plddt}")
            
            if not complex_plddt:
                print(f"BOLTZ: No complex_plddt value found in {file_path}")
                return None
            
            # For BOLTZ, we only use the complex_plddt value for all motif types
            print(f"BOLTZ: Using complex_plddt value for motif {motif_def.get('id')}")
            
            # Extract method and ligand information from the file path
            path_parts = file_path.parts
            
            # Determine if this is an MSA run
            is_msa = any("_with_MSA" in part for part in path_parts)
            method = "boltz_with_MSA" if is_msa else "boltz"
            
            # Extract ligand name from the parent directory
            ligand = file_path.parent.name
            
            # Extract molecule name from the boltz_results directory
            boltz_results_dir = None
            for part in path_parts:
                if part.startswith("boltz_results_"):
                    boltz_results_dir = part
                    break
            
            if boltz_results_dir:
                molecule_name = boltz_results_dir.replace("boltz_results_", "").replace("_with_MSA", "")
            else:
                molecule_name = "unknown"
            
            return {
                "method": method,
                "ligand": ligand,
                "plddt": complex_plddt,
                "motif": motif_def.get("id"),
                "molecule": molecule_name
            }
    
    except Exception as e:
        print(f"Error extracting pLDDT from BOLTZ file {file_path}: {e}")
        return None

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
    
    # Get molecules for this motif
    molecules = motif_def.get("molecules", [])
    if not molecules:
        print(f"Error: No molecules specified for motif '{args.motif}'")
        return
    
    # Get directories from config
    # Handle both old and new configuration structures
    if "directories" in config:
        # Old configuration structure
        chai_dir = Path(config["directories"]["chai_output"])
        boltz_dir = Path(config["directories"]["boltz_output"])
        
        # Use pse_files from command line if provided
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config["directories"]["pse_files"])
    else:
        # New configuration structure
        chai_dir = Path(config.get("chai_output", "OUTPUT/CHAI"))
        boltz_dir = Path(config.get("boltz_output", "OUTPUT/BOLTZ"))
        
        # Use pse_files from command line if provided
        if args.pse_files:
            pse_dir = Path(args.pse_files)
        else:
            pse_dir = Path(config.get("pse_files", "PSE_FILES"))
    
    # Create output directory
    pse_dir.mkdir(exist_ok=True, parents=True)
    
    # Process each molecule
    all_plddt_values = []
    
    for molecule in molecules:
        if not args.quiet:
            print(f"Processing motif {args.motif} for molecule {molecule}")
        
        # Find CHAI JSON files for this molecule
        chai_files = find_chai_json_files(chai_dir, molecule, args.quiet)
        
        # Process each CHAI file
        for file_path in chai_files:
            plddt_data = extract_chai_plddt_for_motif(file_path, motif_def)
            if plddt_data:
                all_plddt_values.append(plddt_data)
        
        # Find BOLTZ JSON files for this molecule
        boltz_files = find_boltz_json_files(boltz_dir, molecule, args.quiet)
        
        # Process each BOLTZ file
        for file_path in boltz_files:
            plddt_data = extract_boltz_plddt_for_motif(file_path, motif_def)
            if plddt_data:
                all_plddt_values.append(plddt_data)
    
    if not all_plddt_values:
        print(f"No pLDDT values generated for motif {args.motif}")
        return
    
    # Write pLDDT values to CSV in the PSE directory
    csv_file = pse_dir / f'motif_plddt_{args.motif}.csv'
    df = pd.DataFrame(all_plddt_values)
    df.to_csv(csv_file, index=False)
    
    if not args.quiet:
        print(f"Wrote {len(all_plddt_values)} pLDDT values to {csv_file}")

if __name__ == "__main__":
    main()
