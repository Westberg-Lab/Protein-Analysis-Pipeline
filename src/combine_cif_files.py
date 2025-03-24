#!/usr/bin/env python3
"""
Script to create PyMOL session (.pse) files for each unique directory name that exists
in both CHAI and BOLTZ outputs. Each .pse file will include:
1. The pred.model_idx_4.cif protein structure from CHAI (with and without MSA)
2. The pred.model_idx_4.cif protein structure from BOLTZ (with and without MSA)
3. The template file from the current working directory

Usage:
    python combine_cif_files.py [--model-idx N] [--template TEMPLATE_FILE]
                               [--chai-output CHAI_DIR] [--boltz-output BOLTZ_DIR]
                               [--pse-files PSE_DIR] [--quiet]

Options:
    --model-idx N        Model index to search for (default: from config)
    --template FILE      Template file (default: from config)
    --chai-output DIR    CHAI output directory (default: from config)
    --boltz-output DIR   BOLTZ output directory (default: from config)
    --pse-files DIR      Output directory for .pse files (default: from config)
    --quiet              Suppress detailed output
"""

import os
import argparse
import pymol
from pymol import cmd
from pathlib import Path
import config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Create PyMOL session files for protein structures.')
    
    # Add common arguments
    parser = config_loader.add_common_args(parser)
    
    return parser.parse_args()

def find_unique_names(chai_dir, boltz_dir, config, quiet=False):
    """Find all unique directory names that exist in both CHAI and BOLTZ outputs."""
    # Get default method values if not present in config
    use_chai = config.get("methods", {}).get("use_chai", True)
    use_boltz = config.get("methods", {}).get("use_boltz", True)
    
    # Get all directory names in CHAI (excluding _with_MSA directories)
    chai_dirs = set()
    if use_chai and chai_dir.exists():
        for d in chai_dir.iterdir():
            if d.is_dir() and not d.name.endswith('_with_MSA'):
                chai_dirs.update(subdir.name for subdir in d.iterdir() if subdir.is_dir())
    
    # Get all directory names in BOLTZ (excluding _with_MSA directories)
    boltz_dirs = set()
    if use_boltz and boltz_dir.exists():
        for d in boltz_dir.iterdir():
            if d.is_dir() and not d.name.endswith('_with_MSA'):
                # Extract base names from boltz_results_ prefix
                boltz_dirs.update(subdir.name.replace('boltz_results_', '') 
                                for subdir in d.iterdir() 
                                if subdir.is_dir() and subdir.name.startswith('boltz_results_'))
    
    # Find common names if both CHAI and BOLTZ are used
    if use_chai and use_boltz:
        unique_names = chai_dirs.intersection(boltz_dirs)
    elif use_chai:
        unique_names = chai_dirs
    elif use_boltz:
        unique_names = boltz_dirs
    else:
        unique_names = set()
    
    if not quiet:
        print(f"Found {len(unique_names)} unique names: {', '.join(unique_names)}")
    
    return unique_names

def sanitize_name(name):
    """Sanitize a name for use in PyMOL by replacing problematic characters."""
    # Replace characters that might cause issues in PyMOL
    sanitized = name.replace('[', '_').replace(']', '_').replace('(', '_').replace(')', '_')
    return sanitized

def find_cif_file(base_dir, name, with_msa, model_idx, quiet=False):
    """Find the CIF file in the specified directory."""
    # Determine the parent directory name
    # For names like KORDshort_WT_SalB, the parent dir is KORDshort_WT
    # For names like KORDshort_SalB, the parent dir is KORDshort
    if name.startswith("KORDshort_WT"):
        parent_dir = "KORDshort_WT"
    else:
        parent_dir = "KORDshort"
    
    # Handle CHAI and BOLTZ differently
    if 'CHAI' in str(base_dir):
        # CHAI files
        if with_msa:
            search_dir = base_dir / f"{parent_dir}_with_MSA" / name
        else:
            search_dir = base_dir / parent_dir / name
        
        if not quiet:
            print(f"    Searching in: {search_dir}")
        
        # Look for the CIF file
        cif_file = search_dir / f"pred.model_idx_{model_idx}.cif"
        
        # Check if file exists
        if cif_file.exists():
            return cif_file
    else:
        # BOLTZ files have a different structure
        if with_msa:
            search_dir = base_dir / parent_dir / f"boltz_results_{name}_with_MSA"
        else:
            search_dir = base_dir / parent_dir / f"boltz_results_{name}"
        
        if not quiet:
            print(f"    Searching in: {search_dir}")
        
        # For BOLTZ, check in predictions/[name]/[name]_model_0.cif
        predictions_dir = search_dir / "predictions"
        
        # First try the expected path with name
        name_dir = predictions_dir / name
        if name_dir.exists() and name_dir.is_dir():
            cif_file = name_dir / f"{name}_model_0.cif"
            if not quiet:
                print(f"    Looking for BOLTZ file: {cif_file}")
            
            if cif_file.exists():
                return cif_file
        
        # If not found, check all subdirectories in predictions/
        if predictions_dir.exists() and predictions_dir.is_dir():
            if not quiet:
                print(f"    Checking all subdirectories in: {predictions_dir}")
            # Check each subdirectory for the CIF file
            for subdir in predictions_dir.iterdir():
                if subdir.is_dir():
                    # Try with subdirectory name
                    cif_file = subdir / f"{subdir.name}_model_0.cif"
                    if not quiet:
                        print(f"    Looking for BOLTZ file: {cif_file}")
                    if cif_file.exists():
                        return cif_file
    
    return None

def create_pse_files(unique_names, chai_dir, boltz_dir, template_file, model_idx, output_dir, config, quiet=False):
    """Create .pse files for each unique name."""
    # Get reference name from template file
    reference_name = template_file.name.split('.')[0]  # Remove file extension
    
    # List to store RMSD values
    rmsd_values = []
    
    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Create a subdirectory for this template
    template_dir = output_dir / reference_name
    template_dir.mkdir(exist_ok=True)
    if not quiet:
        print(f"Created directory for template: {template_dir}")
    
    for name in unique_names:
        if not quiet:
            print(f"Processing {name}...")
        
        # Reset PyMOL
        cmd.reinitialize()
        
        # Load template - let PyMOL assign the default name
        if template_file.exists():
            # Get current objects before loading
            objects_before = set(cmd.get_names())
            
            # Load template file
            cmd.load(str(template_file))
            if not quiet:
                print(f"  Loaded template: {template_file}")
            
            # Get new objects after loading
            objects_after = set(cmd.get_names())
            
            # Find the newly added object(s)
            new_objects = list(objects_after - objects_before)
            
            if not new_objects:
                print(f"  Warning: No objects loaded from template file")
                continue
            
            # If there's only one object, use it as the template
            if len(new_objects) == 1:
                template_obj = new_objects[0]
            else:
                # Find the largest object by atom count
                largest_obj = None
                max_atoms = 0
                for obj in new_objects:
                    atom_count = cmd.count_atoms(obj)
                    if not quiet:
                        print(f"    Object {obj} has {atom_count} atoms")
                    if atom_count > max_atoms:
                        max_atoms = atom_count
                        largest_obj = obj
                
                template_obj = largest_obj
            
            if not quiet:
                print(f"  Using {template_obj} as template for alignment (largest protein)")
            
            # Keep the original name for the template object
        else:
            print(f"  Warning: Template file {template_file} not found")
            continue
        
        # Find and load CIF files
        structures_loaded = 0
        
        # Sanitize the name for PyMOL
        sanitized_name = sanitize_name(name)
        
        # CHAI without MSA
        if config["methods"]["use_chai"]:
            chai_file = find_cif_file(chai_dir, name, False, model_idx, quiet)
            if chai_file:
                structure_name = f'chai_{sanitized_name}'
                cmd.load(str(chai_file), structure_name)
                # Align to template protein
                try:
                    alignment_result = cmd.align(structure_name, template_obj)
                    rmsd = alignment_result[0]  # First element is RMSD
                    rmsd_values.append({
                        'ligand': name,
                        'method': 'chai',
                        'rmsd': rmsd
                    })
                    if not quiet:
                        print(f"  Loaded and aligned CHAI: {chai_file}, RMSD: {rmsd:.4f}")
                    structures_loaded += 1
                except pymol.CmdException as e:
                    print(f"  Error aligning CHAI: {e}")
            elif not quiet:
                print(f"  Warning: CHAI file for {name} (without MSA) not found")
        
        # CHAI with MSA
        if config["methods"]["use_chai"] and config["methods"]["use_msa"]:
            chai_msa_file = find_cif_file(chai_dir, name, True, model_idx, quiet)
            if chai_msa_file:
                structure_name = f'chai_msa_{sanitized_name}'
                cmd.load(str(chai_msa_file), structure_name)
                # Align to template protein
                try:
                    alignment_result = cmd.align(structure_name, template_obj)
                    rmsd = alignment_result[0]  # First element is RMSD
                    rmsd_values.append({
                        'ligand': name,
                        'method': 'chai_with_MSA',
                        'rmsd': rmsd
                    })
                    if not quiet:
                        print(f"  Loaded and aligned CHAI (with MSA): {chai_msa_file}, RMSD: {rmsd:.4f}")
                    structures_loaded += 1
                except pymol.CmdException as e:
                    print(f"  Error aligning CHAI with MSA: {e}")
            elif not quiet:
                print(f"  Warning: CHAI file for {name} (with MSA) not found")
        
        # BOLTZ without MSA
        if config["methods"]["use_boltz"]:
            boltz_file = find_cif_file(boltz_dir, name, False, model_idx, quiet)
            if boltz_file:
                structure_name = f'boltz_{sanitized_name}'
                cmd.load(str(boltz_file), structure_name)
                # Align to template protein
                try:
                    alignment_result = cmd.align(structure_name, template_obj)
                    rmsd = alignment_result[0]  # First element is RMSD
                    rmsd_values.append({
                        'ligand': name,
                        'method': 'boltz',
                        'rmsd': rmsd
                    })
                    if not quiet:
                        print(f"  Loaded and aligned BOLTZ: {boltz_file}, RMSD: {rmsd:.4f}")
                    structures_loaded += 1
                except pymol.CmdException as e:
                    print(f"  Error aligning BOLTZ: {e}")
            elif not quiet:
                print(f"  Warning: BOLTZ file for {name} (without MSA) not found")
        
        # BOLTZ with MSA
        if config["methods"]["use_boltz"] and config["methods"]["use_msa"]:
            boltz_msa_file = find_cif_file(boltz_dir, name, True, model_idx, quiet)
            if boltz_msa_file:
                structure_name = f'boltz_msa_{sanitized_name}'
                cmd.load(str(boltz_msa_file), structure_name)
                # Align to template protein
                try:
                    alignment_result = cmd.align(structure_name, template_obj)
                    rmsd = alignment_result[0]  # First element is RMSD
                    rmsd_values.append({
                        'ligand': name,
                        'method': 'boltz_with_MSA',
                        'rmsd': rmsd
                    })
                    if not quiet:
                        print(f"  Loaded and aligned BOLTZ (with MSA): {boltz_msa_file}, RMSD: {rmsd:.4f}")
                    structures_loaded += 1
                except pymol.CmdException as e:
                    print(f"  Error aligning BOLTZ with MSA: {e}")
            elif not quiet:
                print(f"  Warning: BOLTZ file for {name} (with MSA) not found")
        
        # Set nice visualization
        cmd.hide('everything')
        cmd.show('cartoon')
        
        # Color all template objects the same color (cyan)
        for obj in new_objects:
            cmd.color('green', obj)
        
        # Define a list of distinct colors for the loaded structures
        colors = ['cyan', 'yellow', 'magenta', 'orange', 'pink', 'violet', 'salmon', 'lime', 'blue', 'red']
        color_index = 0
        
        # Color each loaded structure with a different color
        loaded_structures = []
        if config["methods"]["use_chai"] and 'chai_file' in locals() and chai_file:
            loaded_structures.append(f'chai_{sanitized_name}')
        if config["methods"]["use_chai"] and config["methods"]["use_msa"] and 'chai_msa_file' in locals() and chai_msa_file:
            loaded_structures.append(f'chai_msa_{sanitized_name}')
        if config["methods"]["use_boltz"] and 'boltz_file' in locals() and boltz_file:
            loaded_structures.append(f'boltz_{sanitized_name}')
        if config["methods"]["use_boltz"] and config["methods"]["use_msa"] and 'boltz_msa_file' in locals() and boltz_msa_file:
            loaded_structures.append(f'boltz_msa_{sanitized_name}')
        
        for structure in loaded_structures:
            cmd.color(colors[color_index % len(colors)], structure)
            color_index += 1
        
        if structures_loaded > 0:
            cmd.center(template_obj)
            cmd.zoom('all')
            if not quiet:
                print(f"  Colored template objects cyan and each structure with a unique color")
        
            # Save as PSE file if at least one structure was loaded
            if structures_loaded > 0:
                pse_file = template_dir / f"{name}.pse"
                cmd.save(str(pse_file))
                if not quiet:
                    print(f"  Created {pse_file}")
            else:
                print(f"  Skipping {name}.pse - No structures found")

    # Return the RMSD values, reference name, and template directory
    return rmsd_values, reference_name, template_dir

def get_templates(config, args):
    """Get all template files from configuration or command line arguments."""
    # Get templates directory from config
    if "directories" in config and "templates" in config["directories"]:
        templates_dir = Path(config["directories"]["templates"])
    else:
        templates_dir = Path("templates")  # Default to templates directory
    
    if args.template:
        # If template is specified in command line, use it directly
        template_path = Path(args.template)
        # If it's a relative path and doesn't exist, try prepending templates_dir
        if not template_path.is_absolute() and not template_path.exists():
            template_path = templates_dir / template_path
        return [template_path]
    
    # Handle both old and new configuration structures
    if "templates" in config:
        # Old configuration structure
        if "files" in config["templates"]:
            return [templates_dir / template for template in config["templates"]["files"]]
        elif "default_template" in config["templates"]:
            return [templates_dir / config["templates"]["default_template"]]
    else:
        # New configuration structure
        if "files" in config:
            return [templates_dir / template for template in config["files"]]
        elif "default_template" in config:
            return [templates_dir / config["default_template"]]
    
    return []

def process_template(template_file, unique_names, chai_dir, boltz_dir, model_idx, output_dir, csv_dir, config, quiet=False):
    """Process a single template file."""
    if not template_file.exists():
        print(f"Template file {template_file} does not exist.")
        return False
    
    if not quiet:
        print(f"\nProcessing template: {template_file}")
    
    # Create PSE files and get RMSD values, reference name, and template directory
    rmsd_values, reference_name, template_dir = create_pse_files(
        unique_names, chai_dir, boltz_dir, template_file, model_idx, output_dir, config, quiet
    )
    
    if not rmsd_values:
        print(f"No RMSD values generated for template {template_file}.")
        return False
    
    # Write RMSD values to CSV with reference protein name in both locations
    # 1. In the template directory (for backward compatibility)
    pse_csv_file = template_dir / 'rmsd_values.csv'
    with open(pse_csv_file, 'w') as f:
        f.write('ligand,method,rmsd,reference\n')
        for entry in rmsd_values:
            f.write(f"{entry['ligand']},{entry['method']},{entry['rmsd']},{reference_name}\n")
    
    # 2. In the csv directory (for centralized storage)
    csv_file = csv_dir / f'rmsd_values_{reference_name}.csv'
    with open(csv_file, 'w') as f:
        f.write('ligand,method,rmsd,reference\n')
        for entry in rmsd_values:
            f.write(f"{entry['ligand']},{entry['method']},{entry['rmsd']},{reference_name}\n")
    
    if not quiet:
        print(f"RMSD values saved to {pse_csv_file} and {csv_file}")
    
    return True

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    config = config_loader.load_config()
    config = config_loader.update_config_from_args(config, args)
    
    # Get directories from config
    # Handle both old and new configuration structures
    if "directories" in config:
        # Old configuration structure
        chai_dir = Path(config["directories"]["chai_output"])
        boltz_dir = Path(config["directories"]["boltz_output"])
        output_dir = Path(config["directories"]["pse_files"])
        csv_dir = Path(config["directories"]["csv"])
    else:
        # New configuration structure
        chai_dir = Path(config.get("chai_output", "OUTPUT/CHAI"))
        boltz_dir = Path(config.get("boltz_output", "OUTPUT/BOLTZ"))
        output_dir = Path(config.get("pse_files", "PSE_FILES"))
        csv_dir = Path(config.get("csv", "csv"))
    
    # Get model index from config
    if "templates" in config:
        model_idx = config["templates"]["model_idx"]
    else:
        model_idx = config.get("model_idx", 4)
    
    # Find unique names
    unique_names = find_unique_names(chai_dir, boltz_dir, config, args.quiet)
    
    if not unique_names:
        print("No unique names found. Please check if the directories exist and contain data.")
        return
    
    # Create csv directory if it doesn't exist
    if not csv_dir.exists():
        csv_dir.mkdir()
        if not args.quiet:
            print(f"Created directory: {csv_dir}")
    
    # Get all templates from configuration or command line arguments
    templates = get_templates(config, args)
    if not templates:
        print("No templates found in configuration or command line arguments.")
        return
    
    if not args.quiet:
        print(f"Found {len(templates)} templates: {', '.join(str(t) for t in templates)}")
    
    # Process each template
    templates_processed = 0
    for template_file in templates:
        if process_template(template_file, unique_names, chai_dir, boltz_dir, model_idx, output_dir, csv_dir, config, args.quiet):
            templates_processed += 1
    
    if templates_processed > 0:
        print(f"Successfully processed {templates_processed} out of {len(templates)} templates.")
        print(f"All PSE files created successfully!")
    else:
        print("No templates were successfully processed.")

if __name__ == "__main__":
    main()
