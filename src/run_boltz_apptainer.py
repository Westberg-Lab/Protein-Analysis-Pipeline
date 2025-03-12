#!/usr/bin/env python3
"""
Script to run Boltz protein structure prediction using Apptainer containers.

This script processes each YAML file in the BOLTZ_YAML directory and runs
the Boltz prediction tool on it.

Usage:
    python run_boltz_apptainer.py [--input INPUT_DIR] [--output OUTPUT_DIR] 
                                 [--use-msa] [--quiet]

Options:
    --input INPUT_DIR     Input directory containing YAML files (default: BOLTZ_YAML)
    --output OUTPUT_DIR   Output directory for predictions (default: OUTPUT/BOLTZ)
    --use-msa             Use MSA server for predictions
    --quiet               Suppress detailed output
"""

import subprocess
import argparse
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run Boltz protein structure prediction.')
    parser.add_argument('--input', type=str, default='BOLTZ_YAML',
                        help='Input directory containing YAML files (default: BOLTZ_YAML)')
    parser.add_argument('--output', type=str, default='OUTPUT/BOLTZ',
                        help='Output directory for predictions (default: OUTPUT/BOLTZ)')
    parser.add_argument('--use-msa', action='store_true',
                        help='Use MSA server for predictions')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output')
    return parser.parse_args()

def get_msa_config(use_msa):
    """Get MSA configuration based on command line arguments."""
    if use_msa:
        return ["--use_msa_server"], True
    return [], False

def run_apptainer_commands(input_dir='BOLTZ_YAML', output_dir='OUTPUT/BOLTZ', 
                          use_msa=False, quiet=False):
    """Run apptainer commands for each YAML file in input directory."""
    msa_config, using_msa = get_msa_config(use_msa)
    
    # Create output directory
    output_base = Path(output_dir)
    output_base.mkdir(parents=True, exist_ok=True)
    
    for folder in Path(input_dir).iterdir():
        if folder.is_dir():
            # Create subdirectory for each input folder with _with_MSA suffix if using MSA
            folder_name = f"{folder.name}_with_MSA" if using_msa else folder.name
            output_dir = output_base / folder_name
            
            output_dir.mkdir(exist_ok=True)
            
            for yaml_file in folder.glob('*.yaml'):
                # Get base name and add _with_MSA if using MSA
                base_name = yaml_file.stem  # Get filename without extension
                if using_msa:
                    base_name = f"{base_name}_with_MSA"
                
                # Check for typical output directories with boltz_results_ prefix
                expected_dirs = [
                    output_dir / f"boltz_results_{base_name}" / "predictions",
                    output_dir / f"boltz_results_{base_name}" / "processed"
                ]
                
                if all(d.exists() for d in expected_dirs):
                    print(f"Skipping {yaml_file.name} - output directories already exist")
                    continue
                
                cmd = [
                    "apptainer", "run", "--nv",
                    "/emcc/westberg/shared/containers/boltz.sif",
                    str(yaml_file),  # YAML file path as positional argument
                    f"--out_dir=OUTPUT/BOLTZ/{folder_name}"
                ] + msa_config
                
                if not quiet:
                    print(f"Processing: {yaml_file}")
                
                subprocess.run(cmd, capture_output=quiet)
                
                if not quiet:
                    print(f"Completed: {yaml_file}\n")

def main():
    """Main function."""
    args = parse_arguments()
    run_apptainer_commands(
        input_dir=args.input,
        output_dir=args.output,
        use_msa=args.use_msa,
        quiet=args.quiet
    )
    print("All YAML files processed!")

if __name__ == "__main__":
    main()
