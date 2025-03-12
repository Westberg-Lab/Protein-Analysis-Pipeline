#!/usr/bin/env python3
"""
Script to run CHAI protein structure prediction using Apptainer containers.

This script processes each FASTA file in the CHAI_FASTA directory and runs
the CHAI prediction tool on it.

Usage:
    python run_chai_apptainer.py [--input INPUT_DIR] [--output OUTPUT_DIR] 
                                [--use-msa] [--use-msa-dir] [--quiet]

Options:
    --input INPUT_DIR     Input directory containing FASTA files (default: CHAI_FASTA)
    --output OUTPUT_DIR   Output directory for predictions (default: OUTPUT/CHAI)
    --use-msa             Use MSA for predictions
    --use-msa-dir         Use MSA directory (CHAI_MSAs)
    --quiet               Suppress detailed output
"""

import subprocess
import argparse
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run CHAI protein structure prediction.')
    parser.add_argument('--input', type=str, default='CHAI_FASTA',
                        help='Input directory containing FASTA files (default: CHAI_FASTA)')
    parser.add_argument('--output', type=str, default='OUTPUT/CHAI',
                        help='Output directory for predictions (default: OUTPUT/CHAI)')
    parser.add_argument('--use-msa', action='store_true',
                        help='Use MSA for predictions')
    parser.add_argument('--use-msa-dir', action='store_true',
                        help='Use MSA directory (CHAI_MSAs)')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output')
    return parser.parse_args()

def get_msa_config(use_msa, use_msa_dir):
    """Get MSA configuration based on command line arguments."""
    if not use_msa:
        return [], False
    
    if use_msa_dir:
        return ["msa_directory=CHAI_MSAs", "msa_server=1"], True
    
    return ["msa_server=1"], True

def run_apptainer_commands(input_dir='CHAI_FASTA', output_dir='OUTPUT/CHAI', 
                          use_msa=False, use_msa_dir=False, quiet=False):
    """Run apptainer commands for each FASTA file in input directory."""
    msa_configs, using_msa = get_msa_config(use_msa, use_msa_dir)
    
    if not quiet:
        print(f"MSA config is: {msa_configs}")
    
    # Create output directory
    output_base = Path(output_dir)
    output_base.mkdir(parents=True, exist_ok=True)
    
    for folder in Path(input_dir).iterdir():
        if folder.is_dir():
            # Create subdirectory for each input folder with _with_MSA suffix if using MSA
            folder_name = f"{folder.name}_with_MSA" if using_msa else folder.name
            output_dir = output_base / folder_name
            
            output_dir.mkdir(exist_ok=True)
            
            for fasta in folder.glob('*.fasta'):
                # Check if output files already exist for this FASTA
                base_name = fasta.stem  # Get filename without extension
                output_files_exist = False
                
                # Check for typical output files
                expected_files = [
                    output_dir / base_name / "outs.json",
                    output_dir / base_name / "pred.model_idx_0.cif"
                ]
                
                if all(f.exists() for f in expected_files):
                    print(f"Skipping {fasta.name} - output files already exist")
                    continue
                
                cmd = [
                    "apptainer", "run", "--nv",
                    "/emcc/westberg/shared/containers/chai.sif",
                    f"input_paths={fasta}",
                    f"outdir=OUTPUT/CHAI/{folder_name}",
                ] + msa_configs

                if not quiet:
                    print(f"Running subprocess terminal command:\n{cmd}")
                    print(f"Processing: {fasta}")
                
                subprocess.run(cmd, capture_output=quiet)
                
                if not quiet:
                    print(f"Completed: {fasta}\n")

def main():
    """Main function."""
    args = parse_arguments()
    run_apptainer_commands(
        input_dir=args.input,
        output_dir=args.output,
        use_msa=args.use_msa,
        use_msa_dir=args.use_msa_dir,
        quiet=args.quiet
    )
    print("All FASTA files processed!")

if __name__ == "__main__":
    main()
