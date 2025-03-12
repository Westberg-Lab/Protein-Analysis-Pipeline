#!/usr/bin/env python3
"""
Master script to run the entire protein prediction pipeline from start to finish.

This script:
1. Archives previous outputs (using archive_and_clean.py)
2. Runs each script in the pipeline in the correct order
3. Handles errors and provides status updates

Usage:
    python run_pipeline.py [--no-archive] [--skip-step STEP]

Options:
    --no-archive       Delete previous outputs without archiving
    --skip-step STEP   Skip a specific step (chai-fasta, boltz-yaml, chai-run, 
                       boltz-run, combine-cif, rmsd-plot, plddt-plot)
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run the entire protein prediction pipeline.')
    parser.add_argument('--no-archive', action='store_true',
                        help='Delete previous outputs without archiving')
    parser.add_argument('--skip-step', action='append', default=[],
                        choices=['chai-fasta', 'boltz-yaml', 'chai-run', 'boltz-run', 
                                'combine-cif', 'rmsd-plot', 'plddt-plot'],
                        help='Skip a specific step in the pipeline')
    return parser.parse_args()

def log_message(message, level="INFO"):
    """Log a message with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}")

def run_command(command, description):
    """Run a command and handle errors."""
    log_message(f"Starting: {description}")
    log_message(f"Command: {' '.join(command)}")
    
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        log_message(f"Completed: {description}")
        return True
    except subprocess.CalledProcessError as e:
        log_message(f"Error in {description}: {e}", "ERROR")
        log_message(f"Command output: {e.stdout}", "ERROR")
        log_message(f"Command error: {e.stderr}", "ERROR")
        return False

def main():
    """Main function."""
    args = parse_arguments()
    
    # Step 1: Archive previous outputs
    archive_cmd = ["python", "archive_and_clean.py"]
    if args.no_archive:
        archive_cmd.append("--no-archive")
    
    if not run_command(archive_cmd, "Archiving previous outputs"):
        log_message("Failed to archive previous outputs. Exiting.", "ERROR")
        return 1
    
    # Define pipeline steps
    pipeline_steps = [
        {
            "id": "chai-fasta",
            "command": ["python", "generate_chai_fasta.py"],
            "description": "Generating CHAI FASTA files"
        },
        {
            "id": "boltz-yaml",
            "command": ["python", "generate_boltz_yaml.py"],
            "description": "Generating Boltz YAML files"
        },
        {
            "id": "chai-run",
            "command": ["python", "run_chai_apptainer.py"],
            "description": "Running CHAI predictions"
        },
        {
            "id": "boltz-run",
            "command": ["python", "run_boltz_apptainer.py"],
            "description": "Running Boltz predictions"
        },
        {
            "id": "combine-cif",
            "command": ["python", "combine_cif_files.py"],
            "description": "Combining CIF files and creating PyMOL sessions"
        },
        {
            "id": "rmsd-plot",
            "command": ["python", "plot_rmsd_heatmap.py"],
            "description": "Generating RMSD heatmaps"
        },
        {
            "id": "plddt-plot",
            "command": ["python", "plot_plddt_heatmap.py"],
            "description": "Generating pLDDT heatmaps"
        }
    ]
    
    # Run each step in the pipeline
    log_message("Starting pipeline execution")
    
    for step in pipeline_steps:
        if step["id"] in args.skip_step:
            log_message(f"Skipping: {step['description']} (--skip-step {step['id']})")
            continue
        
        if not run_command(step["command"], step["description"]):
            log_message(f"Pipeline failed at step: {step['id']}", "ERROR")
            return 1
    
    log_message("Pipeline completed successfully!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
