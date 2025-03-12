#!/usr/bin/env python3
"""
Master script to run the entire protein prediction pipeline from start to finish.

This script:
1. Archives previous outputs (using archive_and_clean.py)
2. Runs each script in the pipeline in the correct order
3. Handles errors and provides status updates

Usage:
    python run_pipeline.py [--config CONFIG_FILE] [--no-archive] [--skip-step STEP]
                          [--use-chai] [--no-chai] [--use-boltz] [--no-boltz]
                          [--use-msa] [--no-msa] [--use-msa-dir]
                          [--template TEMPLATE_FILE] [--model-idx N]
                          [--quiet]

Options:
    --config CONFIG_FILE  Configuration file (default: pipeline_config.json)
    --no-archive          Delete previous outputs without archiving
    --skip-step STEP      Skip a specific step (chai-fasta, boltz-yaml, chai-run, 
                         boltz-run, combine-cif, rmsd-plot, plddt-plot)
    --use-chai            Use CHAI for predictions
    --no-chai             Do not use CHAI for predictions
    --use-boltz           Use BOLTZ for predictions
    --no-boltz            Do not use BOLTZ for predictions
    --use-msa             Use MSA for predictions
    --no-msa              Do not use MSA for predictions
    --use-msa-dir         Use MSA directory (for CHAI)
    --template FILE       Template file (default: from config)
    --model-idx N         Model index to search for (default: from config)
    --quiet               Suppress detailed output
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from datetime import datetime
import src.config_loader as config_loader

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run the entire protein prediction pipeline.')
    parser.add_argument('--config', type=str, default='pipeline_config.json',
                        help='Configuration file (default: pipeline_config.json)')
    parser.add_argument('--no-archive', action='store_true',
                        help='Delete previous outputs without archiving')
    parser.add_argument('--skip-step', action='append', default=[],
                        choices=['chai-fasta', 'boltz-yaml', 'chai-run', 'boltz-run', 
                                'combine-cif', 'rmsd-plot', 'plddt-plot'],
                        help='Skip a specific step in the pipeline')
    
    # Add common arguments
    parser = config_loader.add_common_args(parser)
    
    return parser.parse_args()

def log_message(message, level="INFO", quiet=False):
    """Log a message with timestamp."""
    if quiet and level == "INFO":
        return
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}")

def run_step(command, description, quiet=False):
    """Run a pipeline step with basic error handling."""
    log_message(f"Running: {description}", quiet=quiet)
    
    if not quiet:
        log_message(f"Command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True, capture_output=quiet)
        log_message(f"Completed: {description}", quiet=quiet)
        return True
    except subprocess.CalledProcessError as e:
        log_message(f"Error in {description}: {e}", level="ERROR", quiet=quiet)
        return False

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    config = config_loader.load_config(args.config)
    config = config_loader.update_config_from_args(config, args)
    
    # Step 1: Archive previous outputs
    archive_cmd = ["python", "src/archive_and_clean.py"]
    if args.no_archive:
        archive_cmd.append("--no-archive")
    
    if not run_step(archive_cmd, "Archiving previous outputs", args.quiet):
        log_message("Failed to archive previous outputs. Exiting.", level="ERROR", quiet=args.quiet)
        return 1
    
    # Define pipeline steps with standardized arguments
    pipeline_steps = []
    
    # Only add CHAI steps if CHAI is enabled
    if config["methods"]["use_chai"]:
        pipeline_steps.append({
            "id": "chai-fasta",
            "command": ["python", "src/generate_chai_fasta.py"],
            "description": "Generating CHAI FASTA files"
        })
        
        # Add CHAI run step with appropriate arguments
        chai_run_cmd = [
            "python", "src/run_chai_apptainer.py",
            f"--input={config['directories']['chai_fasta']}",
            f"--output={config['directories']['chai_output']}"
        ]
        
        # Add MSA options if enabled
        if config["methods"]["use_msa"]:
            chai_run_cmd.append("--use-msa")
            if config["methods"]["use_msa_dir"]:
                chai_run_cmd.append("--use-msa-dir")
        
        # Add quiet option if specified
        if args.quiet:
            chai_run_cmd.append("--quiet")
        
        pipeline_steps.append({
            "id": "chai-run",
            "command": chai_run_cmd,
            "description": "Running CHAI predictions"
        })
    
    # Only add BOLTZ steps if BOLTZ is enabled
    if config["methods"]["use_boltz"]:
        pipeline_steps.append({
            "id": "boltz-yaml",
            "command": ["python", "src/generate_boltz_yaml.py"],
            "description": "Generating Boltz YAML files"
        })
        
        # Add BOLTZ run step with appropriate arguments
        boltz_run_cmd = [
            "python", "src/run_boltz_apptainer.py",
            f"--input={config['directories']['boltz_yaml']}",
            f"--output={config['directories']['boltz_output']}"
        ]
        
        # Add MSA option if enabled
        if config["methods"]["use_msa"]:
            boltz_run_cmd.append("--use-msa")
        
        # Add quiet option if specified
        if args.quiet:
            boltz_run_cmd.append("--quiet")
        
        pipeline_steps.append({
            "id": "boltz-run",
            "command": boltz_run_cmd,
            "description": "Running BOLTZ predictions"
        })
    
    # Add analysis steps with appropriate arguments
    combine_cif_cmd = [
        "python", "src/combine_cif_files.py",
        f"--chai-output={config['directories']['chai_output']}",
        f"--boltz-output={config['directories']['boltz_output']}",
        f"--pse-files={config['directories']['pse_files']}",
        f"--model-idx={config['templates']['model_idx']}"
    ]
    
    # Add template if specified in command line
    if args.template:
        combine_cif_cmd.append(f"--template={args.template}")
    
    # Add method options
    if not config["methods"]["use_chai"]:
        combine_cif_cmd.append("--no-chai")
    if not config["methods"]["use_boltz"]:
        combine_cif_cmd.append("--no-boltz")
    if not config["methods"]["use_msa"]:
        combine_cif_cmd.append("--no-msa")
    
    # Add quiet option if specified
    if args.quiet:
        combine_cif_cmd.append("--quiet")
    
    pipeline_steps.append({
        "id": "combine-cif",
        "command": combine_cif_cmd,
        "description": "Combining CIF files and creating PyMOL sessions"
    })
    
    # Add RMSD plot step
    rmsd_plot_cmd = [
        "python", "src/plot_rmsd_heatmap.py",
        f"--pse-files={config['directories']['pse_files']}",
        f"--plots={config['directories']['plots']}"
    ]
    
    # Add method options
    if not config["methods"]["use_chai"]:
        rmsd_plot_cmd.append("--no-chai")
    if not config["methods"]["use_boltz"]:
        rmsd_plot_cmd.append("--no-boltz")
    if not config["methods"]["use_msa"]:
        rmsd_plot_cmd.append("--no-msa")
    
    # Add quiet option if specified
    if args.quiet:
        rmsd_plot_cmd.append("--quiet")
    
    pipeline_steps.append({
        "id": "rmsd-plot",
        "command": rmsd_plot_cmd,
        "description": "Generating RMSD heatmaps"
    })
    
    # Add pLDDT plot step
    plddt_plot_cmd = [
        "python", "src/plot_plddt_heatmap.py",
        f"--chai-output={config['directories']['chai_output']}",
        f"--boltz-output={config['directories']['boltz_output']}",
        f"--plots={config['directories']['plots']}"
    ]
    
    # Add method options
    if not config["methods"]["use_chai"]:
        plddt_plot_cmd.append("--no-chai")
    if not config["methods"]["use_boltz"]:
        plddt_plot_cmd.append("--no-boltz")
    if not config["methods"]["use_msa"]:
        plddt_plot_cmd.append("--no-msa")
    
    # Add quiet option if specified
    if args.quiet:
        plddt_plot_cmd.append("--quiet")
    
    pipeline_steps.append({
        "id": "plddt-plot",
        "command": plddt_plot_cmd,
        "description": "Generating pLDDT heatmaps"
    })
    
    # Run each step in the pipeline
    log_message("Starting pipeline execution", quiet=args.quiet)
    
    for step in pipeline_steps:
        if step["id"] in args.skip_step:
            log_message(f"Skipping: {step['description']} (--skip-step {step['id']})", quiet=args.quiet)
            continue
        
        if not run_step(step["command"], step["description"], args.quiet):
            log_message(f"Pipeline failed at step: {step['id']}", level="ERROR", quiet=args.quiet)
            return 1
    
    log_message("Pipeline completed successfully!", quiet=args.quiet)
    return 0

if __name__ == "__main__":
    sys.exit(main())
