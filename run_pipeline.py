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
import json
import hashlib
from pathlib import Path
from datetime import datetime
import src.config_loader as config_loader

# Define pipeline steps
PIPELINE_STEPS = [
    'archive',
    'chai-fasta',
    'chai-run',
    'boltz-yaml',
    'boltz-run',
    'combine-cif',
    'rmsd-plot',
    'plddt-plot'
]

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run the entire protein prediction pipeline.')
    parser.add_argument('--config', type=str, default='pipeline_config.json',
                        help='Configuration file (default: pipeline_config.json)')
    parser.add_argument('--no-archive', action='store_true',
                        help='Delete previous outputs without archiving')
    parser.add_argument('--skip-step', action='append', default=[],
                        choices=PIPELINE_STEPS,
                        help='Skip a specific step in the pipeline')
    
    # Add resume options
    parser.add_argument('--resume', action='store_true',
                        help='Resume pipeline from the last failed step')
    parser.add_argument('--force-resume', action='store_true',
                        help='Force resume even if configuration has changed')
    parser.add_argument('--state-file', type=str, default='pipeline_state.json',
                        help='Pipeline state file (default: pipeline_state.json)')
    parser.add_argument('--clean-state', action='store_true',
                        help='Clean the state file before starting')
    
    # Add common arguments
    parser = config_loader.add_common_args(parser)
    
    return parser.parse_args()

def compute_config_hash(config):
    """Compute a hash of the configuration to detect changes."""
    config_str = json.dumps(config, sort_keys=True)
    return hashlib.md5(config_str.encode()).hexdigest()

def read_state_file(state_file_path):
    """Read the pipeline state file."""
    if not Path(state_file_path).exists():
        return {
            "last_run": None,
            "config_hash": None,
            "completed_steps": [],
            "failed_step": None,
            "error_message": None
        }
    
    try:
        with open(state_file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error reading state file: {e}")
        return {
            "last_run": None,
            "config_hash": None,
            "completed_steps": [],
            "failed_step": None,
            "error_message": None
        }

def write_state_file(state_file_path, state):
    """Write the pipeline state file."""
    try:
        with open(state_file_path, 'w') as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"Error writing state file: {e}")

def update_state(state, step_id, success, error_message=None):
    """Update the pipeline state after a step."""
    state["last_run"] = datetime.now().isoformat()
    
    if success:
        if step_id not in state["completed_steps"]:
            state["completed_steps"].append(step_id)
        
        # If this was the failed step, clear it
        if state["failed_step"] == step_id:
            state["failed_step"] = None
            state["error_message"] = None
    else:
        state["failed_step"] = step_id
        state["error_message"] = error_message
    
    return state

def log_message(message, level="INFO", quiet=False, state_file=None, state=None):
    """Log a message with timestamp."""
    if quiet and level == "INFO":
        return
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}")
    
    # If this is an error and we have a state file, update it
    if level == "ERROR" and state_file and state:
        state["error_message"] = message
        write_state_file(state_file, state)

def run_step(command, description, step_id, quiet=False, state_file=None, state=None):
    """Run a pipeline step with basic error handling and state tracking."""
    log_message(f"Running: {description}", quiet=quiet)
    
    if not quiet:
        log_message(f"Command: {' '.join(command)}")
    
    try:
        subprocess.run(command, check=True, capture_output=quiet)
        log_message(f"Completed: {description}", quiet=quiet)
        
        # Update state if tracking is enabled
        if state_file and state:
            state = update_state(state, step_id, True)
            write_state_file(state_file, state)
        
        return True
    except subprocess.CalledProcessError as e:
        error_message = f"Error in {description}: {e}"
        log_message(error_message, level="ERROR", quiet=quiet, state_file=state_file, state=state)
        
        # Update state if tracking is enabled
        if state_file and state:
            state = update_state(state, step_id, False, error_message)
            write_state_file(state_file, state)
        
        return False

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    config = config_loader.load_config(args.config)
    config = config_loader.update_config_from_args(config, args)
    
    # Initialize state tracking if resume is enabled
    state_file = None
    state = None
    
    if args.resume or args.clean_state:
        state_file = args.state_file
        state = read_state_file(state_file)
        
        # Clean state if requested
        if args.clean_state:
            state = {
                "last_run": datetime.now().isoformat(),
                "config_hash": compute_config_hash(config),
                "completed_steps": [],
                "failed_step": None,
                "error_message": None
            }
            write_state_file(state_file, state)
            log_message(f"Cleaned state file: {state_file}", quiet=args.quiet)
        
        # Check if configuration has changed
        if args.resume and state["config_hash"] and state["config_hash"] != compute_config_hash(config):
            if not args.force_resume:
                log_message("Configuration has changed since last run. Use --force-resume to ignore this warning.", 
                           level="ERROR", quiet=args.quiet)
                return 1
            else:
                log_message("Configuration has changed, but continuing due to --force-resume", 
                           level="WARNING", quiet=args.quiet)
        
        # Update config hash
        state["config_hash"] = compute_config_hash(config)
        write_state_file(state_file, state)
        
        if args.resume and state["completed_steps"]:
            log_message(f"Resuming pipeline. Completed steps: {', '.join(state['completed_steps'])}", quiet=args.quiet)
            if state["failed_step"]:
                log_message(f"Last failed step: {state['failed_step']}", quiet=args.quiet)
                log_message(f"Error message: {state['error_message']}", quiet=args.quiet)
    
    # Step 1: Archive previous outputs (if not resuming or if not completed)
    if not args.resume or 'archive' not in state["completed_steps"]:
        if 'archive' not in args.skip_step:
            archive_cmd = ["python", "src/archive_and_clean.py"]
            if args.no_archive:
                archive_cmd.append("--no-archive")
            
            if not run_step(archive_cmd, "Archiving previous outputs", 'archive', args.quiet, state_file, state):
                log_message("Failed to archive previous outputs. Exiting.", level="ERROR", 
                           quiet=args.quiet, state_file=state_file, state=state)
                return 1
        else:
            log_message("Skipping: Archiving previous outputs (--skip-step archive)", quiet=args.quiet)
    
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
        # Skip if explicitly requested or if already completed in a previous run
        if step["id"] in args.skip_step:
            log_message(f"Skipping: {step['description']} (--skip-step {step['id']})", quiet=args.quiet)
            continue
        
        # Skip if already completed in a previous run (when resuming)
        if args.resume and state and step["id"] in state["completed_steps"]:
            log_message(f"Skipping: {step['description']} (already completed in previous run)", quiet=args.quiet)
            continue
        
        # Run the step
        if not run_step(step["command"], step["description"], step["id"], args.quiet, state_file, state):
            log_message(f"Pipeline failed at step: {step['id']}", level="ERROR", 
                       quiet=args.quiet, state_file=state_file, state=state)
            return 1
    
    # If we get here, all steps completed successfully
    log_message("Pipeline completed successfully!", quiet=args.quiet)
    
    # Clear the failed step if we completed successfully
    if state_file and state:
        state["failed_step"] = None
        state["error_message"] = None
        write_state_file(state_file, state)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
