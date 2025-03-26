#!/usr/bin/env python3
"""
Master script to run the entire protein prediction pipeline from start to finish.

This script:
1. Archives previous outputs (using archive_and_clean.py)
2. Runs prediction steps (CHAI and BOLTZ)
3. Runs analysis steps (whole protein and motif-specific)
4. Handles errors and provides status updates

Usage:
    python run_pipeline.py [--config CONFIG_FILE] [--no-archive] [--delete-outputs] [--skip-step STEP]
                          [--use-chai] [--no-chai] [--use-boltz] [--no-boltz]
                          [--use-msa] [--no-msa] [--use-msa-dir]
                          [--template TEMPLATE_FILE] [--model-idx N]
                          [--prediction-runs PRED_IDS] [--analysis-runs ANALYSIS_IDS]
                          [--quiet]

Options:
    --config CONFIG_FILE  Configuration file (default: pipeline_config.json)
    --no-archive          Skip archiving previous outputs (keep existing files)
    --delete-outputs      Delete previous outputs without archiving
    --skip-step STEP      Skip a specific step (chai-fasta, boltz-yaml, chai-run, 
                         boltz-run, combine-cif, rmsd-plot, plddt-plot,
                         motif-align, motif-rmsd, motif-plddt)
    --use-chai            Use CHAI for predictions
    --no-chai             Do not use CHAI for predictions
    --use-boltz           Use BOLTZ for predictions
    --no-boltz            Do not use BOLTZ for predictions
    --use-msa             Use MSA for predictions
    --no-msa              Do not use MSA for predictions
    --use-msa-dir         Use MSA directory (for CHAI)
    --template FILE       Template file (default: from config)
    --model-idx N         Model index to search for (default: from config)
    --prediction-runs IDS Comma-separated list of prediction run IDs to run
    --analysis-runs IDS   Comma-separated list of analysis run IDs to run
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
    'plddt-plot',
    'motif-align',
    'motif-rmsd',
    'motif-plddt'
]

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run the entire protein prediction pipeline.')
    parser.add_argument('--config', type=str, default='pipeline_config.json',
                        help='Configuration file (default: pipeline_config.json)')
    parser.add_argument('--no-archive', action='store_true',
                        help='Skip archiving previous outputs (keep existing files)')
    parser.add_argument('--delete-outputs', action='store_true',
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

def run_prediction_steps(config, args, state_file=None, state=None, run_id=None):
    """Run prediction steps for a specific prediction run."""
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
    
    # Run each step in the pipeline
    run_desc = f" for prediction run '{run_id}'" if run_id else ""
    log_message(f"Starting prediction steps{run_desc}", quiet=args.quiet)
    
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
            return False
    
    # If we get here, all steps completed successfully
    log_message(f"Prediction steps completed successfully{run_desc}!", quiet=args.quiet)
    return True

def run_whole_protein_analysis(config, args, state_file=None, state=None, run_id=None):
    """Run whole protein analysis steps."""
    # Define pipeline steps with standardized arguments
    pipeline_steps = []
    
    # Add analysis steps with appropriate arguments
    model_idx = config.get("templates", {}).get("model_idx", 4)  # Default to 4 if not found
    combine_cif_cmd = [
        "python", "src/combine_cif_files.py",
        f"--chai-output={config['directories']['chai_output']}",
        f"--boltz-output={config['directories']['boltz_output']}",
        f"--pse-files={config['directories']['pse_files']}",
        f"--model-idx={model_idx}"
    ]
    
    # Add template if specified in command line
    if args.template:
        combine_cif_cmd.append(f"--template={args.template}")
    
    # Add method options
    if not config.get("methods", {}).get("use_chai", True):
        combine_cif_cmd.append("--no-chai")
    if not config.get("methods", {}).get("use_boltz", True):
        combine_cif_cmd.append("--no-boltz")
    if not config.get("methods", {}).get("use_msa", True):
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
    if not config.get("methods", {}).get("use_chai", True):
        rmsd_plot_cmd.append("--no-chai")
    if not config.get("methods", {}).get("use_boltz", True):
        rmsd_plot_cmd.append("--no-boltz")
    if not config.get("methods", {}).get("use_msa", True):
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
    if not config.get("methods", {}).get("use_chai", True):
        plddt_plot_cmd.append("--no-chai")
    if not config.get("methods", {}).get("use_boltz", True):
        plddt_plot_cmd.append("--no-boltz")
    if not config.get("methods", {}).get("use_msa", True):
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
    run_desc = f" for analysis run '{run_id}'" if run_id else ""
    log_message(f"Starting whole protein analysis{run_desc}", quiet=args.quiet)
    
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
            return False
    
    # If we get here, all steps completed successfully
    log_message(f"Whole protein analysis completed successfully{run_desc}!", quiet=args.quiet)
    return True

def run_motif_analysis(config, args, state_file=None, state=None, run_id=None, motif_id=None, metrics=None, full_config=None):
    """Run motif-specific analysis steps."""
    # Define pipeline steps with standardized arguments
    pipeline_steps = []
    
    # Get motif definition from full_config if provided, otherwise from config
    if full_config:
        motif_def = config_loader.get_motif_definition(full_config, motif_id)
    else:
        motif_def = config_loader.get_motif_definition(config, motif_id)
        
    if not motif_def:
        log_message(f"Motif '{motif_id}' not found in configuration", level="ERROR", quiet=args.quiet)
        return False
    
    # Create a directory for this analysis run
    pse_base_dir = Path(config["directories"]["pse_files"])
    # Create the base directory if it doesn't exist
    pse_base_dir.mkdir(exist_ok=True, parents=True)
    analysis_run_dir = pse_base_dir / run_id
    analysis_run_dir.mkdir(exist_ok=True)
    
    # Add combine_cif_files step first to create base PSE files
    combine_cif_cmd = [
        "python", "src/combine_cif_files.py",
        f"--chai-output={config['directories']['chai_output']}",
        f"--boltz-output={config['directories']['boltz_output']}",
        f"--pse-files={analysis_run_dir}",  # Use the analysis run directory
        f"--model-idx={motif_def.get('model_idx', config.get('model_idx', 4))}"
    ]
    
    # Add template if specified in motif definition
    if "template" in motif_def:
        combine_cif_cmd.append(f"--template={motif_def['template']}")
    
    # Add molecules if specified in motif definition
    if "molecules" in motif_def and motif_def["molecules"]:
        molecules_str = ",".join(motif_def["molecules"])
        combine_cif_cmd.append(f"--molecules={molecules_str}")
    
    # Add method options
    if not config.get("methods", {}).get("use_chai", True):
        combine_cif_cmd.append("--no-chai")
    if not config.get("methods", {}).get("use_boltz", True):
        combine_cif_cmd.append("--no-boltz")
    if not config.get("methods", {}).get("use_msa", True):
        combine_cif_cmd.append("--no-msa")
    
    # Add quiet option if specified
    if args.quiet:
        combine_cif_cmd.append("--quiet")
    
    pipeline_steps.append({
        "id": f"combine-cif-{motif_id}",
        "command": combine_cif_cmd,
        "description": f"Creating base PSE files for motif {motif_id}"
    })
    
    # Add motif alignment step
    motif_align_cmd = [
        "python", "src/motif_alignment.py",
        f"--motif={motif_id}",
        f"--pse-files={analysis_run_dir}"  # Use the analysis run directory
    ]
    
    # Add quiet option if specified
    if args.quiet:
        motif_align_cmd.append("--quiet")
    
    pipeline_steps.append({
        "id": f"motif-align-{motif_id}",
        "command": motif_align_cmd,
        "description": f"Performing motif-specific alignment for {motif_id}"
    })
    
    # Add motif RMSD plot step if RMSD metric is enabled
    if not metrics or "rmsd" in metrics:
        motif_rmsd_cmd = [
            "python", "src/plot_motif_rmsd.py",
            f"--motif={motif_id}",
            f"--pse-files={analysis_run_dir}"  # Use the analysis run directory
        ]
        
        # Add quiet option if specified
        if args.quiet:
            motif_rmsd_cmd.append("--quiet")
        
        pipeline_steps.append({
            "id": f"motif-rmsd-{motif_id}",
            "command": motif_rmsd_cmd,
            "description": f"Generating motif-specific RMSD heatmap for {motif_id}"
        })
    
    # Add motif pLDDT extraction and plot steps if pLDDT metric is enabled
    if not metrics or "plddt" in metrics:
        motif_plddt_cmd = [
            "python", "src/extract_motif_plddt.py",
            f"--motif={motif_id}",
            f"--pse-files={analysis_run_dir}"  # Use the analysis run directory
        ]
        
        # Add quiet option if specified
        if args.quiet:
            motif_plddt_cmd.append("--quiet")
        
        pipeline_steps.append({
            "id": f"motif-plddt-extract-{motif_id}",
            "command": motif_plddt_cmd,
            "description": f"Extracting motif-specific pLDDT values for {motif_id}"
        })
        
        motif_plddt_plot_cmd = [
            "python", "src/plot_motif_plddt.py",
            f"--motif={motif_id}",
            f"--pse-files={analysis_run_dir}"  # Use the analysis run directory
        ]
        
        # Add quiet option if specified
        if args.quiet:
            motif_plddt_plot_cmd.append("--quiet")
        
        pipeline_steps.append({
            "id": f"motif-plddt-plot-{motif_id}",
            "command": motif_plddt_plot_cmd,
            "description": f"Generating motif-specific pLDDT heatmap for {motif_id}"
        })
    
    # Run each step in the pipeline
    run_desc = f" for analysis run '{run_id}'" if run_id else ""
    log_message(f"Starting motif-specific analysis for {motif_id}{run_desc}", quiet=args.quiet)
    
    for step in pipeline_steps:
        # Skip if explicitly requested
        base_step_id = step["id"].split("-")[0] + "-" + step["id"].split("-")[1]
        if base_step_id in args.skip_step:
            log_message(f"Skipping: {step['description']} (--skip-step {base_step_id})", quiet=args.quiet)
            continue
        
        # Skip if already completed in a previous run (when resuming)
        if args.resume and state and step["id"] in state["completed_steps"]:
            log_message(f"Skipping: {step['description']} (already completed in previous run)", quiet=args.quiet)
            continue
        
        # Run the step
        if not run_step(step["command"], step["description"], step["id"], args.quiet, state_file, state):
            log_message(f"Pipeline failed at step: {step['id']}", level="ERROR", 
                       quiet=args.quiet, state_file=state_file, state=state)
            return False
    
    # If we get here, all steps completed successfully
    log_message(f"Motif-specific analysis for {motif_id} completed successfully{run_desc}!", quiet=args.quiet)
    return True

def run_analysis_steps(config, full_config, args, state_file=None, state=None, run_id=None, analysis_type=None, motif_id=None, metrics=None):
    """Run analysis steps based on the analysis type."""
    if analysis_type == "whole_protein":
        return run_whole_protein_analysis(config, args, state_file, state, run_id)
    elif analysis_type == "motif":
        # Get motif definition
        motif_def = config_loader.get_motif_definition(full_config, motif_id)
        if not motif_def:
            log_message(f"Motif '{motif_id}' not found in configuration", level="ERROR", quiet=args.quiet)
            return False
        
        return run_motif_analysis(config, args, state_file, state, run_id, motif_id, metrics, full_config)
    else:
        log_message(f"Unknown analysis type: {analysis_type}", level="ERROR", quiet=args.quiet)
        return False

def main():
    """Main function."""
    args = parse_arguments()
    
    # Load configuration
    full_config = config_loader.load_config(args.config)
    
    # Apply command-line configuration overrides for backwards compatibility
    if args.enable_config:
        for config_id in args.enable_config:
            for cfg in full_config.get("prediction_runs", []):
                if cfg.get("id") == config_id:
                    cfg["enabled"] = True
    
    if args.disable_config:
        for config_id in args.disable_config:
            for cfg in full_config.get("prediction_runs", []):
                if cfg.get("id") == config_id:
                    cfg["enabled"] = False
    
    # Apply command-line configuration overrides for prediction runs
    if args.enable_prediction:
        for pred_id in args.enable_prediction:
            for pred in full_config.get("prediction_runs", []):
                if pred.get("id") == pred_id:
                    pred["enabled"] = True
    
    if args.disable_prediction:
        for pred_id in args.disable_prediction:
            for pred in full_config.get("prediction_runs", []):
                if pred.get("id") == pred_id:
                    pred["enabled"] = False
    
    # Apply command-line configuration overrides for analysis runs
    if args.enable_analysis:
        for analysis_id in args.enable_analysis:
            for analysis in full_config.get("analysis_runs", []):
                if analysis.get("id") == analysis_id:
                    analysis["enabled"] = True
    
    if args.disable_analysis:
        for analysis_id in args.disable_analysis:
            for analysis in full_config.get("analysis_runs", []):
                if analysis.get("id") == analysis_id:
                    analysis["enabled"] = False
    
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
                "config_hash": compute_config_hash(full_config),
                "completed_steps": [],
                "failed_step": None,
                "error_message": None
            }
            write_state_file(state_file, state)
            log_message(f"Cleaned state file: {state_file}", quiet=args.quiet)
        
        # Check if configuration has changed
        if args.resume and state["config_hash"] and state["config_hash"] != compute_config_hash(full_config):
            if not args.force_resume:
                log_message("Configuration has changed since last run. Use --force-resume to ignore this warning.", 
                           level="ERROR", quiet=args.quiet)
                return 1
            else:
                log_message("Configuration has changed, but continuing due to --force-resume", 
                           level="WARNING", quiet=args.quiet)
        
        # Update config hash
        state["config_hash"] = compute_config_hash(full_config)
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
            if args.delete_outputs:
                archive_cmd.append("--delete-outputs")
            
            if not run_step(archive_cmd, "Archiving previous outputs", 'archive', args.quiet, state_file, state):
                log_message("Failed to archive previous outputs. Exiting.", level="ERROR", 
                           quiet=args.quiet, state_file=state_file, state=state)
                return 1
        else:
            log_message("Skipping: Archiving previous outputs (--skip-step archive)", quiet=args.quiet)
    
    # Step 2: Run prediction steps
    # Get enabled prediction runs
    enabled_prediction_runs = config_loader.get_enabled_prediction_runs(full_config, args.prediction_runs)
    
    if not enabled_prediction_runs:
        log_message("No enabled prediction runs found. Exiting.", level="ERROR", quiet=args.quiet)
        return 1
    
    log_message(f"Found {len(enabled_prediction_runs)} enabled prediction runs", quiet=args.quiet)
    
    # Run each enabled prediction run
    prediction_success = True
    for pred_run in enabled_prediction_runs:
        run_id = pred_run.get("id", "unknown")
        log_message(f"Running prediction run: {run_id} - {pred_run.get('description', '')}", quiet=args.quiet)
        
        # Merge global config with this prediction run
        merged_config = config_loader.deep_merge(full_config.get("global", {}), pred_run)
        
        # Update with command line arguments
        merged_config = config_loader.update_config_from_args(merged_config, args)
        
        # Run the prediction steps with this prediction run
        if not run_prediction_steps(merged_config, args, state_file, state, run_id):
            prediction_success = False
            log_message(f"Pipeline failed for prediction run: {run_id}", level="ERROR", 
                       quiet=args.quiet, state_file=state_file, state=state)
    
    # Step 3: Run analysis steps
    # Get enabled analysis runs
    enabled_analysis_runs = config_loader.get_enabled_analysis_runs(full_config, args.analysis_runs)
    
    if not enabled_analysis_runs:
        log_message("No enabled analysis runs found. Skipping analysis steps.", level="WARNING", quiet=args.quiet)
    else:
        log_message(f"Found {len(enabled_analysis_runs)} enabled analysis runs", quiet=args.quiet)
        
        # Run each enabled analysis run
        analysis_success = True
        for analysis_run in enabled_analysis_runs:
            run_id = analysis_run.get("id", "unknown")
            log_message(f"Running analysis run: {run_id} - {analysis_run.get('description', '')}", quiet=args.quiet)
            
            # Get source prediction runs
            source_prediction_ids = analysis_run.get("source_predictions", [])
            if not source_prediction_ids:
                log_message(f"No source predictions specified for analysis run: {run_id}", level="WARNING", quiet=args.quiet)
                continue
            
            # Get the first source prediction run to use as a base
            source_run = None
            for pred_id in source_prediction_ids:
                source_run = config_loader.get_prediction_run_by_id(full_config, pred_id)
                if source_run:
                    break
            
            if not source_run:
                log_message(f"No valid source prediction runs found for analysis run: {run_id}", level="ERROR", quiet=args.quiet)
                continue
            
            # Merge global config with the source prediction run
            merged_config = config_loader.deep_merge(full_config.get("global", {}), source_run)
            
            # Update with command line arguments
            merged_config = config_loader.update_config_from_args(merged_config, args)
            
            # Run the analysis steps with this analysis run
            analysis_type = analysis_run.get("analysis_type")
            motif_id = analysis_run.get("motif_id")
            metrics = analysis_run.get("metrics")
            
            if not run_analysis_steps(merged_config, full_config, args, state_file, state, run_id, analysis_type, motif_id, metrics):
                analysis_success = False
                log_message(f"Pipeline failed for analysis run: {run_id}", level="ERROR", 
                           quiet=args.quiet, state_file=state_file, state=state)
    
    # Clear the failed step if we completed successfully
    if prediction_success and state_file and state:
        state["failed_step"] = None
        state["error_message"] = None
        write_state_file(state_file, state)
    
    return 0 if prediction_success else 1

if __name__ == "__main__":
    sys.exit(main())
