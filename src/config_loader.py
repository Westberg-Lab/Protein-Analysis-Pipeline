#!/usr/bin/env python3
"""
Configuration loader for the protein prediction pipeline.

This module loads the configuration from the pipeline_config.json file
and provides default values if the file is not found.
"""

import json
import argparse
import copy
from pathlib import Path

def load_config(config_file='pipeline_config.json'):
    """Load configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
            
            # Check if this is the new format with prediction_runs and analysis_runs
            if "global" in config and "prediction_runs" in config:
                return config
            # Check if this is the intermediate format with global and configurations
            elif "global" in config and "configurations" in config:
                # Convert to new format
                return {
                    "global": config["global"],
                    "prediction_runs": config["configurations"],
                    "analysis_runs": [
                        {
                            "id": "whole_protein",
                            "description": "Standard whole protein analysis",
                            "enabled": True,
                            "source_predictions": [cfg["id"] for cfg in config["configurations"] if cfg.get("enabled", True)],
                            "analysis_type": "whole_protein"
                        }
                    ]
                }
            else:
                # Convert old format to new format
                return {
                    "global": config,
                    "prediction_runs": [
                        {
                            "id": "default",
                            "description": "Default configuration",
                            "enabled": True,
                            "methods": config.get("methods", {})
                        }
                    ],
                    "analysis_runs": [
                        {
                            "id": "whole_protein",
                            "description": "Standard whole protein analysis",
                            "enabled": True,
                            "source_predictions": ["default"],
                            "analysis_type": "whole_protein"
                        }
                    ]
                }
    except FileNotFoundError:
        # Return default configuration if file not found
        return {
            "global": {
                "directories": {
                    "chai_fasta": "CHAI_FASTA",
                    "boltz_yaml": "BOLTZ_YAML",
                    "chai_output": "OUTPUT/CHAI",
                    "boltz_output": "OUTPUT/BOLTZ",
                    "pse_files": "PSE_FILES",
                    "plots": "plots",
                    "csv": "csv"
                },
                "templates": {
                    "default_template": "KOr_w_momSalB.cif",
                    "model_idx": 4
                },
                "visualization": {
                    "rmsd_vmin": 0.2,
                    "rmsd_vmax": 6.2
                }
            },
            "prediction_runs": [
                {
                    "id": "standard",
                    "description": "Standard run without MSA",
                    "enabled": True,
                    "methods": {
                        "use_chai": True,
                        "use_boltz": True,
                        "use_msa": False,
                        "use_msa_dir": False
                    }
                },
                {
                    "id": "with_msa",
                    "description": "Run with MSA enabled",
                    "enabled": True,
                    "methods": {
                        "use_chai": True,
                        "use_boltz": True,
                        "use_msa": True,
                        "use_msa_dir": False
                    }
                }
            ],
            "analysis_runs": [
                {
                    "id": "whole_protein",
                    "description": "Standard whole protein analysis",
                    "enabled": True,
                    "source_predictions": ["standard", "with_msa"],
                    "analysis_type": "whole_protein"
                }
            ]
        }

def deep_merge(base, override):
    """
    Deep merge two dictionaries.
    
    Values in override will overwrite values in base.
    For dictionaries, the merge is recursive.
    """
    result = copy.deepcopy(base)
    
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    
    return result

def get_merged_config(config, prediction_id=None):
    """
    Get a merged configuration by combining global settings with a specific prediction run.
    
    If prediction_id is None, returns the first enabled prediction run.
    If prediction_id is specified, returns that specific prediction run.
    """
    global_config = config.get("global", {})
    prediction_runs = config.get("prediction_runs", [])
    
    if not prediction_runs:
        return global_config
    
    if prediction_id:
        # Find the specific prediction run
        for run in prediction_runs:
            if run.get("id") == prediction_id:
                return deep_merge(global_config, run)
        
        # If not found, use the first enabled prediction run
        print(f"Warning: Prediction run '{prediction_id}' not found. Using first enabled prediction run.")
    
    # Get the first enabled prediction run
    for run in prediction_runs:
        if run.get("enabled", True):
            return deep_merge(global_config, run)
    
    # If no enabled prediction runs, use the first one
    return deep_merge(global_config, prediction_runs[0])

def get_enabled_prediction_runs(config, prediction_ids=None):
    """
    Get a list of enabled prediction runs.
    
    If prediction_ids is provided, only include those specific prediction runs.
    """
    prediction_runs = config.get("prediction_runs", [])
    
    if prediction_ids:
        # Filter to only the specified prediction runs
        prediction_id_list = prediction_ids.split(',')
        filtered_runs = [run for run in prediction_runs if run.get("id") in prediction_id_list]
        if not filtered_runs:
            print(f"Warning: No prediction runs found with IDs: {prediction_ids}. Using all enabled prediction runs.")
            filtered_runs = [run for run in prediction_runs if run.get("enabled", True)]
        return filtered_runs
    else:
        # Return all enabled prediction runs
        return [run for run in prediction_runs if run.get("enabled", True)]

def get_enabled_analysis_runs(config, analysis_ids=None):
    """
    Get a list of enabled analysis runs.
    
    If analysis_ids is provided, only include those specific analysis runs.
    """
    analysis_runs = config.get("analysis_runs", [])
    
    if analysis_ids:
        # Filter to only the specified analysis runs
        analysis_id_list = analysis_ids.split(',')
        filtered_runs = [run for run in analysis_runs if run.get("id") in analysis_id_list]
        if not filtered_runs:
            print(f"Warning: No analysis runs found with IDs: {analysis_ids}. Using all enabled analysis runs.")
            filtered_runs = [run for run in analysis_runs if run.get("enabled", True)]
        return filtered_runs
    else:
        # Return all enabled analysis runs
        return [run for run in analysis_runs if run.get("enabled", True)]

def get_motif_definition(config, motif_id):
    """
    Get a motif definition by ID.
    
    Returns None if the motif is not found.
    """
    motifs = config.get("global", {}).get("motifs", {}).get("definitions", [])
    for motif in motifs:
        if motif.get("id") == motif_id:
            return motif
    return None

def get_motif_template(config, motif_id, templates_dir=None):
    """
    Get the template file for a motif by ID.
    
    Returns None if the motif is not found or has no template.
    
    Args:
        config: The configuration dictionary
        motif_id: The ID of the motif
        templates_dir: Optional directory to prepend to relative template paths
    
    Returns:
        Path or None: Path to the template file if found, None otherwise
    """
    motif_def = get_motif_definition(config, motif_id)
    if not motif_def or "template" not in motif_def:
        return None
    
    template_path = Path(motif_def["template"])
    
    # If it's a relative path and templates_dir is provided, prepend templates_dir
    if templates_dir and not template_path.is_absolute():
        template_path = Path(templates_dir) / template_path
    
    return template_path

def get_prediction_run_by_id(config, prediction_id):
    """
    Get a prediction run by ID.
    
    Returns None if the prediction run is not found.
    """
    prediction_runs = config.get("prediction_runs", [])
    for run in prediction_runs:
        if run.get("id") == prediction_id:
            return run
    return None

def update_config_from_args(config, args):
    """Update configuration with command-line arguments."""
    # Create a dictionary from args
    args_dict = vars(args)
    
    # Update directories
    if "directories" in config:
        for key in config["directories"]:
            arg_key = key.replace('_', '-')
            if arg_key in args_dict and args_dict[arg_key] is not None:
                config["directories"][key] = args_dict[arg_key]
    
    # Update methods
    if "methods" in config:
        for key in config["methods"]:
            arg_key = key.replace('_', '-')
            if arg_key in args_dict:
                config["methods"][key] = args_dict[arg_key]
    
    # Update templates
    if "templates" in config:
        if "model_idx" in args_dict and args_dict["model_idx"] is not None:
            config["templates"]["model_idx"] = args_dict["model_idx"]
        if "template" in args_dict and args_dict["template"] is not None:
            if "default_template" in config["templates"]:
                config["templates"]["default_template"] = args_dict["template"]
            elif "files" in config["templates"]:
                config["templates"]["files"] = [args_dict["template"]]
    
    return config

def add_common_args(parser, exclude=None):
    """Add common arguments to an ArgumentParser."""
    exclude = exclude or []
    
    # Directory arguments
    if 'chai-fasta' not in exclude:
        parser.add_argument('--chai-fasta', type=str,
                            help=f'CHAI FASTA directory (default: CHAI_FASTA)')
    if 'boltz-yaml' not in exclude:
        parser.add_argument('--boltz-yaml', type=str,
                            help=f'BOLTZ YAML directory (default: BOLTZ_YAML)')
    if 'chai-output' not in exclude:
        parser.add_argument('--chai-output', type=str,
                            help=f'CHAI output directory (default: OUTPUT/CHAI)')
    if 'boltz-output' not in exclude:
        parser.add_argument('--boltz-output', type=str,
                            help=f'BOLTZ output directory (default: OUTPUT/BOLTZ)')
    if 'pse-files' not in exclude:
        parser.add_argument('--pse-files', type=str,
                            help=f'PSE files directory (default: PSE_FILES)')
    if 'plots' not in exclude:
        parser.add_argument('--plots', type=str,
                            help=f'Plots directory (default: plots)')
    if 'csv' not in exclude:
        parser.add_argument('--csv', type=str,
                            help=f'CSV files directory (default: csv)')
    
    # Method arguments
    if 'use-chai' not in exclude:
        parser.add_argument('--use-chai', action='store_true',
                            help='Use CHAI for predictions')
        parser.add_argument('--no-chai', dest='use_chai', action='store_false',
                            help='Do not use CHAI for predictions')
    if 'use-boltz' not in exclude:
        parser.add_argument('--use-boltz', action='store_true',
                            help='Use BOLTZ for predictions')
        parser.add_argument('--no-boltz', dest='use_boltz', action='store_false',
                            help='Do not use BOLTZ for predictions')
    if 'use-msa' not in exclude:
        parser.add_argument('--use-msa', action='store_true',
                            help='Use MSA for predictions')
        parser.add_argument('--no-msa', dest='use_msa', action='store_false',
                            help='Do not use MSA for predictions')
    if 'use-msa-dir' not in exclude:
        parser.add_argument('--use-msa-dir', action='store_true',
                            help='Use MSA directory (for CHAI)')
    
    # Template arguments
    if 'template' not in exclude:
        parser.add_argument('--template', type=str,
                            help=f'Template file (default: from config)')
    if 'model-idx' not in exclude:
        parser.add_argument('--model-idx', type=int,
                            help=f'Model index to search for (default: from config)')
    
    # Configuration arguments
    if 'prediction-runs' not in exclude:
        parser.add_argument('--prediction-runs', type=str,
                            help='Comma-separated list of prediction run IDs to run (default: all enabled)')
    if 'analysis-runs' not in exclude:
        parser.add_argument('--analysis-runs', type=str,
                            help='Comma-separated list of analysis run IDs to run (default: all enabled)')
    if 'enable-prediction' not in exclude:
        parser.add_argument('--enable-prediction', action='append', default=[],
                            help='Enable a specific prediction run (can be used multiple times)')
    if 'disable-prediction' not in exclude:
        parser.add_argument('--disable-prediction', action='append', default=[],
                            help='Disable a specific prediction run (can be used multiple times)')
    if 'enable-analysis' not in exclude:
        parser.add_argument('--enable-analysis', action='append', default=[],
                            help='Enable a specific analysis run (can be used multiple times)')
    if 'disable-analysis' not in exclude:
        parser.add_argument('--disable-analysis', action='append', default=[],
                            help='Disable a specific analysis run (can be used multiple times)')
    
    # Backwards compatibility
    if 'configs' not in exclude:
        parser.add_argument('--configs', type=str,
                            help='Comma-separated list of configuration IDs to run (default: all enabled)')
    if 'enable-config' not in exclude:
        parser.add_argument('--enable-config', action='append', default=[],
                            help='Enable a specific configuration (can be used multiple times)')
    if 'disable-config' not in exclude:
        parser.add_argument('--disable-config', action='append', default=[],
                            help='Disable a specific configuration (can be used multiple times)')
    
    # Motif arguments
    if 'motif' not in exclude:
        parser.add_argument('--motif', type=str,
                            help='Motif ID to use for analysis')
    
    # Other arguments
    if 'quiet' not in exclude:
        parser.add_argument('--quiet', action='store_true',
                            help='Suppress detailed output')
    
    return parser
