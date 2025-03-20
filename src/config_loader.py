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
            
            # Check if this is the new format with global and configurations
            if "global" in config and "configurations" in config: # instead just remove this backwards compatibility
                return config
            else:
                # Convert old format to new format
                return {
                    "global": config,
                    "configurations": [
                        {
                            "id": "default",
                            "description": "Default configuration",
                            "enabled": True,
                            "methods": config.get("methods", {})
                        }
                    ]
                }
    except FileNotFoundError: # instead remove this and just dont run if it cannot find the config file.
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
            "configurations": [
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

def get_merged_config(config, config_id=None):
    """
    Get a merged configuration by combining global settings with a specific configuration.
    
    If config_id is None, returns the first enabled configuration.
    If config_id is specified, returns that specific configuration.
    """
    global_config = config.get("global", {})
    configurations = config.get("configurations", [])
    
    if not configurations:
        return global_config
    
    if config_id:
        # Find the specific configuration
        for cfg in configurations:
            if cfg.get("id") == config_id:
                return deep_merge(global_config, cfg)
        
        # If not found, use the first enabled configuration
        print(f"Warning: Configuration '{config_id}' not found. Using first enabled configuration.")
    
    # Get the first enabled configuration
    for cfg in configurations:
        if cfg.get("enabled", True):
            return deep_merge(global_config, cfg)
    
    # If no enabled configurations, use the first one
    return deep_merge(global_config, configurations[0])

def get_enabled_configurations(config, config_ids=None):
    """
    Get a list of enabled configurations.
    
    If config_ids is provided, only include those specific configurations.
    """
    configurations = config.get("configurations", [])
    
    if config_ids:
        # Filter to only the specified configurations
        config_id_list = config_ids.split(',')
        filtered_configs = [cfg for cfg in configurations if cfg.get("id") in config_id_list]
        if not filtered_configs:
            print(f"Warning: No configurations found with IDs: {config_ids}. Using all enabled configurations.")
            filtered_configs = [cfg for cfg in configurations if cfg.get("enabled", True)]
        return filtered_configs
    else:
        # Return all enabled configurations
        return [cfg for cfg in configurations if cfg.get("enabled", True)]

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

def add_common_args(parser):
    """Add common arguments to an ArgumentParser."""
    # Directory arguments
    parser.add_argument('--chai-fasta', type=str,
                        help=f'CHAI FASTA directory (default: CHAI_FASTA)')
    parser.add_argument('--boltz-yaml', type=str,
                        help=f'BOLTZ YAML directory (default: BOLTZ_YAML)')
    parser.add_argument('--chai-output', type=str,
                        help=f'CHAI output directory (default: OUTPUT/CHAI)')
    parser.add_argument('--boltz-output', type=str,
                        help=f'BOLTZ output directory (default: OUTPUT/BOLTZ)')
    parser.add_argument('--pse-files', type=str,
                        help=f'PSE files directory (default: PSE_FILES)')
    parser.add_argument('--plots', type=str,
                        help=f'Plots directory (default: plots)')
    parser.add_argument('--csv', type=str,
                        help=f'CSV files directory (default: csv)')
    
    # Method arguments
    parser.add_argument('--use-chai', action='store_true',
                        help='Use CHAI for predictions')
    parser.add_argument('--no-chai', dest='use_chai', action='store_false',
                        help='Do not use CHAI for predictions')
    parser.add_argument('--use-boltz', action='store_true',
                        help='Use BOLTZ for predictions')
    parser.add_argument('--no-boltz', dest='use_boltz', action='store_false',
                        help='Do not use BOLTZ for predictions')
    parser.add_argument('--use-msa', action='store_true',
                        help='Use MSA for predictions')
    parser.add_argument('--no-msa', dest='use_msa', action='store_false',
                        help='Do not use MSA for predictions')
    parser.add_argument('--use-msa-dir', action='store_true',
                        help='Use MSA directory (for CHAI)')
    
    # Template arguments
    parser.add_argument('--template', type=str,
                        help=f'Template file (default: from config)')
    parser.add_argument('--model-idx', type=int,
                        help=f'Model index to search for (default: from config)')
    
    # Configuration arguments
    parser.add_argument('--configs', type=str,
                        help='Comma-separated list of configuration IDs to run (default: all enabled)')
    parser.add_argument('--enable-config', action='append', default=[],
                        help='Enable a specific configuration (can be used multiple times)')
    parser.add_argument('--disable-config', action='append', default=[],
                        help='Disable a specific configuration (can be used multiple times)')
    
    # Other arguments
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output')
    
    return parser
