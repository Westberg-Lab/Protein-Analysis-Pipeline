#!/usr/bin/env python3
"""
Configuration loader for the protein prediction pipeline.

This module loads the configuration from the pipeline_config.json file
and provides default values if the file is not found.
"""

import json
import argparse
from pathlib import Path

def load_config(config_file='pipeline_config.json'):
    """Load configuration from JSON file."""
    try:
        with open(config_file, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        # Return default configuration if file not found
        return {
            "directories": {
                "chai_fasta": "CHAI_FASTA",
                "boltz_yaml": "BOLTZ_YAML",
                "chai_output": "OUTPUT/CHAI",
                "boltz_output": "OUTPUT/BOLTZ",
                "pse_files": "PSE_FILES",
                "plots": "plots",
                "csv": "csv"
            },
            "methods": {
                "use_chai": True,
                "use_boltz": True,
                "use_msa": True,
                "use_msa_dir": False
            },
            "templates": {
                "default_template": "KOr_w_momSalB.cif",
                "model_idx": 4
            },
            "visualization": {
                "rmsd_vmin": 0.2,
                "rmsd_vmax": 6.2
            }
        }

def update_config_from_args(config, args):
    """Update configuration with command-line arguments."""
    # Create a dictionary from args
    args_dict = vars(args)
    
    # Update directories
    for key in config["directories"]:
        arg_key = key.replace('_', '-')
        if arg_key in args_dict and args_dict[arg_key] is not None:
            config["directories"][key] = args_dict[arg_key]
    
    # Update methods
    for key in config["methods"]:
        arg_key = key.replace('_', '-')
        if arg_key in args_dict:
            config["methods"][key] = args_dict[arg_key]
    
    # Update templates
    if "model_idx" in args_dict and args_dict["model_idx"] is not None:
        config["templates"]["model_idx"] = args_dict["model_idx"]
    if "template" in args_dict and args_dict["template"] is not None:
        config["templates"]["default_template"] = args_dict["template"]
    
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
                        help=f'Template file (default: KOr_w_momSalB.cif)')
    parser.add_argument('--model-idx', type=int,
                        help=f'Model index to search for (default: 4)')
    
    # Other arguments
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress detailed output')
    
    return parser
