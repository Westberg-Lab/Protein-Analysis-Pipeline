# Protein Prediction Pipeline: Detailed Script Documentation

This document provides detailed documentation for each Python script in the protein prediction pipeline.

## Table of Contents

1. [run_pipeline.py](#run_pipelinepy)
2. [config_loader.py](#config_loaderpy)
3. [archive_and_clean.py](#archive_and_cleanpy)
4. [generate_chai_fasta.py](#generate_chai_fastapy)
5. [generate_boltz_yaml.py](#generate_boltz_yamlpy)
6. [run_chai_apptainer.py](#run_chai_apptainerpy)
7. [run_boltz_apptainer.py](#run_boltz_apptainerpy)
8. [combine_cif_files.py](#combine_cif_filespy)
9. [plot_rmsd_heatmap.py](#plot_rmsd_heatmappy)
10. [plot_plddt_heatmap.py](#plot_plddt_heatmappy)

---

## run_pipeline.py

### Purpose
The master script that orchestrates the entire protein prediction pipeline from start to finish.

### Functionality
- Archives previous outputs
- Runs each script in the pipeline in the correct order
- Handles errors and provides status updates
- Supports resuming from failures

### Command-line Arguments
```
python run_pipeline.py [--config CONFIG_FILE] [--no-archive] [--skip-step STEP]
                      [--use-chai] [--no-chai] [--use-boltz] [--no-boltz]
                      [--use-msa] [--no-msa] [--use-msa-dir]
                      [--template TEMPLATE_FILE] [--model-idx N]
                      [--quiet] [--resume] [--force-resume] [--state-file STATE_FILE]
                      [--clean-state]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `compute_config_hash()`: Computes a hash of the configuration to detect changes
- `read_state_file()`: Reads the pipeline state file
- `write_state_file()`: Writes the pipeline state file
- `update_state()`: Updates the pipeline state after a step
- `log_message()`: Logs a message with timestamp
- `run_step()`: Runs a pipeline step with error handling and state tracking
- `main()`: Main function that orchestrates the pipeline

### Dependencies
- `config_loader.py`: For loading and updating configuration
- All other scripts in the pipeline

### Example Usage
```bash
# Run the entire pipeline with default settings
python run_pipeline.py

# Run with a specific configuration file
python run_pipeline.py --config custom_config.json

# Skip the archiving step
python run_pipeline.py --no-archive

# Skip specific steps
python run_pipeline.py --skip-step chai-fasta --skip-step boltz-yaml

# Use only CHAI for predictions
python run_pipeline.py --use-chai --no-boltz

# Resume from a previous run
python run_pipeline.py --resume
```

---

## config_loader.py

### Purpose
Loads and manages configuration for the protein prediction pipeline.

### Functionality
- Loads configuration from a JSON file
- Provides default values if the file is not found
- Updates configuration with command-line arguments
- Adds common arguments to an ArgumentParser

### Key Functions
- `load_config(config_file='pipeline_config.json')`: Loads configuration from a JSON file
- `update_config_from_args(config, args)`: Updates configuration with command-line arguments
- `add_common_args(parser)`: Adds common arguments to an ArgumentParser

### Default Configuration
```json
{
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
    "use_chai": true,
    "use_boltz": true,
    "use_msa": true,
    "use_msa_dir": false
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
```

### Example Usage
```python
# Load configuration
config = config_loader.load_config()

# Update configuration with command-line arguments
config = config_loader.update_config_from_args(config, args)

# Add common arguments to an ArgumentParser
parser = argparse.ArgumentParser()
parser = config_loader.add_common_args(parser)
```

---

## archive_and_clean.py

### Purpose
Archives previous outputs and creates fresh directories for new runs.

### Functionality
- Creates a timestamped archive directory
- Moves previous output directories and files to the archive
- Copies configuration files (molecules.json and pipeline_config.json) to the archive
- Creates fresh empty directories for a new run
- Can delete previous outputs without archiving if requested

### Command-line Arguments
```
python archive_and_clean.py [--no-archive]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `create_archive_directory()`: Creates a timestamped archive directory
- `is_file_empty()`: Checks if a file is empty (zero bytes)
- `is_dir_empty()`: Checks if a directory is empty or contains only empty files and directories
- `archive_directories()`: Moves non-empty directories to the archive directory
- `archive_files()`: Moves non-empty files to the archive directory
- `copy_config_files()`: Copies configuration files to the archive directory
- `delete_directories()`: Deletes directories without archiving
- `delete_files()`: Deletes files without archiving
- `create_project_directories()`: Creates all necessary project directories
- `main()`: Main function that orchestrates the archiving and cleaning process

### Directories and Files Handled
- Directories: "CHAI_FASTA", "BOLTZ_YAML", "OUTPUT", "PSE_FILES", "plots", "csv"
- Files to archive: "rmsd_values.csv", "plddt_values.csv", "rmsd_heatmap.png", "plddt_heatmap.png"
- Configuration files to copy: "molecules.json", "pipeline_config.json"

### Example Usage
```bash
# Archive previous outputs and create fresh directories
python archive_and_clean.py

# Delete previous outputs without archiving
python archive_and_clean.py --no-archive
```

---

## generate_chai_fasta.py

### Purpose
Generates FASTA files for protein combinations based on molecule definitions.

### Functionality
- Loads molecule definitions from a JSON file
- Generates FASTA files for combinations of molecules
- Creates a directory structure for organizing the FASTA files

### Key Functions
- `load_molecules(json_file="molecules.json")`: Loads molecule definitions from a JSON file
- `generate_fasta_files(output_base_dir="CHAI_FASTA")`: Generates FASTA files for combinations

### Input Format
The script expects a JSON file with the following structure:
```json
{
  "molecule_1": [
    ["protein", "molecule_name_1", "SEQUENCE_1"],
    ["protein", "molecule_name_2", "SEQUENCE_2"],
    ...
  ],
  "molecule_2": [
    ["protein", "molecule_name_3", "SEQUENCE_3"],
    ["protein", "molecule_name_4", "SEQUENCE_4"],
    ...
  ]
}
```

### Output Format
The script generates FASTA files with the following format:
```
>protein|name=molecule_name_1
SEQUENCE_1
>protein|name=molecule_name_3
SEQUENCE_3
```

### Example Usage
```bash
# Generate FASTA files with default settings
python generate_chai_fasta.py
```

---

## generate_boltz_yaml.py

### Purpose
Generates YAML files for protein combinations based on molecule definitions, following the Boltz-1 schema.

### Functionality
- Loads molecule definitions from a JSON file
- Generates YAML files for combinations of molecules
- Creates a directory structure for organizing the YAML files

### Key Functions
- `load_molecules(json_file="molecules.json")`: Loads molecule definitions from a JSON file
- `generate_yaml_files(output_base_dir="BOLTZ_YAML", use_msa=False)`: Generates YAML files for combinations

### Input Format
The script expects a JSON file with the following structure:
```json
{
  "molecule_1": [
    ["protein", "molecule_name_1", "SEQUENCE_1"],
    ["protein", "molecule_name_2", "SEQUENCE_2"],
    ...
  ],
  "molecule_2": [
    ["protein", "molecule_name_3", "SEQUENCE_3"],
    ["protein", "molecule_name_4", "SEQUENCE_4"],
    ...
  ]
}
```

### Output Format
The script generates YAML files with the following format:
```yaml
sequences:
  - protein:
      id: A
      sequence: SEQUENCE_1
      msa: empty
  - protein:
      id: B
      sequence: SEQUENCE_3
      msa: empty
```

### Example Usage
```bash
# Generate YAML files with default settings
python generate_boltz_yaml.py
```

---

## run_chai_apptainer.py

### Purpose
Runs CHAI protein structure prediction using Apptainer containers.

### Functionality
- Processes each FASTA file in the input directory
- Runs the CHAI prediction tool on each FASTA file
- Supports MSA-based predictions
- Skips files that have already been processed

### Command-line Arguments
```
python run_chai_apptainer.py [--input INPUT_DIR] [--output OUTPUT_DIR] 
                            [--use-msa] [--use-msa-dir] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_msa_config(use_msa, use_msa_dir)`: Gets MSA configuration based on command-line arguments
- `run_apptainer_commands(input_dir='CHAI_FASTA', output_dir='OUTPUT/CHAI', use_msa=False, use_msa_dir=False, quiet=False)`: Runs apptainer commands for each FASTA file in the input directory
- `main()`: Main function that orchestrates the CHAI prediction process

### Apptainer Command
```bash
apptainer run --nv /emcc/westberg/shared/containers/chai.sif input_paths=<fasta_file> outdir=<output_dir> [msa_server=1] [msa_directory=CHAI_MSAs]
```

### Example Usage
```bash
# Run CHAI predictions with default settings
python run_chai_apptainer.py

# Run CHAI predictions with MSA
python run_chai_apptainer.py --use-msa

# Run CHAI predictions with MSA directory
python run_chai_apptainer.py --use-msa --use-msa-dir

# Run CHAI predictions with custom input and output directories
python run_chai_apptainer.py --input custom_fasta --output custom_output
```

---

## run_boltz_apptainer.py

### Purpose
Runs BOLTZ protein structure prediction using Apptainer containers.

### Functionality
- Processes each YAML file in the input directory
- Runs the BOLTZ prediction tool on each YAML file
- Supports MSA-based predictions
- Skips files that have already been processed

### Command-line Arguments
```
python run_boltz_apptainer.py [--input INPUT_DIR] [--output OUTPUT_DIR] 
                             [--use-msa] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_msa_config(use_msa)`: Gets MSA configuration based on command-line arguments
- `run_apptainer_commands(input_dir='BOLTZ_YAML', output_dir='OUTPUT/BOLTZ', use_msa=False, quiet=False)`: Runs apptainer commands for each YAML file in the input directory
- `main()`: Main function that orchestrates the BOLTZ prediction process

### Apptainer Command
```bash
apptainer run --nv /emcc/westberg/shared/containers/boltz.sif <yaml_file> --out_dir=<output_dir> [--use_msa_server]
```

### Example Usage
```bash
# Run BOLTZ predictions with default settings
python run_boltz_apptainer.py

# Run BOLTZ predictions with MSA
python run_boltz_apptainer.py --use-msa

# Run BOLTZ predictions with custom input and output directories
python run_boltz_apptainer.py --input custom_yaml --output custom_output
```

---

## combine_cif_files.py

### Purpose
Creates PyMOL session (.pse) files for each unique directory name that exists in both CHAI and BOLTZ outputs.

### Functionality
- Finds protein structures from both CHAI and BOLTZ outputs
- Loads structures into PyMOL
- Aligns structures to templates
- Calculates RMSD values
- Creates PyMOL session files (.pse)
- Saves RMSD values to CSV files

### Command-line Arguments
```
python combine_cif_files.py [--model-idx N] [--template TEMPLATE_FILE]
                           [--chai-output CHAI_DIR] [--boltz-output BOLTZ_DIR]
                           [--pse-files PSE_DIR] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `find_unique_names(chai_dir, boltz_dir, config, quiet=False)`: Finds all unique directory names that exist in both CHAI and BOLTZ outputs
- `sanitize_name(name)`: Sanitizes a name for use in PyMOL by replacing problematic characters
- `find_cif_file(base_dir, name, with_msa, model_idx, quiet=False)`: Finds the CIF file in the specified directory
- `create_pse_files(unique_names, chai_dir, boltz_dir, template_file, model_idx, output_dir, config, quiet=False)`: Creates .pse files for each unique name
- `get_templates(config, args)`: Gets all template files from configuration or command-line arguments
- `process_template(template_file, unique_names, chai_dir, boltz_dir, model_idx, output_dir, csv_dir, config, quiet=False)`: Processes a single template file
- `main()`: Main function that orchestrates the process

### Output Files
- PyMOL session files (.pse) for each unique name
- CSV files with RMSD values

### Example Usage
```bash
# Create PyMOL session files with default settings
python combine_cif_files.py

# Create PyMOL session files with a specific template
python combine_cif_files.py --template custom_template.cif

# Create PyMOL session files with a specific model index
python combine_cif_files.py --model-idx 0
```

---

## plot_rmsd_heatmap.py

### Purpose
Generates a heatmap visualization of RMSD values.

### Functionality
- Reads RMSD values from CSV files
- Creates a heatmap with methods on the x-axis, ligand names on the y-axis, and RMSD values as cell values
- Saves the heatmap as a PNG file
- Saves the data to a CSV file for further analysis

### Command-line Arguments
```
python plot_rmsd_heatmap.py [--input INPUT_CSV] [--output OUTPUT_PNG]
                           [--reference REFERENCE] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_templates(config)`: Gets all template files from configuration
- `find_rmsd_csv_files(pse_dir, template_name=None)`: Finds all rmsd_values.csv files in PSE_FILES subdirectories
- `read_rmsd_values(csv_file, quiet=False)`: Reads RMSD values from a CSV file
- `create_heatmap(df, output_file, config, quiet=False)`: Creates a heatmap visualization of RMSD values
- `process_template(template_name, config, args, plots_dir, csv_dir)`: Processes RMSD data for a specific template
- `main()`: Main function that orchestrates the process

### Output Files
- PNG file with the heatmap visualization
- CSV file with the RMSD values

### Example Usage
```bash
# Generate a heatmap with default settings
python plot_rmsd_heatmap.py

# Generate a heatmap with a specific input CSV file
python plot_rmsd_heatmap.py --input custom_rmsd.csv

# Generate a heatmap with a specific output PNG file
python plot_rmsd_heatmap.py --output custom_heatmap.png

# Generate a heatmap for a specific reference
python plot_rmsd_heatmap.py --reference KOr_w_momSalB
```

---

## plot_plddt_heatmap.py

### Purpose
Generates a heatmap visualization of pLDDT (predicted local distance difference test) values.

### Functionality
- Recursively searches through the subfolders within OUTPUT/CHAI and OUTPUT/BOLTZ to find all the JSON files containing pLDDT values
- Extracts the complex_plddt values
- Creates a heatmap visualization
- Saves the heatmap as a PNG file
- Saves the data to a CSV file for further analysis

### Command-line Arguments
```
python plot_plddt_heatmap.py [--chai-output CHAI_DIR] [--boltz-output BOLTZ_DIR]
                            [--output OUTPUT_PNG] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `find_chai_json_files(root_dir, quiet=False)`: Finds all outs.json files in the CHAI output directory
- `find_boltz_json_files(root_dir, quiet=False)`: Finds all confidence_*_model_0.json files in the BOLTZ output directory
- `get_chai_identifier_from_path(file_path, chai_dir)`: Extracts identifier from CHAI file path
- `get_boltz_identifier_from_path(file_path, boltz_dir)`: Extracts identifier from BOLTZ file path
- `extract_chai_plddt(file_path)`: Extracts complex_plddt from a CHAI outs.json file
- `extract_boltz_plddt(file_path)`: Extracts complex_plddt from a BOLTZ confidence_*_model_0.json file
- `organize_data(chai_files, boltz_files, chai_dir, boltz_dir, config)`: Extracts and organizes complex_plddt values from both CHAI and BOLTZ files
- `create_heatmap(data_dict, output_file, config, quiet=False)`: Creates a heatmap visualization of the complex_plddt values
- `main()`: Main function that orchestrates the process

### Output Files
- PNG file with the heatmap visualization
- CSV file with the pLDDT values

### Example Usage
```bash
# Generate a heatmap with default settings
python plot_plddt_heatmap.py

# Generate a heatmap with a specific output PNG file
python plot_plddt_heatmap.py --output custom_heatmap.png

# Generate a heatmap with custom CHAI and BOLTZ output directories
python plot_plddt_heatmap.py --chai-output custom_chai --boltz-output custom_boltz
