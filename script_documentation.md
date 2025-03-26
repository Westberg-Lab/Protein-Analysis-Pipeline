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
11. [motif_alignment.py](#motif_alignmentpy)
12. [extract_motif_plddt.py](#extract_motif_plddtpy)
13. [plot_motif_rmsd.py](#plot_motif_rmsdpy)
14. [plot_motif_plddt.py](#plot_motif_plddtpy)

---

## run_pipeline.py

### Purpose
The master script that orchestrates the entire protein prediction pipeline from start to finish.

### Functionality
- Archives previous outputs
- Runs prediction steps (CHAI and BOLTZ)
- Runs analysis steps (whole protein and motif-specific)
- Handles errors and provides status updates
- Supports resuming from failures
- Supports running multiple prediction and analysis runs in a single pipeline execution

### Command-line Arguments
```
python run_pipeline.py [--config CONFIG_FILE] [--no-archive] [--delete-outputs] [--skip-step STEP]
                      [--use-chai] [--no-chai] [--use-boltz] [--no-boltz]
                      [--use-msa] [--no-msa] [--use-msa-dir]
                      [--template TEMPLATE_FILE] [--model-idx N]
                      [--prediction-runs PRED_IDS] [--analysis-runs ANALYSIS_IDS]
                      [--enable-prediction PRED_ID] [--disable-prediction PRED_ID]
                      [--enable-analysis ANALYSIS_ID] [--disable-analysis ANALYSIS_ID]
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
- `run_prediction_steps()`: Runs prediction steps for a specific prediction run
- `run_whole_protein_analysis()`: Runs whole protein analysis steps
- `run_motif_analysis()`: Runs motif-specific analysis steps
- `run_analysis_steps()`: Runs analysis steps based on the analysis type
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

# Run only specific prediction runs
python run_pipeline.py --prediction-runs standard,with_msa

# Run only specific analysis runs
python run_pipeline.py --analysis-runs whole_protein_hM1D_analysis,binding_pocket_hM1D_analysis

# Enable specific prediction and analysis runs
python run_pipeline.py --enable-prediction standard --disable-prediction with_msa --enable-analysis whole_protein_hM1D_analysis

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
- Supports multiple prediction and analysis runs in a single file
- Provides default values if the file is not found
- Updates configuration with command-line arguments
- Adds common arguments to an ArgumentParser
- Retrieves motif definitions and templates

### Key Functions
- `load_config(config_file='pipeline_config.json')`: Loads configuration from a JSON file
- `deep_merge(base, override)`: Deep merges two dictionaries
- `get_merged_config(config, prediction_id=None)`: Gets a merged configuration by combining global settings with a specific prediction run
- `get_enabled_prediction_runs(config, prediction_ids=None)`: Gets a list of enabled prediction runs
- `get_enabled_analysis_runs(config, analysis_ids=None)`: Gets a list of enabled analysis runs
- `get_motif_definition(config, motif_id)`: Gets a motif definition by ID
- `get_motif_template(config, motif_id, templates_dir=None)`: Gets the template file for a motif by ID
- `get_motif_template_and_model_idx(config, motif_id, templates_dir=None)`: Gets the template file and model_idx for a motif by ID
- `get_prediction_run_by_id(config, prediction_id)`: Gets a prediction run by ID
- `update_config_from_args(config, args)`: Updates configuration with command-line arguments
- `add_common_args(parser, exclude=None)`: Adds common arguments to an ArgumentParser

### Configuration Structure
The configuration file has three main sections:

1. **Global Settings**: Common settings shared across all runs
2. **Prediction Runs**: Multiple specific prediction configurations
3. **Analysis Runs**: Multiple specific analysis configurations

```json
{
  "global": {
    "directories": {
      "chai_fasta": "CHAI_FASTA",
      "boltz_yaml": "BOLTZ_YAML",
      "chai_output": "OUTPUT/CHAI",
      "boltz_output": "OUTPUT/BOLTZ",
      "pse_files": "PSE_FILES",
      "plots": "plots",
      "csv": "csv",
      "templates": "templates"
    },
    "model_idx": 4,
    "visualization": {
      "rmsd_vmin": 0.2,
      "rmsd_vmax": 6.2
    },
    "motifs": {
      "definitions": [
        {
          "id": "binding_pocket_hM1D",
          "description": "Ligand binding pocket residues for hM1Dshort_8E9X",
          "molecules": ["hM1Dshort_8E9X"],
          "residues": [112, 113, 116, 117, 120, 164, 196, 199, 200, 203, 204, 413, 416, 417, 420, 439, 442, 443],
          "chain": "A",
          "color": "magenta",
          "template": "templates/hM1D_template_8E9X.cif",
          "model_idx": 4
        },
        {
          "id": "whole_protein_hM1D",
          "description": "Complete protein structure for hM1Dshort_8E9X",
          "molecules": ["hM1Dshort_8E9X"],
          "whole_protein": true,
          "chain": "A",
          "color": "cyan",
          "template": "templates/hM1D_template_8E9X.cif",
          "model_idx": 4
        }
      ]
    }
  },
  "prediction_runs": [
    {
      "id": "standard",
      "description": "Standard run without MSA",
      "enabled": true,
      "methods": {
        "use_chai": true,
        "use_boltz": true,
        "use_msa": false,
        "use_msa_dir": false
      }
    },
    {
      "id": "with_msa",
      "description": "Run with MSA enabled",
      "enabled": true,
      "methods": {
        "use_chai": true,
        "use_boltz": true,
        "use_msa": true,
        "use_msa_dir": false
      }
    }
  ],
  "analysis_runs": [
    {
      "id": "whole_protein_hM1D_analysis",
      "description": "Whole protein analysis for hM1D",
      "enabled": true,
      "source_predictions": ["standard", "with_msa"],
      "analysis_type": "motif",
      "motif_id": "whole_protein_hM1D",
      "metrics": ["rmsd", "plddt"]
    },
    {
      "id": "binding_pocket_hM1D_analysis",
      "description": "Binding pocket motif analysis for hM1D",
      "enabled": true,
      "source_predictions": ["standard", "with_msa"],
      "analysis_type": "motif",
      "motif_id": "binding_pocket_hM1D",
      "metrics": ["rmsd", "plddt"]
    }
  ]
}
```

### Example Usage
```python
# Load configuration
full_config = config_loader.load_config()

# Get enabled prediction runs
enabled_prediction_runs = config_loader.get_enabled_prediction_runs(full_config, args.prediction_runs)

# Get enabled analysis runs
enabled_analysis_runs = config_loader.get_enabled_analysis_runs(full_config, args.analysis_runs)

# Get a motif definition
motif_def = config_loader.get_motif_definition(full_config, "binding_pocket_hM1D")

# Get a motif template
template_path = config_loader.get_motif_template(full_config, "binding_pocket_hM1D")

# Get a motif template and model_idx
template_path, model_idx = config_loader.get_motif_template_and_model_idx(full_config, "binding_pocket_hM1D")
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
- Supports molecule-specific templates

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
- `get_molecule_specific_template(molecule_name, full_config, templates_dir)`: Gets the specific template for a molecule based on motif definitions
- `get_templates(config, args, full_config=None)`: Gets all template files from configuration or command-line arguments
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
```

---

## motif_alignment.py

### Purpose
Performs motif-specific alignment and RMSD calculation.

### Functionality
- Takes protein structures and motif definitions from the configuration
- Extracts the specified motif regions (e.g., binding pockets, active sites)
- Performs alignment on just those regions
- Calculates motif-specific RMSD values
- Creates specialized PyMOL session files that highlight the motif regions
- Supports both specific residue motifs and whole protein motifs

### Command-line Arguments
```
python motif_alignment.py [--motif MOTIF_ID] [--template TEMPLATE_FILE]
                         [--model-idx N] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_motif_definition(config, motif_id)`: Gets a motif definition by ID
- `find_cif_files(motif_def, config, quiet=False)`: Finds all CIF files for the molecules in the motif definition
- `perform_motif_alignment(motif_def, cif_files, template_file, model_idx, output_dir, quiet=False)`: Performs motif-specific alignment and RMSD calculation
- `main()`: Main function that orchestrates the process

### Output Files
- PyMOL session files (.pse) for each motif
- CSV files with motif-specific RMSD values

### Example Usage
```bash
# Perform motif-specific alignment with default settings
python motif_alignment.py --motif binding_pocket_hM1D

# Perform motif-specific alignment with a specific template
python motif_alignment.py --motif binding_pocket_hM1D --template custom_template.cif

# Perform motif-specific alignment with a specific model index
python motif_alignment.py --motif binding_pocket_hM1D --model-idx 0
```

---

## extract_motif_plddt.py

### Purpose
Extracts pLDDT values for specific motifs.

### Functionality
- Finds all JSON files with pLDDT values
- Extracts pLDDT values for the specified motif residues
- Calculates average pLDDT values for the motif regions
- Supports both specific residue motifs and whole protein motifs

### Command-line Arguments
```
python extract_motif_plddt.py [--motif MOTIF_ID] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_motif_definition(config, motif_id)`: Gets a motif definition by ID
- `find_json_files(motif_def, config, quiet=False)`: Finds all JSON files with pLDDT values for the molecules in the motif definition
- `extract_motif_plddt(motif_def, json_files, output_dir, quiet=False)`: Extracts pLDDT values for the specified motif residues
- `main()`: Main function that orchestrates the process

### Output Files
- CSV files with motif-specific pLDDT values

### Example Usage
```bash
# Extract motif-specific pLDDT values with default settings
python extract_motif_plddt.py --motif binding_pocket_hM1D
```

---

## plot_motif_rmsd.py

### Purpose
Generates heatmap visualizations of motif-specific RMSD values.

### Functionality
- Reads motif-specific RMSD values from CSV files
- Creates a heatmap with methods on the x-axis, ligand names on the y-axis, and RMSD values as cell values
- Saves the heatmap as a PNG file
- Saves the data to a CSV file for further analysis

### Command-line Arguments
```
python plot_motif_rmsd.py [--motif MOTIF_ID] [--quiet]
```

### Key Functions
- `parse_arguments()`: Parses command-line arguments
- `get_motif_definition(config, motif_id)`: Gets a motif definition by ID
- `find_rmsd_csv_files(motif_id, quiet=False)`: Finds all rmsd_values.csv files for the specified motif
- `read_rmsd_values(csv_file, quiet=False)`: Reads RMSD values from a CSV file
- `create_heatmap(df, output_file, config, quiet=False)`: Creates a heatmap visualization of RMSD values
- `main()`: Main function that orchestrates the process
