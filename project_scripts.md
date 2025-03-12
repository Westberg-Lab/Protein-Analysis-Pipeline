# Protein Prediction Project Scripts

## Pipeline Management Scripts

### `run_pipeline.py`
Master script to run the entire protein prediction pipeline from start to finish.
- Archives previous outputs using `archive_and_clean.py`
- Runs each script in the pipeline in the correct order
- Handles errors and provides status updates
- Supports skipping specific steps with command-line options
- Logs progress with timestamps

### `archive_and_clean.py`
Script to archive previous outputs and create fresh directories for new runs.
- Creates a timestamped archive directory
- Moves previous output directories and files to the archive
- Optionally deletes outputs without archiving
- Creates fresh empty directories for a new run

## Data Generation Scripts

### `generate_chai_fasta.py`
Generates FASTA files for protein sequences and their combinations based on definitions in `molecules.json`. These files are used as input for the CHAI protein structure prediction tool.
- Creates a directory structure in `CHAI_FASTA/` with subdirectories for each molecule
- Generates individual FASTA files for single molecules and combination files for molecule pairs
- Each FASTA file contains proper headers with molecule type and name information

### `generate_boltz_yaml.py`
Generates YAML configuration files for the Boltz-1 protein structure prediction tool based on definitions in `molecules.json`.
- Creates a directory structure in `BOLTZ_YAML/` with subdirectories for each molecule
- Generates YAML files following the Boltz-1 schema for both single molecules and combinations
- Configures protein and ligand entities with appropriate IDs and sequences/SMILES

## Prediction Execution Scripts

### `run_chai_apptainer.py`
Executes the CHAI protein structure prediction tool using Apptainer containers for each generated FASTA file.
- Prompts the user for MSA (Multiple Sequence Alignment) configuration options
- Creates an output directory structure in `OUTPUT/CHAI/`
- Processes each FASTA file in the `CHAI_FASTA/` directory
- Skips files that have already been processed
- Runs the CHAI container with appropriate parameters

### `run_boltz_apptainer.py`
Executes the Boltz protein structure prediction tool using Apptainer containers for each generated YAML file.
- Prompts the user for MSA server usage
- Creates an output directory structure in `OUTPUT/BOLTZ/`
- Processes each YAML file in the `BOLTZ_YAML/` directory
- Skips files where output directories already exist
- Runs the Boltz container with appropriate parameters

## Analysis and Visualization Scripts

### `combine_cif_files.py`
Creates PyMOL session (.pse) files that combine and align protein structures from different prediction methods.
- Finds unique directory names that exist in both CHAI and BOLTZ outputs
- For each unique name, loads the corresponding CIF files from CHAI and BOLTZ (with and without MSA)
- Aligns all structures to a template protein
- Calculates RMSD (Root Mean Square Deviation) values for each alignment
- Saves PyMOL sessions with colored visualizations
- Exports RMSD values to CSV files for further analysis

### `plot_rmsd_heatmap.py`
Generates heatmap visualizations of RMSD values to compare the quality of different prediction methods.
- Reads RMSD values from CSV files in the `PSE_FILES/` subdirectories
- Creates heatmaps with methods on the x-axis and ligands on the y-axis
- Uses a color gradient to visualize RMSD values (green for good alignments, red for poor alignments)
- Saves heatmaps to the `plots/` directory
- Supports filtering by reference protein

### `plot_plddt_heatmap.py`
Generates heatmap visualizations of pLDDT (predicted Local Distance Difference Test) values to assess prediction confidence.
- Recursively searches through `OUTPUT/CHAI/` and `OUTPUT/BOLTZ/` directories to find JSON files with pLDDT values
- Extracts complex_plddt values from both CHAI and BOLTZ outputs
- Creates a heatmap with methods on the x-axis and ligands on the y-axis
- Uses a color gradient to visualize pLDDT values (green for high confidence, red for low confidence)
- Saves the heatmap and exports the data to a CSV file
