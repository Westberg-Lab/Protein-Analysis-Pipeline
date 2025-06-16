# Protein Analysis Pipeline

A simple pipeline for protein structure prediction and analysis using CHAI and BOLTZ tools.

## What It Does

This pipeline:
1. Predicts protein structures using CHAI and BOLTZ tools
2. Analyzes the predictions (whole proteins and specific motifs)
3. Generates visualizations (PyMOL sessions and heatmaps)

## How to Use It

### Basic Usage

```bash
# Run the entire pipeline with default settings
python run_pipeline.py
```

### Common Options

```bash
# Run with a specific configuration file
python run_pipeline.py --config custom_config.json

# Skip archiving previous outputs
python run_pipeline.py --no-archive

# Skip specific steps
python run_pipeline.py --skip-step chai-fasta --skip-step boltz-yaml

# Resume from a previous run
python run_pipeline.py --resume
```

## Configuration

The pipeline is configured using `pipeline_config.json`, which defines:

1. **Prediction Runs**: Different prediction configurations (with/without MSA)
2. **Analysis Runs**: What analyses to perform on the predictions
3. **Motifs**: Specific regions of proteins to analyze (e.g., binding pockets)

Here's a simple example of the configuration file with all parameters:
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
          "residues": [85, 86, 89, 90, 93, 137, 169, 172, 173, 176, 177, 231, 234, 238, 235, 257, 260, 261],
          "chain": "A",
          "color": "magenta",
          "template": "templates/hMD_template_8E9X.cif"
        }
      ]
    }
  },
  "prediction_runs": [
    {
      "id": "standard",
      "description": "Standard prediction without MSA",
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
      "description": "Prediction with MSA enabled",
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
      "id": "binding_pocket_analysis",
      "description": "Binding pocket motif analysis",
      "enabled": true,
      "source_predictions": ["standard", "with_msa"],
      "analysis_type": "motif",
      "motif_id": "binding_pocket_hM1D",
      "metrics": ["rmsd", "plddt"]
    }
  ]
}
```

Molecules are defined in `molecules.json`. Here's a simple example:

```json
{
  "molecule_1": [
    ["protein", "KORDshort", "SEQUENCE..."],
    ["protein", "KORDshort_WT", "SEQUENCE..."]
  ],
  "molecule_2": [
    ["ligand", "SalA", "SMILES..."],
    ["ligand", "SalB", "SMILES..."],
    [],
    ["protein", "Arodyn1-6", "SEQUENCE..."],
    ["protein", "Arodyn1-7", "SEQUENCE..."]
  ]
}
```

Each entry follows the format: `[type, name, sequence]`

## Output Files

- **PyMOL Sessions**: `PSE_FILES/` directory
- **Heatmaps**: `plots/` directory
- **CSV Data**: `csv/` directory

## Requirements

- Python 3.6+
- PyMOL
- Apptainer (formerly Singularity)
- CHAI and BOLTZ containers
- Pandas, Matplotlib, Seaborn
