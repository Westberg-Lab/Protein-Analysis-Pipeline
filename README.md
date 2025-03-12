# CHAI and Boltz-1 Workflow Scripts

Scripts for generating and processing files with CHAI and Boltz-1.

## CHAI Workflow

### Quick Start

1. Generate FASTA files:
```bash
python generate_fasta_files.py
```
Creates FASTA files in CHAI_FASTA/ combining KORDshort variants with ligands.

2. Process with CHAI:
```bash
python run_chai_apptainer.py
```
Prompts:
- "Use MSA? (y/n)"
  - No: Uses msa_server=0
  - Yes: Goes to next prompt
- "Use MSA directory? (y/n)"
  - No: Uses msa_server=1
  - Yes: Uses msa_directory=MSAs + msa_server=1

### Directory Structure
```
.
├── CHAI_FASTA/           # Input FASTA files
├── MSAs/                 # MSA files (if using)
└── OUTPUT/CHAI/          # Results
    ├── (MoleculeNameHere)/       # Without MSA
    └── (MoleculeNameHere)_with_MSA/  # With MSA
```

## Boltz-1 Workflow

### Quick Start

1. Generate YAML files:
```bash
python generate_boltz_yaml.py
```
Creates YAML files in BOLTZ_YAML/ combining KORDshort variants with ligands.

2. Process with Boltz-1:
```bash
python run_boltz_apptainer.py
```
Prompts:
- "Use MSA server? (y/n)"
  - No: Runs without MSA
  - Yes: Adds --use_msa_server flag

### Directory Structure
```
.
├── BOLTZ_YAML/          # Input YAML files
└── OUTPUT/BOLTZ/        # Results
    └── (MoleculeNameHere)/  # Output for each molecule
```

### YAML File Format
```yaml
sequences:
  - protein:
      id: chain_A
      sequence: SEQUENCE_HERE
  - ligand:
      id: chain_B
      smiles: 'SMILES_HERE'
```
