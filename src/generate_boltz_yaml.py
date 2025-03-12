import os
import yaml
import json
import string

def load_molecules(json_file="molecules.json"):
    """Load molecule definitions from JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
        return data["molecule_1"], data["molecule_2"]

def generate_yaml_files(output_base_dir="BOLTZ_YAML", use_msa=False):
    """Generate YAML files for combinations following Boltz-1 schema."""
    molecule_1, molecule_2 = load_molecules()
    
    for mol1 in molecule_1:
        if not mol1:
            continue
            
        mol1_type, mol1_name, mol1_seq = mol1
        mol_dir = os.path.join(output_base_dir, mol1_name)
        os.makedirs(mol_dir, exist_ok=True)
        
        # Create YAML files for all valid entries
        for mol2 in molecule_2:
            ENTITY_LIST = []
            
            # Add first molecule (always first letter)
            if mol1_type == "protein":
                entity = {"protein": {"id": "A", "sequence": mol1_seq}}
                # Add MSA field
                entity["protein"]["msa"] = "empty"  # Use empty for explicit single sequence mode
            else:  # ligand
                entity = {"ligand": {"id": "A", "smiles": mol1_seq}}
            ENTITY_LIST.append(entity)
            
            if not mol2:
                # If mol2 is empty, create single molecule YAML file
                filename = f"{mol_dir}/{mol1_name}.yaml"
            else:
                mol2_type, mol2_name, mol2_seq = mol2
                # Add second molecule (always second letter)
                if mol2_type == "protein":
                    entity = {"protein": {"id": "B", "sequence": mol2_seq}}
                    # Add MSA field
                    entity["protein"]["msa"] = "empty"  # Use empty for explicit single sequence mode
                else:  # ligand
                    entity = {"ligand": {"id": "B", "smiles": mol2_seq}}
                ENTITY_LIST.append(entity)
                filename = f"{mol_dir}/{mol1_name}_{mol2_name}.yaml"
            
            yaml_data = {"sequences": ENTITY_LIST}
            
            # Write YAML file
            with open(filename, "w") as f:
                yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)

if __name__ == "__main__":
    # Ask user if they want to use MSA
    generate_yaml_files()
    print("YAML files generated successfully!")
