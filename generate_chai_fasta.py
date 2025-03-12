import os
import json

def load_molecules(json_file="molecules.json"):
    """Load molecule definitions from JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
        return data["molecule_1"], data["molecule_2"]

def generate_fasta_files(output_base_dir="CHAI_FASTA"):
    """Generate FASTA files for combinations."""
    molecule_1, molecule_2 = load_molecules()
    
    for mol1 in molecule_1:
        if not mol1:
            continue
            
        mol1_type, mol1_name, mol1_seq = mol1
        mol_dir = os.path.join(output_base_dir, mol1_name)
        os.makedirs(mol_dir, exist_ok=True)
        
        # Create combined files for all valid entries
        for mol2 in molecule_2:
            if not mol2:
                # If this specific entry is empty, create single fasta file
                with open(f"{mol_dir}/{mol1_name}.fasta", "w") as f:
                    f.write(f">{mol1_type}|name={mol1_name}\n{mol1_seq}\n")
            else:
                mol2_type, mol2_name, mol2_seq = mol2
                with open(f"{mol_dir}/{mol1_name}_{mol2_name}.fasta", "w") as f:
                    f.write(f">{mol1_type}|name={mol1_name}\n{mol1_seq}\n")
                    f.write(f">{mol2_type}|name={mol2_name}\n{mol2_seq}\n")

if __name__ == "__main__":
    generate_fasta_files()
    print("FASTA files generated successfully!")
