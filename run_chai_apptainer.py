import subprocess
from pathlib import Path

def get_msa_config():
    """Get MSA configuration through simple prompts."""
    use_msa = input("Do you want to use MSA? (y/n): ").lower().strip() == 'y'
    
    if not use_msa:
        return [], False
    
    use_msa_dir = input("Do you want to use MSA directory? (y/n): ").lower().strip() == 'y'
    
    if use_msa_dir:
        return ["msa_directory=CHAI_MSAs", "msa_server=1"], True
    
    return ["msa_server=1"], True

def run_apptainer_commands():
    """Run apptainer commands for each FASTA file in CHAI_FASTA subdirectories."""
    msa_configs, using_msa = get_msa_config()
    print(f"MSA config is: {msa_configs}")
    
    # Create output directory
    output_base = Path('OUTPUT/CHAI')
    output_base.mkdir(parents=True, exist_ok=True)
    
    for folder in Path('CHAI_FASTA').iterdir():
        if folder.is_dir():
            # Create subdirectory for each input folder with _with_MSA suffix if using MSA
            folder_name = f"{folder.name}_with_MSA" if using_msa else folder.name
            output_dir = output_base / folder_name
            
            output_dir.mkdir(exist_ok=True)
            
            for fasta in folder.glob('*.fasta'):
                # Check if output files already exist for this FASTA
                base_name = fasta.stem  # Get filename without extension
                output_files_exist = False
                
                # Check for typical output files
                expected_files = [
                    output_dir / base_name / "outs.json",
                    output_dir / base_name / "pred.model_idx_0.cif"
                ]
                
                if all(f.exists() for f in expected_files):
                    print(f"Skipping {fasta.name} - output files already exist")
                    continue
                
                cmd = [
                    "apptainer", "run", "--nv",
                    "/emcc/westberg/shared/containers/chai.sif",
                    f"input_paths={fasta}",
                    f"outdir=OUTPUT/CHAI/{folder_name}",
                ] + msa_configs

                print(f"Running subprocess terminal command:\n{cmd}")
                
                print(f"Processing: {fasta}")
                subprocess.run(cmd)
                print(f"Completed: {fasta}\n")
    

if __name__ == "__main__":
    run_apptainer_commands()
    print("All FASTA files processed!")
