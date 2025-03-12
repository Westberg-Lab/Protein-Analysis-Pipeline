import subprocess
from pathlib import Path

def get_msa_config():
    """Get MSA configuration through simple prompt."""
    use_msa = input("Do you want to use MSA server? (y/n): ").lower().strip() == 'y'
    return ["--use_msa_server"], use_msa  # Return both the config and the boolean

def run_apptainer_commands():
    """Run apptainer commands for each YAML file in BOLTZ_YAML subdirectories."""
    msa_config, using_msa = get_msa_config()
    
    # Create output directory
    output_base = Path('OUTPUT/BOLTZ')
    output_base.mkdir(parents=True, exist_ok=True)
    
    for folder in Path('BOLTZ_YAML').iterdir():
        if folder.is_dir():
            # Create subdirectory for each input folder with _with_MSA suffix if using MSA
            folder_name = f"{folder.name}_with_MSA" if using_msa else folder.name
            output_dir = output_base / folder_name
            
            output_dir.mkdir(exist_ok=True)
            
            for yaml_file in folder.glob('*.yaml'):
                # Get base name and add _with_MSA if using MSA
                base_name = yaml_file.stem  # Get filename without extension
                if using_msa:
                    base_name = f"{base_name}_with_MSA"
                
                # Check for typical output directories with boltz_results_ prefix
                expected_dirs = [
                    output_dir / f"boltz_results_{base_name}" / "predictions",
                    output_dir / f"boltz_results_{base_name}" / "processed"
                ]
                
                if all(d.exists() for d in expected_dirs):
                    print(f"Skipping {yaml_file.name} - output directories already exist")
                    continue
                
                cmd = [
                    "apptainer", "run", "--nv",
                    "/emcc/westberg/shared/containers/boltz.sif",
                    str(yaml_file),  # YAML file path as positional argument
                    f"--out_dir=OUTPUT/BOLTZ/{folder_name}"
                ] + msa_config
                
                print(f"Processing: {yaml_file}")
                subprocess.run(cmd)
                print(f"Completed: {yaml_file}\n")

if __name__ == "__main__":
    run_apptainer_commands()
    print("All YAML files processed!")
