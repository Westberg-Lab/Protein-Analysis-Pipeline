#!/usr/bin/env python3
"""
Script to archive previous outputs and create fresh directories for new runs.

This script:
1. Creates a timestamped archive directory
2. Moves previous output directories and files to the archive
3. Creates fresh empty directories for a new run

Usage:
    python archive_and_clean.py [--no-archive]

Options:
    --no-archive    Delete previous outputs without archiving
"""

import os
import shutil
import argparse
from datetime import datetime
from pathlib import Path

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Archive previous outputs and create fresh directories.')
    parser.add_argument('--no-archive', action='store_true',
                        help='Delete previous outputs without archiving')
    return parser.parse_args()

def create_archive_directory():
    """Create a timestamped archive directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = Path(f"archive_{timestamp}")
    archive_dir.mkdir(exist_ok=True)
    print(f"Created archive directory: {archive_dir}")
    return archive_dir

def archive_directories(archive_dir, dirs_to_archive):
    """Move directories to the archive directory."""
    for dir_path in dirs_to_archive:
        dir_path = Path(dir_path)
        if dir_path.exists():
            # Create the destination directory in the archive
            dest_dir = archive_dir / dir_path
            dest_dir.parent.mkdir(parents=True, exist_ok=True)
            
            # Move the directory to the archive
            shutil.move(str(dir_path), str(dest_dir))
            print(f"Archived: {dir_path} -> {dest_dir}")
        else:
            print(f"Skipped (not found): {dir_path}")

def archive_files(archive_dir, files_to_archive):
    """Move files to the archive directory."""
    for file_path in files_to_archive:
        file_path = Path(file_path)
        if file_path.exists():
            # Move the file to the archive
            shutil.move(str(file_path), str(archive_dir / file_path.name))
            print(f"Archived: {file_path} -> {archive_dir / file_path.name}")
        else:
            print(f"Skipped (not found): {file_path}")

def delete_directories(dirs_to_delete):
    """Delete directories without archiving."""
    for dir_path in dirs_to_delete:
        dir_path = Path(dir_path)
        if dir_path.exists():
            shutil.rmtree(dir_path)
            print(f"Deleted: {dir_path}")
        else:
            print(f"Skipped (not found): {dir_path}")

def delete_files(files_to_delete):
    """Delete files without archiving."""
    for file_path in files_to_delete:
        file_path = Path(file_path)
        if file_path.exists():
            file_path.unlink()
            print(f"Deleted: {file_path}")
        else:
            print(f"Skipped (not found): {file_path}")

def create_fresh_directories(dirs_to_create):
    """Create fresh empty directories."""
    for dir_path in dirs_to_create:
        dir_path = Path(dir_path)
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {dir_path}")

def main():
    """Main function."""
    args = parse_arguments()
    
    # Define directories and files to archive/delete
    dirs_to_handle = [
        "CHAI_FASTA",
        "BOLTZ_YAML",
        "OUTPUT",
        "PSE_FILES",
        "plots"
    ]
    
    files_to_handle = [
        "rmsd_values.csv",
        "plddt_values.csv",
        "rmsd_heatmap.png",
        "plddt_heatmap.png"
    ]
    
    if args.no_archive:
        # Delete without archiving
        print("Deleting previous outputs without archiving...")
        delete_directories(dirs_to_handle)
        delete_files(files_to_handle)
    else:
        # Create archive directory
        archive_dir = create_archive_directory()
        
        # Archive directories and files
        print("Archiving previous outputs...")
        archive_directories(archive_dir, dirs_to_handle)
        archive_files(archive_dir, files_to_handle)
    
    # Create fresh directories
    print("\nCreating fresh directories for new run...")
    create_fresh_directories(dirs_to_handle)
    
    print("\nClean-up complete! Ready for a new run with fresh directories.")

if __name__ == "__main__":
    main()
