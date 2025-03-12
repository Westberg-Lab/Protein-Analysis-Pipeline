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

def is_file_empty(file_path):
    """Check if a file is empty (zero bytes)."""
    return file_path.stat().st_size == 0

def is_dir_empty(dir_path):
    """
    Check if a directory is empty or contains only empty files and directories.
    
    A directory is considered empty if:
    1. It has no files or subdirectories, or
    2. It only contains empty files and empty subdirectories
    """
    if not any(dir_path.iterdir()):
        # Directory has no files or subdirectories
        return True
    
    # Check each item in the directory
    for item in dir_path.iterdir():
        if item.is_file():
            if not is_file_empty(item):
                # Found a non-empty file
                return False
        elif item.is_dir():
            if not is_dir_empty(item):
                # Found a non-empty subdirectory
                return False
    
    # All files and subdirectories are empty
    return True

def archive_directories(archive_dir, dirs_to_archive):
    """Move non-empty directories to the archive directory."""
    for dir_path in dirs_to_archive:
        dir_path = Path(dir_path)
        if dir_path.exists():
            # Skip empty directories
            if is_dir_empty(dir_path):
                print(f"Skipped (empty): {dir_path}")
                continue
            
            # Create the destination directory in the archive
            dest_dir = archive_dir / dir_path
            dest_dir.parent.mkdir(parents=True, exist_ok=True)
            
            # Move the directory to the archive
            shutil.move(str(dir_path), str(dest_dir))
            print(f"Archived: {dir_path} -> {dest_dir}")
        else:
            print(f"Skipped (not found): {dir_path}")

def archive_files(archive_dir, files_to_archive):
    """Move non-empty files to the archive directory."""
    for file_path in files_to_archive:
        file_path = Path(file_path)
        if file_path.exists():
            # Skip empty files
            if is_file_empty(file_path):
                print(f"Skipped (empty): {file_path}")
                continue
            
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

def create_project_directories():
    """Create all necessary project directories."""
    dirs = ["CHAI_FASTA", "BOLTZ_YAML", "OUTPUT", "PSE_FILES", "plots", "csv"]
    for dir_path in dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {dir_path}")
    
    # Create subdirectories that are commonly needed
    Path("OUTPUT/CHAI").mkdir(parents=True, exist_ok=True)
    Path("OUTPUT/BOLTZ").mkdir(parents=True, exist_ok=True)

def main():
    """Main function."""
    args = parse_arguments()
    
    # Define directories and files to archive/delete
    dirs_to_handle = [
        "CHAI_FASTA",
        "BOLTZ_YAML",
        "OUTPUT",
        "PSE_FILES",
        "plots",
        "csv"
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
    create_project_directories()
    
    print("\nClean-up complete! Ready for a new run with fresh directories.")

if __name__ == "__main__":
    main()
