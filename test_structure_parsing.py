"""
Test the structure parsing to identify the issue.
"""

import os
import json
from pathlib import Path


def test_single_structure():
    """Test parsing a single structure to identify the issue."""
    
    base_dir = Path("../../../Data/dystrophin/folds_2025_10_22_19_55")
    test_dir = base_dir / "dystroseq_exons_1_9_protein"
    
    print("Testing single structure parsing...")
    print(f"Directory: {test_dir}")
    
    # Check if directory exists
    if not test_dir.exists():
        print(f"Directory does not exist: {test_dir}")
        return
    
    # List files in directory
    print("\nFiles in directory:")
    for file in test_dir.iterdir():
        print(f"  {file.name}")
    
    # Find CIF file
    cif_files = list(test_dir.glob("*_model_0.cif"))
    if not cif_files:
        print("No CIF files found!")
        return
    
    cif_file = cif_files[0]
    print(f"\nUsing CIF file: {cif_file}")
    
    # Parse CIF file
    print("\nParsing CIF file...")
    atoms = []
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    in_atom_site = False
    atom_count = 0
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            continue
        
        elif line.startswith('ATOM') and in_atom_site:
            atom_count += 1
            if atom_count <= 5:  # Show first 5 atoms
                print(f"  Atom line {atom_count}: {line}")
                
                # Parse the line
                fields = line.split()
                if len(fields) >= 15:
                    try:
                        atom_data = {
                            'atom_id': fields[3],
                            'residue': fields[5],
                            'seq_id': int(fields[8]),
                            'x': float(fields[10]),
                            'y': float(fields[11]),
                            'z': float(fields[12]),
                            'b_factor': float(fields[14])
                        }
                        atoms.append(atom_data)
                        print(f"    Parsed: {atom_data}")
                    except (ValueError, IndexError) as e:
                        print(f"    Error parsing: {e}")
        
        elif line.startswith('#') and in_atom_site:
            break
    
    print(f"\nTotal atoms parsed: {len(atoms)}")
    
    # Check confidence data
    print("\nChecking confidence data...")
    confidence_files = list(test_dir.glob("*_full_data_*.json"))
    if confidence_files:
        print(f"Found confidence file: {confidence_files[0]}")
        with open(confidence_files[0], 'r') as f:
            data = json.load(f)
        
        atom_plddts = data.get('atom_plddts', [])
        print(f"Confidence data length: {len(atom_plddts)}")
        if atom_plddts:
            print(f"First 5 confidence scores: {atom_plddts[:5]}")
    else:
        print("No confidence files found!")


def main():
    """Main test function."""
    test_single_structure()


if __name__ == "__main__":
    main()
