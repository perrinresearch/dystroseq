"""
Create clean PDB and CIF files by properly parsing the original AlphaFold structures.
"""

import os
import glob
import json
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import re


def extract_exon_range(dirname: str) -> Tuple[Optional[int], Optional[int]]:
    """Extract exon range from directory name."""
    parts = dirname.split('_')
    if len(parts) >= 4 and parts[1] == 'exons':
        try:
            start_exon = int(parts[2])
            end_exon = int(parts[3])
            return start_exon, end_exon
        except ValueError:
            pass
    
    # Fallback: try to extract numbers from the string
    numbers = re.findall(r'\d+', dirname)
    if len(numbers) >= 2:
        return int(numbers[0]), int(numbers[1])
    
    return None, None


def load_confidence_data(structure_dir: Path) -> List[float]:
    """Load confidence data for a structure."""
    confidence_files = list(structure_dir.glob("*_full_data_*.json"))
    if not confidence_files:
        return []
    
    # Load the first confidence file
    with open(confidence_files[0], 'r') as f:
        data = json.load(f)
    
    return data.get('atom_plddts', [])


def parse_cif_file_properly(cif_file: Path) -> List[Dict]:
    """Parse a CIF file properly and extract atomic coordinates."""
    atoms = []
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    in_atom_site = False
    atom_columns = {}
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            # Parse column definitions
            if 'label_atom_id' in line:
                atom_columns['atom_id'] = len(atom_columns)
            elif 'label_comp_id' in line:
                atom_columns['residue'] = len(atom_columns)
            elif 'label_seq_id' in line:
                atom_columns['seq_id'] = len(atom_columns)
            elif 'Cartn_x' in line:
                atom_columns['x'] = len(atom_columns)
            elif 'Cartn_y' in line:
                atom_columns['y'] = len(atom_columns)
            elif 'Cartn_z' in line:
                atom_columns['z'] = len(atom_columns)
            elif 'B_iso_or_equiv' in line:
                atom_columns['b_factor'] = len(atom_columns)
        
        elif line.startswith('ATOM') and in_atom_site:
            # Parse atom data - format: ATOM id type atom_id . residue chain entity seq_id ? x y z occupancy b_factor auth_seq_id auth_chain model
            fields = line.split()
            if len(fields) >= 15:
                try:
                    atom_data = {
                        'atom_id': fields[3],  # label_atom_id
                        'residue': fields[5],  # label_comp_id
                        'seq_id': int(fields[8]),  # label_seq_id
                        'x': float(fields[10]),  # Cartn_x
                        'y': float(fields[11]),  # Cartn_y
                        'z': float(fields[12]),  # Cartn_z
                        'b_factor': float(fields[14])  # B_iso_or_equiv
                    }
                    atoms.append(atom_data)
                except (ValueError, IndexError):
                    continue
        
        elif line.startswith('#') and in_atom_site:
            # End of atom site section
            break
    
    return atoms


def find_structures(base_dir: str) -> List[Dict]:
    """Find all structure directories and their data."""
    base_path = Path(base_dir)
    structures = []
    
    structure_dirs = [d for d in os.listdir(base_path) 
                     if os.path.isdir(base_path / d) and d.startswith('dystroseq_exons_')]
    
    for structure_dir in sorted(structure_dirs):
        print(f"Processing {structure_dir}...")
        
        # Extract exon range
        start_exon, end_exon = extract_exon_range(structure_dir)
        if start_exon is None or end_exon is None:
            print(f"Warning: Could not extract exon range from {structure_dir}")
            continue
        
        # Find CIF files
        cif_files = list((base_path / structure_dir).glob("*_model_0.cif"))
        if not cif_files:
            print(f"Warning: No model_0.cif files found in {structure_dir}")
            continue
        
        # Load confidence data
        confidence_data = load_confidence_data(base_path / structure_dir)
        
        structures.append({
            'name': structure_dir,
            'start_exon': start_exon,
            'end_exon': end_exon,
            'cif_file': cif_files[0],
            'confidence_data': confidence_data
        })
    
    return structures


def create_clean_pdb_file(atoms: List[Dict], confidence_scores: List[float], output_file: str):
    """Create a clean PDB file with proper formatting."""
    print(f"Creating clean PDB file: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("HEADER    COMBINED DYSTROPHIN STRUCTURE                   01-JAN-24   COMB\n")
        f.write("TITLE     COMBINED ALPHAFOLD STRUCTURES FOR DYSTROPHIN\n")
        f.write("REMARK   1 CONFIDENCE SCORES STORED IN B-FACTOR COLUMN\n")
        f.write("REMARK   2 GENERATED BY DYSTROSEQ STRUCTURE COMBINER\n")
        f.write("REMARK   3 CLEAN FORMAT FOR VMD COMPATIBILITY\n\n")
        
        # Write atoms with proper PDB format
        for i, atom in enumerate(atoms):
            # Get confidence score for this atom
            confidence = 50.0  # Default confidence
            if i < len(confidence_scores):
                confidence = confidence_scores[i]
            
            # Get element symbol from atom name
            element = atom['atom_id'][0] if atom['atom_id'] else 'C'
            
            # Format according to PDB specification
            # Columns: 1-6, 7-11, 13-16, 17, 18-20, 22, 23-26, 31-38, 39-46, 47-54, 55-60, 61-66, 77-78
            pdb_line = f"ATOM  {i+1:5d}  {atom['atom_id']:4s} {atom['residue']:3s} A{atom['seq_id']:4d}    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  {1.00:4.2f}{confidence:6.2f}           {element:2s}\n"
            f.write(pdb_line)
        
        # Write END
        f.write("END\n")
    
    print(f"Clean PDB file created: {output_file}")


def create_clean_cif_file(atoms: List[Dict], confidence_scores: List[float], output_file: str):
    """Create a clean CIF file with proper formatting."""
    print(f"Creating clean CIF file: {output_file}")
    
    with open(output_file, 'w') as f:
        # Write CIF header
        f.write("data_combined_dystrophin\n")
        f.write("# Combined AlphaFold structures for dystrophin\n")
        f.write("# Generated by DystroSeq structure combiner\n")
        f.write("# Confidence scores stored in B-factor column\n")
        f.write("# Clean format for VMD compatibility\n\n")
        
        # Write cell parameters
        f.write("_cell.length_a    200.0\n")
        f.write("_cell.length_b    200.0\n")
        f.write("_cell.length_c    200.0\n")
        f.write("_cell.angle_alpha 90.0\n")
        f.write("_cell.angle_beta  90.0\n")
        f.write("_cell.angle_gamma 90.0\n\n")
        
        # Write space group
        f.write("_space_group.name_H-M 'P 1'\n")
        f.write("_space_group.IT_number 1\n\n")
        
        # Write atom site loop
        f.write("loop_\n")
        f.write("_atom_site.group_PDB\n")
        f.write("_atom_site.id\n")
        f.write("_atom_site.type_symbol\n")
        f.write("_atom_site.label_atom_id\n")
        f.write("_atom_site.label_alt_id\n")
        f.write("_atom_site.label_comp_id\n")
        f.write("_atom_site.label_asym_id\n")
        f.write("_atom_site.label_entity_id\n")
        f.write("_atom_site.label_seq_id\n")
        f.write("_atom_site.pdbx_PDB_ins_code\n")
        f.write("_atom_site.Cartn_x\n")
        f.write("_atom_site.Cartn_y\n")
        f.write("_atom_site.Cartn_z\n")
        f.write("_atom_site.occupancy\n")
        f.write("_atom_site.B_iso_or_equiv\n")
        f.write("_atom_site.pdbx_formal_charge\n")
        f.write("_atom_site.auth_seq_id\n")
        f.write("_atom_site.auth_comp_id\n")
        f.write("_atom_site.auth_atom_id\n")
        f.write("_atom_site.auth_asym_id\n")
        f.write("_atom_site.pdbx_PDB_model_num\n")
        
        # Write atom data
        for i, atom in enumerate(atoms):
            # Get confidence score for this atom
            confidence = 50.0  # Default confidence
            if i < len(confidence_scores):
                confidence = confidence_scores[i]
            
            # Get element symbol from atom name
            element = atom['atom_id'][0] if atom['atom_id'] else 'C'
            
            # Write CIF atom line
            f.write(f"ATOM {i+1:5d} {element:2s} {atom['atom_id']:4s} . {atom['residue']:3s} A 1 {atom['seq_id']:4d} ? "
                   f"{atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f} 1.00 {confidence:6.2f} ? "
                   f"{atom['seq_id']:4d} {atom['residue']:3s} {atom['atom_id']:4s} A 1\n")
    
    print(f"Clean CIF file created: {output_file}")


def main():
    """Main function to create clean structure files."""
    base_dir = "../../../Data/dystrophin/folds_2025_10_22_19_55"
    
    print("=" * 60)
    print("CREATING CLEAN STRUCTURE FILES")
    print("=" * 60)
    
    # Find all structures
    print("Finding AlphaFold structures...")
    structures = find_structures(base_dir)
    
    if not structures:
        print("No structures found!")
        return
    
    print(f"Found {len(structures)} structures:")
    for struct in structures:
        print(f"  {struct['name']}: exons {struct['start_exon']}-{struct['end_exon']}")
    
    # Combine structures
    print("\nCombining structures...")
    combined_atoms = []
    combined_confidence = []
    
    # Sort structures by start exon
    sorted_structures = sorted(structures, key=lambda x: x['start_exon'])
    
    for i, struct in enumerate(sorted_structures):
        print(f"Processing {struct['name']}...")
        
        # Parse CIF file
        atoms = parse_cif_file_properly(struct['cif_file'])
        confidence = struct['confidence_data']
        
        # Adjust sequence IDs to avoid conflicts
        if combined_atoms:
            max_seq_id = max(atom['seq_id'] for atom in combined_atoms)
            for atom in atoms:
                atom['seq_id'] += max_seq_id
        
        combined_atoms.extend(atoms)
        combined_confidence.extend(confidence)
        
        print(f"  Added {len(atoms)} atoms")
    
    print(f"\nCombined structure has {len(combined_atoms)} atoms")
    
    # Create clean files
    print("\nCreating clean structure files...")
    create_clean_pdb_file(combined_atoms, combined_confidence, "dystrophin_combined_clean.pdb")
    create_clean_cif_file(combined_atoms, combined_confidence, "dystrophin_combined_clean.cif")
    
    print(f"\nClean structure files created:")
    print(f"  - dystrophin_combined_clean.pdb (VMD-compatible PDB format)")
    print(f"  - dystrophin_combined_clean.cif (Proper CIF format)")
    print(f"  - Total atoms: {len(combined_atoms)}")
    print(f"  - Confidence scores: {len(combined_confidence)}")
    print(f"  - These files should now load correctly in VMD and other molecular viewers.")


if __name__ == "__main__":
    main()
