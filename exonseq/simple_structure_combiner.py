"""
Simplified structure combiner for AlphaFold CIF files.
Creates a combined structure with confidence scores in B-factor column.
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


def parse_cif_file(cif_file: Path) -> List[Dict]:
    """Parse a CIF file and extract atomic coordinates."""
    atoms = []
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    in_atom_site = False
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            continue
        
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


def combine_structures_sequential(structures: List[Dict]) -> Tuple[List[Dict], List[float]]:
    """Combine structures sequentially, handling overlaps."""
    print("Combining structures sequentially...")
    
    if not structures:
        return [], []
    
    # Sort by start exon
    sorted_structures = sorted(structures, key=lambda x: x['start_exon'])
    
    combined_atoms = []
    combined_confidence = []
    
    # Start with the first structure
    first_struct = sorted_structures[0]
    first_atoms = parse_cif_file(first_struct['cif_file'])
    combined_atoms.extend(first_atoms)
    combined_confidence.extend(first_struct['confidence_data'])
    
    print(f"Starting with {first_struct['name']} (exons {first_struct['start_exon']}-{first_struct['end_exon']})")
    
    # Add subsequent structures
    for i in range(1, len(sorted_structures)):
        current_struct = sorted_structures[i]
        prev_struct = sorted_structures[i-1]
        
        print(f"Adding {current_struct['name']} (exons {current_struct['start_exon']}-{current_struct['end_exon']})")
        
        # Check for overlap
        overlap_start = max(prev_struct['start_exon'], current_struct['start_exon'])
        overlap_end = min(prev_struct['end_exon'], current_struct['end_exon'])
        
        if overlap_start <= overlap_end:
            print(f"  Overlap detected in exons {overlap_start}-{overlap_end}")
            # For now, we'll add all atoms and let the user handle overlaps manually
            # In a more sophisticated approach, you'd align and merge overlapping regions
        
        # Add atoms from current structure
        current_atoms = parse_cif_file(current_struct['cif_file'])
        current_confidence = current_struct['confidence_data']
        
        # Adjust sequence IDs to avoid conflicts
        max_seq_id = max(atom['seq_id'] for atom in combined_atoms) if combined_atoms else 0
        
        for atom in current_atoms:
            # Adjust sequence ID to continue from previous structure
            atom['seq_id'] += max_seq_id
            combined_atoms.append(atom)
        
        # Add confidence data
        combined_confidence.extend(current_confidence)
    
    print(f"Combined structure has {len(combined_atoms)} atoms")
    return combined_atoms, combined_confidence


def save_combined_structure(atoms: List[Dict], confidence_scores: List[float], output_file: str):
    """Save the combined structure to a CIF file."""
    print(f"Saving combined structure to {output_file}")
    
    with open(output_file, 'w') as f:
        # Write CIF header
        f.write("data_combined_dystrophin\n")
        f.write("# Combined AlphaFold structures for dystrophin\n")
        f.write("# Generated by DystroSeq structure combiner\n")
        f.write("# Confidence scores stored in B-factor column\n\n")
        
        # Write cell parameters (placeholder)
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
            
            f.write(f"ATOM {i+1:5d} {element:2s} {atom['atom_id']:4s} . {atom['residue']:3s} A 1 {atom['seq_id']:4d} ? "
                   f"{atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f} 1.00 {confidence:6.2f} ? "
                   f"{atom['seq_id']:4d} {atom['residue']:3s} {atom['atom_id']:4s} A 1\n")


def create_pdb_file(atoms: List[Dict], confidence_scores: List[float], output_file: str):
    """Create a PDB file from the combined structure."""
    print(f"Creating PDB file: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("HEADER    COMBINED DYSTROPHIN STRUCTURE                   01-JAN-24   COMB\n")
        f.write("TITLE     COMBINED ALPHAFOLD STRUCTURES FOR DYSTROPHIN\n")
        f.write("REMARK   1 CONFIDENCE SCORES STORED IN B-FACTOR COLUMN\n")
        f.write("REMARK   2 GENERATED BY DYSTROSEQ STRUCTURE COMBINER\n\n")
        
        # Write atom data
        for i, atom in enumerate(atoms):
            # Get confidence score for this atom
            confidence = 50.0  # Default confidence
            if i < len(confidence_scores):
                confidence = confidence_scores[i]
            
            # Get element symbol from atom name
            element = atom['atom_id'][0] if atom['atom_id'] else 'C'
            
            # Format PDB line
            f.write(f"ATOM  {i+1:5d}  {atom['atom_id']:4s} {atom['residue']:3s} A{atom['seq_id']:4d}    "
                   f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00{confidence:6.2f}           {element:2s}\n")
        
        f.write("END\n")


def main():
    """Main function to combine AlphaFold structures."""
    base_dir = "../../../Data/dystrophin/folds_2025_10_22_19_55"
    
    print("=" * 60)
    print("DYSTROSEQ STRUCTURE COMBINER")
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
    combined_atoms, combined_confidence = combine_structures_sequential(structures)
    
    if not combined_atoms:
        print("Failed to combine structures!")
        return
    
    # Save combined structure
    print(f"\nSaving combined structure...")
    save_combined_structure(combined_atoms, combined_confidence, "combined_dystrophin.cif")
    create_pdb_file(combined_atoms, combined_confidence, "combined_dystrophin.pdb")
    
    print(f"\nStructure combination complete!")
    print(f"  Total atoms: {len(combined_atoms)}")
    print(f"  Confidence scores: {len(combined_confidence)}")
    print(f"  Output files: combined_dystrophin.cif, combined_dystrophin.pdb")
    print(f"  Confidence scores are stored in the B-factor column")


if __name__ == "__main__":
    main()
