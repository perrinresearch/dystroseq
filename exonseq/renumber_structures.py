"""
Renumber AlphaFold structure files according to the residue numbering mapping.
"""

import os
import re
from pathlib import Path
from typing import Dict, List, Tuple


def parse_residue_numbering(filename: str) -> Dict[str, Tuple[int, int]]:
    """Parse the residue numbering file to get start/end positions for each structure."""
    print(f"Loading residue numbering from {filename}")
    
    numbering_map = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for line in lines[1:]:  # Skip header
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                filename = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                numbering_map[filename] = (start, end)
                print(f"  {filename}: residues {start}-{end}")
    
    return numbering_map


def parse_cif_file(cif_file: Path) -> List[Dict]:
    """Parse a CIF file and extract atomic coordinates and residue info."""
    print(f"  Parsing {cif_file.name}...")
    
    atoms = []
    
    with open(cif_file, 'r') as f:
        lines = f.readlines()
    
    in_atom_site = False
    atom_site_labels = []
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('loop_'):
            in_atom_site = False
            atom_site_labels = []
        
        if line.startswith('_atom_site.'):
            in_atom_site = True
            atom_site_labels.append(line.split('.')[1])
            continue
        
        if in_atom_site and not line.startswith(('_', '#', 'data_')) and line:
            # This is an atom data line
            fields = line.split()
            if len(fields) >= len(atom_site_labels):
                atom_data = {}
                for i, label in enumerate(atom_site_labels):
                    atom_data[label] = fields[i]
                
                try:
                    # Extract relevant fields
                    atom = {
                        'atom_id': atom_data.get('label_atom_id', '').strip(),
                        'type_symbol': atom_data.get('type_symbol', '').strip(),
                        'residue': atom_data.get('label_comp_id', '').strip(),
                        'chain': atom_data.get('label_asym_id', '').strip(),
                        'seq_id': int(atom_data.get('label_seq_id', 0)),
                        'x': float(atom_data.get('Cartn_x', 0.0)),
                        'y': float(atom_data.get('Cartn_y', 0.0)),
                        'z': float(atom_data.get('Cartn_z', 0.0)),
                        'occupancy': float(atom_data.get('occupancy', 1.0)),
                        'b_factor': float(atom_data.get('B_iso_or_equiv', 0.0)),
                        'element': atom_data.get('type_symbol', '').strip()
                    }
                    atoms.append(atom)
                except (ValueError, IndexError):
                    continue
        elif line.startswith('#') and in_atom_site:
            # End of atom site section
            in_atom_site = False
            atom_site_labels = []
    
    print(f"    Parsed {len(atoms)} atoms")
    return atoms


def renumber_cif_atoms(atoms: List[Dict], start_residue: int) -> List[Dict]:
    """Renumber atoms to start from the specified residue number."""
    print(f"    Renumbering atoms to start from residue {start_residue}")
    
    # Find the minimum original sequence ID
    min_original_id = min(atom['seq_id'] for atom in atoms)
    
    # Renumber all atoms
    for atom in atoms:
        # Calculate offset from minimum original ID
        offset = atom['seq_id'] - min_original_id
        # Set new sequence ID
        atom['seq_id'] = start_residue + offset
    
    return atoms


def write_renumbered_cif(atoms: List[Dict], output_file: Path):
    """Write renumbered atoms to a new CIF file."""
    print(f"    Writing renumbered structure to {output_file.name}")
    
    with open(output_file, 'w') as f:
        # Write CIF header
        f.write("data_renumbered_structure\n")
        f.write("# Renumbered AlphaFold structure\n")
        f.write("# Original residue numbering updated to match full-length dystrophin\n")
        f.write("\n")
        
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
        f.write("_atom_site.auth_asym_id\n")
        f.write("_atom_site.auth_atom_id\n")
        f.write("_atom_site.pdbx_PDB_model_num\n")
        
        # Write atom data
        for i, atom in enumerate(atoms, 1):
            f.write(f"ATOM {i} {atom['type_symbol']} {atom['atom_id']} . {atom['residue']} {atom['chain']} 1 {atom['seq_id']} ? {atom['x']:.3f} {atom['y']:.3f} {atom['z']:.3f} {atom['occupancy']:.2f} {atom['b_factor']:.2f} ? {atom['seq_id']} {atom['residue']} {atom['chain']} {atom['atom_id']} 1\n")
    
    print(f"    Successfully wrote {len(atoms)} atoms to {output_file.name}")


def process_structure_file(cif_file: Path, output_dir: Path, numbering_map: Dict[str, Tuple[int, int]]):
    """Process a single structure file and renumber it."""
    filename = cif_file.name
    
    if filename not in numbering_map:
        print(f"  Warning: No numbering information found for {filename}")
        return
    
    start_residue, end_residue = numbering_map[filename]
    print(f"  Processing {filename} -> residues {start_residue}-{end_residue}")
    
    # Parse the CIF file
    atoms = parse_cif_file(cif_file)
    if not atoms:
        print(f"    Error: No atoms found in {filename}")
        return
    
    # Renumber atoms
    renumbered_atoms = renumber_cif_atoms(atoms, start_residue)
    
    # Create output filename
    output_filename = f"renumbered_{filename}"
    output_file = output_dir / output_filename
    
    # Write renumbered structure
    write_renumbered_cif(renumbered_atoms, output_file)
    
    # Verify the renumbering
    min_seq_id = min(atom['seq_id'] for atom in renumbered_atoms)
    max_seq_id = max(atom['seq_id'] for atom in renumbered_atoms)
    print(f"    Renumbered range: {min_seq_id}-{max_seq_id}")
    
    return output_file


def main():
    """Main function to renumber all structure files."""
    print("=" * 60)
    print("RENUMBERING ALPHAFOLD STRUCTURES")
    print("=" * 60)
    
    # Load residue numbering
    numbering_map = parse_residue_numbering("dystrophin_residue_numbering.txt")
    
    # Set up paths
    base_dir = Path("../../../Data/dystrophin/folds_2025_10_22_19_55/")
    output_dir = Path("renumbered_structures")
    output_dir.mkdir(exist_ok=True)
    
    print(f"\nOutput directory: {output_dir}")
    
    # Process each structure directory
    processed_files = []
    
    for structure_dir in base_dir.iterdir():
        if structure_dir.is_dir() and structure_dir.name.startswith('dystroseq_exons_'):
            print(f"\nProcessing {structure_dir.name}...")
            
            # Find the model_0.cif file
            cif_file = structure_dir / f"fold_{structure_dir.name}_model_0.cif"
            
            if cif_file.exists():
                output_file = process_structure_file(cif_file, output_dir, numbering_map)
                if output_file:
                    processed_files.append(output_file)
            else:
                print(f"  Warning: {cif_file.name} not found")
    
    # Summary
    print(f"\n" + "=" * 60)
    print("RENUMBERING COMPLETE")
    print("=" * 60)
    print(f"Processed {len(processed_files)} structure files")
    print(f"Output directory: {output_dir}")
    
    print(f"\nRenumbered files:")
    for file in processed_files:
        print(f"  {file.name}")
    
    print(f"\nResidue numbering summary:")
    for filename, (start, end) in numbering_map.items():
        print(f"  {filename}: {start}-{end} ({end-start+1} residues)")
    
    print(f"\nAll structures have been renumbered to match the full-length dystrophin sequence!")
    print(f"Files are ready for structural alignment and combination.")


if __name__ == "__main__":
    main()
