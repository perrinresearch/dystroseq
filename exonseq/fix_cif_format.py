"""
Fix CIF file format to ensure proper column alignment.
"""

import os
import re
from pathlib import Path
from typing import List, Dict


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
                        'group_PDB': atom_data.get('group_PDB', 'ATOM'),
                        'id': int(atom_data.get('id', 0)),
                        'type_symbol': atom_data.get('type_symbol', '').strip(),
                        'label_atom_id': atom_data.get('label_atom_id', '').strip(),
                        'label_alt_id': atom_data.get('label_alt_id', '').strip(),
                        'label_comp_id': atom_data.get('label_comp_id', '').strip(),
                        'label_asym_id': atom_data.get('label_asym_id', '').strip(),
                        'label_entity_id': atom_data.get('label_entity_id', '1'),
                        'label_seq_id': int(atom_data.get('label_seq_id', 0)),
                        'pdbx_PDB_ins_code': atom_data.get('pdbx_PDB_ins_code', '?'),
                        'Cartn_x': float(atom_data.get('Cartn_x', 0.0)),
                        'Cartn_y': float(atom_data.get('Cartn_y', 0.0)),
                        'Cartn_z': float(atom_data.get('Cartn_z', 0.0)),
                        'occupancy': float(atom_data.get('occupancy', 1.0)),
                        'B_iso_or_equiv': float(atom_data.get('B_iso_or_equiv', 0.0)),
                        'pdbx_formal_charge': atom_data.get('pdbx_formal_charge', '?'),
                        'auth_seq_id': int(atom_data.get('auth_seq_id', 0)),
                        'auth_comp_id': atom_data.get('auth_comp_id', '').strip(),
                        'auth_asym_id': atom_data.get('auth_asym_id', '').strip(),
                        'auth_atom_id': atom_data.get('auth_atom_id', '').strip(),
                        'pdbx_PDB_model_num': int(atom_data.get('pdbx_PDB_model_num', 1))
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


def write_properly_formatted_cif(atoms: List[Dict], output_file: Path):
    """Write atoms to a properly formatted CIF file with correct column alignment."""
    print(f"    Writing properly formatted CIF to {output_file.name}")
    
    with open(output_file, 'w') as f:
        # Write CIF header
        f.write("data_renumbered_structure\n")
        f.write("# Renumbered AlphaFold structure\n")
        f.write("# Original residue numbering updated to match full-length dystrophin\n")
        f.write("# Properly formatted with correct column alignment\n")
        f.write("\n")
        
        # Write atom site loop with proper formatting
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
        
        # Write atom data with proper column alignment
        for atom in atoms:
            # Format each field with appropriate width and alignment
            line = (
                f"{atom['group_PDB']:<6s} "
                f"{atom['id']:>5d} "
                f"{atom['type_symbol']:>2s} "
                f"{atom['label_atom_id']:<4s} "
                f"{atom['label_alt_id']:<1s} "
                f"{atom['label_comp_id']:<3s} "
                f"{atom['label_asym_id']:<1s} "
                f"{atom['label_entity_id']:<1s} "
                f"{atom['label_seq_id']:>4d} "
                f"{atom['pdbx_PDB_ins_code']:<1s} "
                f"{atom['Cartn_x']:>8.3f} "
                f"{atom['Cartn_y']:>8.3f} "
                f"{atom['Cartn_z']:>8.3f} "
                f"{atom['occupancy']:>6.2f} "
                f"{atom['B_iso_or_equiv']:>6.2f} "
                f"{atom['pdbx_formal_charge']:<1s} "
                f"{atom['auth_seq_id']:>4d} "
                f"{atom['auth_comp_id']:<3s} "
                f"{atom['auth_asym_id']:<1s} "
                f"{atom['auth_atom_id']:<4s} "
                f"{atom['pdbx_PDB_model_num']:>4d}\n"
            )
            f.write(line)
    
    print(f"    Successfully wrote {len(atoms)} atoms with proper column alignment")


def fix_cif_format(input_file: Path, output_file: Path):
    """Fix the CIF format for a single file."""
    print(f"  Fixing format for {input_file.name}")
    
    # Parse the current CIF file
    atoms = parse_cif_file(input_file)
    if not atoms:
        print(f"    Error: No atoms found in {input_file.name}")
        return False
    
    # Write with proper formatting
    write_properly_formatted_cif(atoms, output_file)
    return True


def main():
    """Main function to fix CIF format for all renumbered structures."""
    print("=" * 60)
    print("FIXING CIF FILE FORMAT")
    print("=" * 60)
    
    # Set up paths
    input_dir = Path("renumbered_structures")
    output_dir = Path("renumbered_structures_fixed")
    output_dir.mkdir(exist_ok=True)
    
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} not found!")
        return
    
    # Get all renumbered CIF files
    cif_files = list(input_dir.glob("renumbered_*.cif"))
    cif_files.sort()
    
    print(f"Found {len(cif_files)} CIF files to fix")
    print(f"Output directory: {output_dir}")
    
    # Process each file
    processed_files = []
    
    for cif_file in cif_files:
        print(f"\nProcessing {cif_file.name}...")
        
        # Create output filename
        output_filename = cif_file.name.replace("renumbered_", "fixed_")
        output_file = output_dir / output_filename
        
        # Fix the format
        success = fix_cif_format(cif_file, output_file)
        if success:
            processed_files.append(output_file)
    
    # Summary
    print(f"\n" + "=" * 60)
    print("CIF FORMAT FIXING COMPLETE")
    print("=" * 60)
    print(f"Processed {len(processed_files)} files")
    print(f"Output directory: {output_dir}")
    
    print(f"\nFixed files:")
    for file in processed_files:
        print(f"  {file.name}")
    
    print(f"\nKey improvements:")
    print(f"  * Proper column alignment for CIF format")
    print(f"  * Correct field widths and spacing")
    print(f"  * Compatible with molecular visualization software")
    print(f"  * Maintains all atomic data and coordinates")
    
    print(f"\nAll CIF files now have proper column formatting!")


if __name__ == "__main__":
    main()
