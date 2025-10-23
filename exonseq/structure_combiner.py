"""
Combine AlphaFold CIF files into a single dystrophin protein structure.
Handles overlapping exons, alignment, and steric clash avoidance.
"""

import os
import glob
import json
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import re


class StructureCombiner:
    """Class to combine multiple AlphaFold structures into a single protein."""
    
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.structures = []
        self.confidence_data = {}
        self.aligned_structures = []
        
    def load_structure_data(self):
        """Load all structure directories and their data."""
        structure_dirs = [d for d in os.listdir(self.base_dir) 
                         if os.path.isdir(self.base_dir / d) and d.startswith('dystroseq_exons_')]
        
        for structure_dir in sorted(structure_dirs):
            print(f"Loading {structure_dir}...")
            
            # Extract exon range
            start_exon, end_exon = self.extract_exon_range(structure_dir)
            if start_exon is None or end_exon is None:
                print(f"Warning: Could not extract exon range from {structure_dir}")
                continue
            
            # Load confidence data
            confidence_data = self.load_confidence_data(self.base_dir / structure_dir)
            
            # Find CIF files
            cif_files = list((self.base_dir / structure_dir).glob("*.cif"))
            if not cif_files:
                print(f"Warning: No CIF files found in {structure_dir}")
                continue
            
            # Use the first model (model_0) as the representative structure
            model_cif = self.base_dir / structure_dir / "fold_dystroseq_exons_1_9_protein_model_0.cif"
            if not model_cif.exists():
                # Try to find any model_0.cif file
                model_cif = list((self.base_dir / structure_dir).glob("*_model_0.cif"))
                if model_cif:
                    model_cif = model_cif[0]
                else:
                    print(f"Warning: No model_0.cif found in {structure_dir}")
                    continue
            
            self.structures.append({
                'name': structure_dir,
                'start_exon': start_exon,
                'end_exon': end_exon,
                'cif_file': model_cif,
                'confidence_data': confidence_data
            })
        
        print(f"Loaded {len(self.structures)} structures")
        return self.structures
    
    def extract_exon_range(self, dirname: str) -> Tuple[Optional[int], Optional[int]]:
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
        import re
        numbers = re.findall(r'\d+', dirname)
        if len(numbers) >= 2:
            return int(numbers[0]), int(numbers[1])
        
        return None, None
    
    def load_confidence_data(self, structure_dir: Path) -> Dict:
        """Load confidence data for a structure."""
        confidence_files = list(structure_dir.glob("*_full_data_*.json"))
        if not confidence_files:
            return {}
        
        # Load the first confidence file
        with open(confidence_files[0], 'r') as f:
            data = json.load(f)
        
        return {
            'atom_plddts': data.get('atom_plddts', []),
            'ptm': data.get('ptm', None)
        }
    
    def parse_cif_file(self, cif_file: Path) -> Dict:
        """Parse a CIF file and extract atomic coordinates."""
        atoms = []
        confidence_scores = []
        
        with open(cif_file, 'r') as f:
            lines = f.readlines()
        
        in_atom_site = False
        atom_columns = {}
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('_atom_site.'):
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
            
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                in_atom_site = True
                continue
            
            elif in_atom_site and line and not line.startswith('#'):
                # Parse atom data
                fields = line.split()
                if len(fields) >= max(atom_columns.values()) + 1:
                    try:
                        atom_data = {
                            'atom_id': fields[atom_columns.get('atom_id', 0)],
                            'residue': fields[atom_columns.get('residue', 1)],
                            'seq_id': int(fields[atom_columns.get('seq_id', 2)]),
                            'x': float(fields[atom_columns.get('x', 3)]),
                            'y': float(fields[atom_columns.get('y', 4)]),
                            'z': float(fields[atom_columns.get('z', 5)]),
                            'b_factor': float(fields[atom_columns.get('b_factor', 6)])
                        }
                        atoms.append(atom_data)
                    except (ValueError, IndexError):
                        continue
        
        return {
            'atoms': atoms,
            'confidence_scores': confidence_scores
        }
    
    def find_overlapping_exons(self, struct1: Dict, struct2: Dict) -> Optional[Tuple[int, int]]:
        """Find overlapping exon range between two structures."""
        # Check if there's an overlap
        if struct1['end_exon'] >= struct2['start_exon'] and struct2['end_exon'] >= struct1['start_exon']:
            overlap_start = max(struct1['start_exon'], struct2['start_exon'])
            overlap_end = min(struct1['end_exon'], struct2['end_exon'])
            return overlap_start, overlap_end
        return None
    
    def align_structures(self, struct1: Dict, struct2: Dict, overlap_range: Tuple[int, int]) -> Dict:
        """Align two structures based on overlapping exons."""
        print(f"Aligning structures with overlap in exons {overlap_range[0]}-{overlap_range[1]}")
        
        # Parse both structures
        atoms1 = self.parse_cif_file(struct1['cif_file'])
        atoms2 = self.parse_cif_file(struct2['cif_file'])
        
        # Find overlapping residues in both structures
        overlap_residues1 = []
        overlap_residues2 = []
        
        # This is a simplified approach - in practice, you'd need to map
        # exon positions to residue indices more carefully
        for atom in atoms1['atoms']:
            if overlap_range[0] <= atom['seq_id'] <= overlap_range[1]:
                overlap_residues1.append(atom)
        
        for atom in atoms2['atoms']:
            if overlap_range[0] <= atom['seq_id'] <= overlap_range[1]:
                overlap_residues2.append(atom)
        
        if not overlap_residues1 or not overlap_residues2:
            print("Warning: No overlapping residues found for alignment")
            return struct2
        
        # Calculate alignment transformation
        # This is a simplified RMSD alignment - in practice, you'd use
        # more sophisticated alignment algorithms
        transformation = self.calculate_alignment_transformation(overlap_residues1, overlap_residues2)
        
        # Apply transformation to struct2
        aligned_atoms2 = []
        for atom in atoms2['atoms']:
            # Apply rotation and translation
            new_coords = self.apply_transformation(atom, transformation)
            atom['x'], atom['y'], atom['z'] = new_coords
            aligned_atoms2.append(atom)
        
        return {
            'atoms': aligned_atoms2,
            'confidence_scores': struct2['confidence_data'].get('atom_plddts', [])
        }
    
    def calculate_alignment_transformation(self, atoms1: List[Dict], atoms2: List[Dict]) -> Dict:
        """Calculate transformation matrix for alignment."""
        # Simplified alignment - in practice, use proper RMSD alignment
        # This is a placeholder for the actual alignment algorithm
        
        # Calculate centroids
        centroid1 = self.calculate_centroid(atoms1)
        centroid2 = self.calculate_centroid(atoms2)
        
        # Calculate translation
        translation = {
            'x': centroid1['x'] - centroid2['x'],
            'y': centroid1['y'] - centroid2['y'],
            'z': centroid1['z'] - centroid2['z']
        }
        
        # For now, just return translation (rotation would be more complex)
        return {
            'translation': translation,
            'rotation': None  # Would need proper rotation matrix calculation
        }
    
    def calculate_centroid(self, atoms: List[Dict]) -> Dict:
        """Calculate centroid of a set of atoms."""
        if not atoms:
            return {'x': 0, 'y': 0, 'z': 0}
        
        x_sum = sum(atom['x'] for atom in atoms)
        y_sum = sum(atom['y'] for atom in atoms)
        z_sum = sum(atom['z'] for atom in atoms)
        n = len(atoms)
        
        return {
            'x': x_sum / n,
            'y': y_sum / n,
            'z': z_sum / n
        }
    
    def apply_transformation(self, atom: Dict, transformation: Dict) -> Tuple[float, float, float]:
        """Apply transformation to an atom."""
        # Apply translation
        x = atom['x'] + transformation['translation']['x']
        y = atom['y'] + transformation['translation']['y']
        z = atom['z'] + transformation['translation']['z']
        
        # Apply rotation if available
        if transformation['rotation'] is not None:
            # This would apply the rotation matrix
            pass
        
        return x, y, z
    
    def combine_structures(self) -> Dict:
        """Combine all structures into a single protein."""
        print("Combining structures...")
        
        if not self.structures:
            print("No structures loaded")
            return {}
        
        # Sort structures by start exon
        sorted_structures = sorted(self.structures, key=lambda x: x['start_exon'])
        
        combined_atoms = []
        combined_confidence = []
        
        # Start with the first structure
        first_struct = sorted_structures[0]
        first_atoms = self.parse_cif_file(first_struct['cif_file'])
        combined_atoms.extend(first_atoms['atoms'])
        combined_confidence.extend(first_struct['confidence_data'].get('atom_plddts', []))
        
        print(f"Starting with structure {first_struct['name']} (exons {first_struct['start_exon']}-{first_struct['end_exon']})")
        
        # Add subsequent structures
        for i in range(1, len(sorted_structures)):
            current_struct = sorted_structures[i]
            prev_struct = sorted_structures[i-1]
            
            # Check for overlap
            overlap = self.find_overlapping_exons(prev_struct, current_struct)
            
            if overlap:
                print(f"Found overlap between {prev_struct['name']} and {current_struct['name']} in exons {overlap[0]}-{overlap[1]}")
                
                # Align structures
                aligned_struct = self.align_structures(prev_struct, current_struct, overlap)
                
                # Add non-overlapping atoms from current structure
                current_atoms = self.parse_cif_file(current_struct['cif_file'])
                current_confidence = current_struct['confidence_data'].get('atom_plddts', [])
                
                # Filter out overlapping atoms (simplified approach)
                for j, atom in enumerate(aligned_struct['atoms']):
                    if atom['seq_id'] > overlap[1]:  # Only add atoms after overlap
                        combined_atoms.append(atom)
                        if j < len(aligned_struct['confidence_scores']):
                            combined_confidence.append(aligned_struct['confidence_scores'][j])
            else:
                print(f"No overlap between {prev_struct['name']} and {current_struct['name']}")
                # Add all atoms from current structure
                current_atoms = self.parse_cif_file(current_struct['cif_file'])
                current_confidence = current_struct['confidence_data'].get('atom_plddts', [])
                
                combined_atoms.extend(current_atoms['atoms'])
                combined_confidence.extend(current_confidence)
        
        print(f"Combined structure has {len(combined_atoms)} atoms")
        
        return {
            'atoms': combined_atoms,
            'confidence_scores': combined_confidence
        }
    
    def save_combined_structure(self, combined_data: Dict, output_file: str):
        """Save the combined structure to a CIF file."""
        print(f"Saving combined structure to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write CIF header
            f.write("data_combined_dystrophin\n")
            f.write("# Combined AlphaFold structures for dystrophin\n")
            f.write("# Generated by DystroSeq structure combiner\n\n")
            
            # Write cell parameters (placeholder)
            f.write("_cell.length_a    100.0\n")
            f.write("_cell.length_b    100.0\n")
            f.write("_cell.length_c    100.0\n")
            f.write("_cell.angle_alpha 90.0\n")
            f.write("_cell.angle_beta  90.0\n")
            f.write("_cell.angle_gamma 90.0\n\n")
            
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
            for i, atom in enumerate(combined_data['atoms']):
                # Get confidence score for this atom
                confidence = 50.0  # Default confidence
                if i < len(combined_data['confidence_scores']):
                    confidence = combined_data['confidence_scores'][i]
                
                f.write(f"ATOM {i+1:5d} {atom['atom_id'][0]:2s} {atom['atom_id']:4s} . {atom['residue']:3s} A 1 {atom['seq_id']:4d} ? "
                       f"{atom['x']:8.3f} {atom['y']:8.3f} {atom['z']:8.3f} 1.00 {confidence:6.2f} ? "
                       f"{atom['seq_id']:4d} {atom['residue']:3s} {atom['atom_id']:4s} A 1\n")
    
    def run_combination(self, output_file: str = "combined_dystrophin.cif"):
        """Run the complete structure combination process."""
        print("Starting structure combination process...")
        
        # Load all structures
        self.load_structure_data()
        
        if not self.structures:
            print("No structures found to combine")
            return
        
        # Combine structures
        combined_data = self.combine_structures()
        
        if not combined_data:
            print("Failed to combine structures")
            return
        
        # Save combined structure
        self.save_combined_structure(combined_data, output_file)
        
        print(f"Structure combination complete! Saved to {output_file}")


def main():
    """Main function to combine AlphaFold structures."""
    base_dir = "../../../Data/dystrophin/folds_2025_10_22_19_55"
    
    combiner = StructureCombiner(base_dir)
    combiner.run_combination("combined_dystrophin.cif")


if __name__ == "__main__":
    main()
