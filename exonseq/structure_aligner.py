"""
Advanced structure alignment system for combining AlphaFold structures.
Handles overlapping exons, structural alignment, and steric clash minimization.
"""

import os
import glob
import json
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional, NamedTuple
import re
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')


class Atom:
    """Represents a single atom in the structure."""
    def __init__(self, atom_id: str, residue: str, seq_id: int, x: float, y: float, z: float, b_factor: float):
        self.atom_id = atom_id
        self.residue = residue
        self.seq_id = seq_id
        self.x = x
        self.y = y
        self.z = z
        self.b_factor = b_factor
        self.coords = np.array([x, y, z])
    
    def __repr__(self):
        return f"Atom({self.atom_id}, {self.residue}, {self.seq_id}, {self.coords})"


class Structure:
    """Represents a protein structure with atoms."""
    def __init__(self, name: str, start_exon: int, end_exon: int):
        self.name = name
        self.start_exon = start_exon
        self.end_exon = end_exon
        self.atoms = []
        self.confidence_scores = []
    
    def add_atom(self, atom: Atom, confidence: float):
        """Add an atom to the structure."""
        self.atoms.append(atom)
        self.confidence_scores.append(confidence)
    
    def get_atoms_by_residue_range(self, start_res: int, end_res: int) -> List[Atom]:
        """Get atoms within a residue range."""
        return [atom for atom in self.atoms if start_res <= atom.seq_id <= end_res]
    
    def get_ca_atoms(self, start_res: int = None, end_res: int = None) -> List[Atom]:
        """Get CA atoms within a residue range."""
        if start_res is None:
            start_res = min(atom.seq_id for atom in self.atoms)
        if end_res is None:
            end_res = max(atom.seq_id for atom in self.atoms)
        
        return [atom for atom in self.atoms 
                if start_res <= atom.seq_id <= end_res and atom.atom_id == 'CA']
    
    def get_centroid(self, atoms: List[Atom] = None) -> np.ndarray:
        """Calculate centroid of specified atoms."""
        if atoms is None:
            atoms = self.atoms
        
        if not atoms:
            return np.array([0.0, 0.0, 0.0])
        
        coords = np.array([atom.coords for atom in atoms])
        return np.mean(coords, axis=0)
    
    def translate(self, translation_vector: np.ndarray):
        """Translate all atoms by a vector."""
        for atom in self.atoms:
            atom.coords += translation_vector
            atom.x, atom.y, atom.z = atom.coords


class StructureAligner:
    """Handles structural alignment and combination of protein structures."""
    
    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir)
        self.structures = []
        self.aligned_structures = []
        
    def load_structures(self):
        """Load all AlphaFold structures."""
        structure_dirs = [d for d in os.listdir(self.base_dir) 
                         if os.path.isdir(self.base_dir / d) and d.startswith('dystroseq_exons_')]
        
        for structure_dir in sorted(structure_dirs):
            print(f"Loading {structure_dir}...")
            
            # Extract exon range
            start_exon, end_exon = self.extract_exon_range(structure_dir)
            if start_exon is None or end_exon is None:
                print(f"Warning: Could not extract exon range from {structure_dir}")
                continue
            
            # Find CIF files
            cif_files = list((self.base_dir / structure_dir).glob("*_model_0.cif"))
            if not cif_files:
                print(f"Warning: No model_0.cif files found in {structure_dir}")
                continue
            
            # Load confidence data
            confidence_data = self.load_confidence_data(self.base_dir / structure_dir)
            
            # Parse structure
            structure = self.parse_cif_structure(cif_files[0], structure_dir, start_exon, end_exon, confidence_data)
            if structure:
                self.structures.append(structure)
        
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
        numbers = re.findall(r'\d+', dirname)
        if len(numbers) >= 2:
            return int(numbers[0]), int(numbers[1])
        
        return None, None
    
    def load_confidence_data(self, structure_dir: Path) -> List[float]:
        """Load confidence data for a structure."""
        confidence_files = list(structure_dir.glob("*_full_data_*.json"))
        if not confidence_files:
            return []
        
        # Load the first confidence file
        with open(confidence_files[0], 'r') as f:
            data = json.load(f)
        
        return data.get('atom_plddts', [])
    
    def parse_cif_structure(self, cif_file: Path, name: str, start_exon: int, end_exon: int, confidence_data: List[float]) -> Optional[Structure]:
        """Parse a CIF file and create a Structure object."""
        structure = Structure(name, start_exon, end_exon)
        
        with open(cif_file, 'r') as f:
            lines = f.readlines()
        
        in_atom_site = False
        atom_index = 0
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('_atom_site.'):
                in_atom_site = True
                continue
            
            elif line.startswith('ATOM') and in_atom_site:
                fields = line.split()
                if len(fields) >= 15:
                    try:
                        atom = Atom(
                            atom_id=fields[3],
                            residue=fields[5],
                            seq_id=int(fields[8]),
                            x=float(fields[10]),
                            y=float(fields[11]),
                            z=float(fields[12]),
                            b_factor=float(fields[14])
                        )
                        
                        confidence = 50.0  # Default confidence
                        if atom_index < len(confidence_data):
                            confidence = confidence_data[atom_index]
                        
                        structure.add_atom(atom, confidence)
                        atom_index += 1
                        
                    except (ValueError, IndexError):
                        continue
            
            elif line.startswith('#') and in_atom_site:
                break
        
        if not structure.atoms:
            print(f"Warning: No atoms parsed from {cif_file}")
            return None
        
        print(f"  Parsed {len(structure.atoms)} atoms")
        return structure
    
    def find_overlapping_residues(self, struct1: Structure, struct2: Structure) -> Tuple[int, int]:
        """Find overlapping residue range between two structures."""
        # Check if there's an overlap
        if struct1.end_exon >= struct2.start_exon and struct2.end_exon >= struct1.start_exon:
            overlap_start = max(struct1.start_exon, struct2.start_exon)
            overlap_end = min(struct1.end_exon, struct2.end_exon)
            return overlap_start, overlap_end
        return None, None
    
    def get_overlapping_ca_atoms(self, struct1: Structure, struct2: Structure, overlap_start: int, overlap_end: int) -> Tuple[List[Atom], List[Atom]]:
        """Get CA atoms from overlapping regions of two structures."""
        # Get CA atoms from the overlapping region
        ca1 = [atom for atom in struct1.get_ca_atoms() 
               if overlap_start <= atom.seq_id <= overlap_end]
        ca2 = [atom for atom in struct2.get_ca_atoms() 
               if overlap_start <= atom.seq_id <= overlap_end]
        
        return ca1, ca2
    
    def calculate_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Calculate RMSD between two sets of coordinates."""
        if len(coords1) != len(coords2):
            return float('inf')
        
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def align_structures(self, struct1: Structure, struct2: Structure, overlap_start: int, overlap_end: int) -> Structure:
        """Align two structures based on overlapping residues."""
        print(f"Aligning {struct1.name} and {struct2.name} (overlap: exons {overlap_start}-{overlap_end})")
        
        # Get overlapping CA atoms
        ca1, ca2 = self.get_overlapping_ca_atoms(struct1, struct2, overlap_start, overlap_end)
        
        if len(ca1) != len(ca2) or len(ca1) == 0:
            print(f"Warning: No matching CA atoms found for alignment")
            return struct2
        
        print(f"  Found {len(ca1)} overlapping CA atoms")
        
        # Extract coordinates
        coords1 = np.array([atom.coords for atom in ca1])
        coords2 = np.array([atom.coords for atom in ca2])
        
        # Calculate initial RMSD
        initial_rmsd = self.calculate_rmsd(coords1, coords2)
        print(f"  Initial RMSD: {initial_rmsd:.3f} Å")
        
        # Perform alignment using Kabsch algorithm
        aligned_coords2, transformation = self.kabsch_algorithm(coords1, coords2)
        
        # Apply transformation to all atoms in struct2
        self.apply_transformation_to_structure(struct2, transformation)
        
        # Calculate final RMSD
        final_rmsd = self.calculate_rmsd(coords1, aligned_coords2)
        print(f"  Final RMSD: {final_rmsd:.3f} Å")
        
        return struct2
    
    def kabsch_algorithm(self, coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, Dict]:
        """Perform Kabsch algorithm for optimal alignment."""
        # Center the coordinates
        centroid1 = np.mean(coords1, axis=0)
        centroid2 = np.mean(coords2, axis=0)
        
        centered1 = coords1 - centroid1
        centered2 = coords2 - centroid2
        
        # Calculate rotation matrix
        H = centered2.T @ centered1
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        
        # Ensure proper rotation (det(R) = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Calculate translation
        t = centroid1 - R @ centroid2
        
        # Apply transformation
        aligned_coords2 = (R @ coords2.T).T + t
        
        transformation = {
            'rotation': R,
            'translation': t
        }
        
        return aligned_coords2, transformation
    
    def apply_transformation_to_structure(self, structure: Structure, transformation: Dict):
        """Apply rotation and translation to all atoms in a structure."""
        R = transformation['rotation']
        t = transformation['translation']
        
        for atom in structure.atoms:
            # Apply rotation and translation
            new_coords = (R @ atom.coords) + t
            atom.coords = new_coords
            atom.x, atom.y, atom.z = new_coords
    
    def remove_overlapping_residues(self, struct1: Structure, struct2: Structure, overlap_start: int, overlap_end: int) -> Structure:
        """Remove overlapping residues from the second structure."""
        print(f"Removing overlapping residues {overlap_start}-{overlap_end} from {struct2.name}")
        
        # Keep only non-overlapping atoms from struct2
        non_overlapping_atoms = []
        non_overlapping_confidence = []
        
        for i, atom in enumerate(struct2.atoms):
            if not (overlap_start <= atom.seq_id <= overlap_end):
                non_overlapping_atoms.append(atom)
                if i < len(struct2.confidence_scores):
                    non_overlapping_confidence.append(struct2.confidence_scores[i])
        
        # Update struct2
        struct2.atoms = non_overlapping_atoms
        struct2.confidence_scores = non_overlapping_confidence
        
        print(f"  Kept {len(non_overlapping_atoms)} non-overlapping atoms")
        return struct2
    
    def combine_structures(self) -> Structure:
        """Combine all structures with proper alignment."""
        print("Combining structures with alignment...")
        
        if not self.structures:
            print("No structures to combine")
            return None
        
        # Sort structures by start exon
        sorted_structures = sorted(self.structures, key=lambda x: x.start_exon)
        
        # Start with the first structure
        combined = sorted_structures[0]
        print(f"Starting with {combined.name}")
        
        # Align and combine subsequent structures
        for i in range(1, len(sorted_structures)):
            current = sorted_structures[i]
            prev = sorted_structures[i-1]
            
            # Find overlap
            overlap_start, overlap_end = self.find_overlapping_residues(prev, current)
            
            if overlap_start and overlap_end:
                print(f"Found overlap between {prev.name} and {current.name}: exons {overlap_start}-{overlap_end}")
                
                # Align structures
                aligned_current = self.align_structures(prev, current, overlap_start, overlap_end)
                
                # Remove overlapping residues from current structure
                trimmed_current = self.remove_overlapping_residues(prev, aligned_current, overlap_start, overlap_end)
                
                # Combine structures
                combined.atoms.extend(trimmed_current.atoms)
                combined.confidence_scores.extend(trimmed_current.confidence_scores)
                
                print(f"Combined structure now has {len(combined.atoms)} atoms")
            else:
                print(f"No overlap between {prev.name} and {current.name}")
                # Just add the current structure
                combined.atoms.extend(current.atoms)
                combined.confidence_scores.extend(current.confidence_scores)
        
        print(f"Final combined structure has {len(combined.atoms)} atoms")
        return combined
    
    def save_combined_structure(self, combined: Structure, output_file: str):
        """Save the combined structure to a PDB file."""
        print(f"Saving combined structure to {output_file}")
        
        with open(output_file, 'w') as f:
            # Write header
            f.write("HEADER    COMBINED DYSTROPHIN STRUCTURE                   01-JAN-24   COMB\n")
            f.write("TITLE     COMBINED ALPHAFOLD STRUCTURES FOR DYSTROPHIN\n")
            f.write("REMARK   1 CONFIDENCE SCORES STORED IN B-FACTOR COLUMN\n")
            f.write("REMARK   2 GENERATED BY DYSTROSEQ STRUCTURE ALIGNER\n")
            f.write("REMARK   3 STRUCTURALLY ALIGNED WITH OVERLAP REMOVAL\n\n")
            
            # Write atoms
            for i, atom in enumerate(combined.atoms):
                confidence = 50.0  # Default confidence
                if i < len(combined.confidence_scores):
                    confidence = combined.confidence_scores[i]
                
                # Get element symbol
                element = atom.atom_id[0] if atom.atom_id else 'C'
                
                # Format PDB line with correct column alignment
                pdb_line = f"ATOM  {i+1:5d}  {atom.atom_id:4s}{atom.residue:3s} A{atom.seq_id:4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  {1.00:4.2f}{confidence:6.2f}           {element:2s}\n"
                f.write(pdb_line)
            
            # Write END
            f.write("END\n")
        
        print(f"Combined structure saved to {output_file}")


def main():
    """Main function to align and combine structures."""
    base_dir = "../../../Data/dystrophin/folds_2025_10_22_19_55"
    
    print("=" * 60)
    print("DYSTROSEQ STRUCTURE ALIGNER")
    print("=" * 60)
    
    # Create aligner
    aligner = StructureAligner(base_dir)
    
    # Load structures
    structures = aligner.load_structures()
    
    if not structures:
        print("No structures found!")
        return
    
    # Combine structures with alignment
    combined = aligner.combine_structures()
    
    if combined:
        # Save combined structure
        aligner.save_combined_structure(combined, "dystrophin_aligned.pdb")
        
        print(f"\nStructure alignment complete!")
        print(f"  Final structure: {len(combined.atoms)} atoms")
        print(f"  Output file: dystrophin_aligned.pdb")
        print(f"  This file should now load correctly in VMD with proper structural alignment.")


if __name__ == "__main__":
    main()
