"""
Validate and renumber the combined structure against the full-length dystrophin sequence.
"""

import os
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SequenceValidator:
    """Validates and renumbers protein structures against reference sequences."""
    
    def __init__(self):
        self.reference_sequence = None
        self.structure_sequence = None
        
    def load_reference_sequence(self, fasta_file: str) -> str:
        """Load the reference protein sequence from FASTA file."""
        print(f"Loading reference sequence from {fasta_file}")
        
        if not os.path.exists(fasta_file):
            print(f"Error: Reference FASTA file not found: {fasta_file}")
            return None
        
        try:
            # Try to find the dystrophin protein sequence
            for record in SeqIO.parse(fasta_file, "fasta"):
                if "dystrophin" in record.description.lower() or "dmd" in record.description.lower():
                    self.reference_sequence = str(record.seq)
                    print(f"Found dystrophin sequence: {len(self.reference_sequence)} amino acids")
                    return self.reference_sequence
            
            # If no specific dystrophin sequence found, use the first one
            record = next(SeqIO.parse(fasta_file, "fasta"))
            self.reference_sequence = str(record.seq)
            print(f"Using first sequence: {len(self.reference_sequence)} amino acids")
            return self.reference_sequence
            
        except Exception as e:
            print(f"Error loading reference sequence: {e}")
            return None
    
    def extract_structure_sequence(self, pdb_file: str) -> Tuple[str, List[Dict]]:
        """Extract amino acid sequence from PDB file."""
        print(f"Extracting sequence from {pdb_file}")
        
        if not os.path.exists(pdb_file):
            print(f"Error: PDB file not found: {pdb_file}")
            return None, []
        
        # Dictionary to store residue information
        residues = {}
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        # Parse PDB line
                        atom_name = line[12:16].strip()
                        residue_name = line[17:20].strip()
                        seq_id = int(line[22:26].strip())
                        
                        # Only process CA atoms to get unique residues
                        if atom_name == 'CA':
                            if seq_id not in residues:
                                residues[seq_id] = {
                                    'residue_name': residue_name,
                                    'seq_id': seq_id,
                                    'atoms': []
                                }
                            residues[seq_id]['atoms'].append(line.strip())
                            
                    except (ValueError, IndexError):
                        continue
        
        # Sort by sequence ID and extract amino acid sequence
        sorted_residues = sorted(residues.items())
        sequence = ""
        residue_info = []
        
        for seq_id, residue_data in sorted_residues:
            aa_code = self.three_to_one_letter(residue_data['residue_name'])
            sequence += aa_code
            residue_info.append({
                'seq_id': seq_id,
                'residue_name': residue_data['residue_name'],
                'aa_code': aa_code,
                'atoms': residue_data['atoms']
            })
        
        self.structure_sequence = sequence
        print(f"Extracted sequence: {len(sequence)} amino acids")
        return sequence, residue_info
    
    def three_to_one_letter(self, three_letter: str) -> str:
        """Convert three-letter amino acid code to one-letter code."""
        aa_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return aa_map.get(three_letter.upper(), 'X')
    
    def find_sequence_match(self, structure_seq: str, reference_seq: str) -> Tuple[int, int, float]:
        """Find the best match between structure and reference sequences."""
        print("Finding sequence alignment...")
        
        if not structure_seq or not reference_seq:
            return None, None, 0.0
        
        best_match_start = 0
        best_match_end = 0
        best_score = 0.0
        
        # Try different starting positions in the reference sequence
        for start_pos in range(len(reference_seq) - len(structure_seq) + 1):
            ref_segment = reference_seq[start_pos:start_pos + len(structure_seq)]
            
            # Calculate match score
            matches = sum(1 for a, b in zip(structure_seq, ref_segment) if a == b)
            score = matches / len(structure_seq)
            
            if score > best_score:
                best_score = score
                best_match_start = start_pos
                best_match_end = start_pos + len(structure_seq)
        
        print(f"Best match: positions {best_match_start+1}-{best_match_end} (score: {best_score:.3f})")
        return best_match_start, best_match_end, best_score
    
    def renumber_structure(self, pdb_file: str, output_file: str, start_pos: int, end_pos: int):
        """Renumber the structure to match the reference sequence."""
        print(f"Renumbering structure to match reference positions {start_pos+1}-{end_pos}")
        
        with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('ATOM'):
                    try:
                        # Parse current line
                        atom_name = line[12:16].strip()
                        residue_name = line[17:20].strip()
                        old_seq_id = int(line[22:26].strip())
                        
                        # Calculate new sequence ID
                        new_seq_id = start_pos + old_seq_id
                        
                        # Reconstruct the line with new sequence ID
                        new_line = (line[:22] + f"{new_seq_id:4d}" + line[26:])
                        f_out.write(new_line)
                        
                    except (ValueError, IndexError):
                        f_out.write(line)
                else:
                    f_out.write(line)
        
        print(f"Renumbered structure saved to {output_file}")
    
    def validate_sequence_match(self, pdb_file: str, reference_seq: str, start_pos: int, end_pos: int) -> bool:
        """Validate that the renumbered structure matches the reference sequence."""
        print("Validating sequence match...")
        
        # Extract sequence from renumbered structure
        structure_seq, _ = self.extract_structure_sequence(pdb_file)
        
        # Get corresponding reference segment
        ref_segment = reference_seq[start_pos:end_pos]
        
        # Compare sequences
        matches = sum(1 for a, b in zip(structure_seq, ref_segment) if a == b)
        total = len(structure_seq)
        match_percentage = (matches / total) * 100
        
        print(f"Sequence validation:")
        print(f"  Structure length: {len(structure_seq)}")
        print(f"  Reference length: {len(ref_segment)}")
        print(f"  Matches: {matches}/{total} ({match_percentage:.1f}%)")
        
        if match_percentage < 95:
            print("WARNING: Low sequence match! Check for errors.")
            return False
        else:
            print("Sequence match validation passed!")
            return True
    
    def process_structure(self, pdb_file: str, reference_fasta: str, output_file: str):
        """Complete process: load reference, extract structure, align, renumber, and validate."""
        print("=" * 60)
        print("SEQUENCE VALIDATION AND RENUMBERING")
        print("=" * 60)
        
        # Load reference sequence
        reference_seq = self.load_reference_sequence(reference_fasta)
        if not reference_seq:
            return False
        
        # Extract structure sequence
        structure_seq, residue_info = self.extract_structure_sequence(pdb_file)
        if not structure_seq:
            return False
        
        # Find best match
        start_pos, end_pos, score = self.find_sequence_match(structure_seq, reference_seq)
        if start_pos is None:
            print("Error: Could not find sequence match!")
            return False
        
        # Renumber structure
        self.renumber_structure(pdb_file, output_file, start_pos, end_pos)
        
        # Validate the result
        success = self.validate_sequence_match(output_file, reference_seq, start_pos, end_pos)
        
        if success:
            print(f"\nSuccess! Renumbered structure saved to {output_file}")
            print(f"Structure now matches reference positions {start_pos+1}-{end_pos}")
        
        return success


def main():
    """Main function to validate and renumber the structure."""
    # Look for reference FASTA files
    possible_fasta_files = [
        "dystrophin_dystroseq_analysis.fasta",
        "dystrophin_analysis_fixed/protein.fasta",
        "dystrophin_analysis/protein.fasta"
    ]
    
    reference_fasta = None
    for fasta_file in possible_fasta_files:
        if os.path.exists(fasta_file):
            reference_fasta = fasta_file
            break
    
    if not reference_fasta:
        print("Error: No reference FASTA file found!")
        print("Please provide a reference dystrophin protein sequence.")
        return
    
    # Process the aligned structure
    validator = SequenceValidator()
    success = validator.process_structure(
        "dystrophin_aligned.pdb",
        reference_fasta,
        "dystrophin_renumbered.pdb"
    )
    
    if success:
        print("\nStructure validation and renumbering complete!")
    else:
        print("\nStructure validation failed. Please check the sequences.")


if __name__ == "__main__":
    main()
