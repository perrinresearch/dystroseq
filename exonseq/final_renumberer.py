"""
Final structure renumbering to match the full-length dystrophin sequence.
"""

import os
from pathlib import Path
from typing import List, Dict, Tuple


def three_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_map.get(three_letter.upper(), 'X')


def load_reference_sequence(fasta_file: str) -> str:
    """Load the reference protein sequence from FASTA file."""
    print(f"Loading reference sequence from {fasta_file}")
    
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    reference_seq = ""
    for line in lines:
        if not line.startswith('>') and not line.startswith('#'):
            reference_seq += line.strip()
    
    print(f"Reference sequence length: {len(reference_seq)} amino acids")
    return reference_seq


def extract_structure_sequence(pdb_file: str) -> Tuple[str, List[Dict]]:
    """Extract amino acid sequence from PDB file."""
    print(f"Extracting sequence from {pdb_file}")
    
    residues = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    seq_id = int(line[22:26].strip())
                    
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
        aa_code = three_to_one_letter(residue_data['residue_name'])
        sequence += aa_code
        residue_info.append({
            'seq_id': seq_id,
            'residue_name': residue_data['residue_name'],
            'aa_code': aa_code,
            'atoms': residue_data['atoms']
        })
    
    print(f"Extracted sequence: {len(sequence)} amino acids")
    return sequence, residue_info


def find_sequence_match(structure_seq: str, reference_seq: str) -> Tuple[int, int, float]:
    """Find the best match between structure and reference sequences."""
    print("Finding sequence alignment...")
    
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


def renumber_structure(pdb_file: str, output_file: str, start_pos: int):
    """Renumber the structure to match the reference sequence."""
    print(f"Renumbering structure to match reference positions {start_pos+1}-{start_pos+413}")
    
    with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                try:
                    # Parse current line
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    old_seq_id = int(line[22:26].strip())
                    
                    # Calculate new sequence ID (1-based)
                    new_seq_id = start_pos + old_seq_id
                    
                    # Reconstruct the line with new sequence ID
                    new_line = (line[:22] + f"{new_seq_id:4d}" + line[26:])
                    f_out.write(new_line)
                    
                except (ValueError, IndexError):
                    f_out.write(line)
            else:
                f_out.write(line)
    
    print(f"Renumbered structure saved to {output_file}")


def validate_sequence_match(pdb_file: str, reference_seq: str, start_pos: int) -> bool:
    """Validate that the renumbered structure matches the reference sequence."""
    print("Validating sequence match...")
    
    # Extract sequence from renumbered structure
    structure_seq, _ = extract_structure_sequence(pdb_file)
    
    # Get corresponding reference segment
    ref_segment = reference_seq[start_pos:start_pos + len(structure_seq)]
    
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


def main():
    """Main function to renumber the structure."""
    print("=" * 60)
    print("FINAL STRUCTURE RENUMBERING")
    print("=" * 60)
    
    # Load reference sequence
    reference_seq = load_reference_sequence("dystrophin_dystroseq_analysis.fasta")
    
    # Extract structure sequence
    structure_seq, residue_info = extract_structure_sequence("dystrophin_aligned.pdb")
    
    # Find best match
    start_pos, end_pos, score = find_sequence_match(structure_seq, reference_seq)
    
    if start_pos is None:
        print("Error: Could not find sequence match!")
        return
    
    print(f"Structure matches reference positions {start_pos+1}-{end_pos}")
    print(f"Match score: {score:.3f}")
    
    # Renumber structure
    renumber_structure("dystrophin_aligned.pdb", "dystrophin_final.pdb", start_pos)
    
    # Validate the result
    success = validate_sequence_match("dystrophin_final.pdb", reference_seq, start_pos)
    
    if success:
        print(f"\nSuccess! Final structure saved to dystrophin_final.pdb")
        print(f"Structure now matches reference positions {start_pos+1}-{end_pos}")
        print(f"Residue numbering: {start_pos+1} to {end_pos}")
    else:
        print(f"\nValidation failed. Please check the sequences.")


if __name__ == "__main__":
    main()
