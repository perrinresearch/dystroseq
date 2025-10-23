"""
Debug the sequence matching issue by examining the sequences in detail.
"""

import os
from pathlib import Path


def analyze_sequences():
    """Analyze the sequences to understand the mismatch."""
    
    print("=" * 60)
    print("SEQUENCE ANALYSIS AND DEBUGGING")
    print("=" * 60)
    
    # Load reference sequence
    print("\n1. Loading reference sequence...")
    with open("dystrophin_dystroseq_analysis.fasta", 'r') as f:
        lines = f.readlines()
    
    reference_seq = ""
    for line in lines:
        if not line.startswith('>') and not line.startswith('#'):
            reference_seq += line.strip()
    
    print(f"Reference sequence length: {len(reference_seq)}")
    print(f"First 50 amino acids: {reference_seq[:50]}")
    print(f"Last 50 amino acids: {reference_seq[-50:]}")
    
    # Extract structure sequence
    print("\n2. Extracting structure sequence...")
    structure_seq = ""
    residues = {}
    
    with open("dystrophin_aligned.pdb", 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    seq_id = int(line[22:26].strip())
                    
                    if atom_name == 'CA':
                        if seq_id not in residues:
                            residues[seq_id] = residue_name
                            
                except (ValueError, IndexError):
                    continue
    
    # Sort by sequence ID and build sequence
    sorted_residues = sorted(residues.items())
    for seq_id, residue_name in sorted_residues:
        aa_code = three_to_one_letter(residue_name)
        structure_seq += aa_code
    
    print(f"Structure sequence length: {len(structure_seq)}")
    print(f"First 50 amino acids: {structure_seq[:50]}")
    print(f"Last 50 amino acids: {structure_seq[-50:]}")
    
    # Find the best alignment manually
    print("\n3. Finding best alignment...")
    best_start = 0
    best_score = 0
    best_match = ""
    
    for start_pos in range(len(reference_seq) - len(structure_seq) + 1):
        ref_segment = reference_seq[start_pos:start_pos + len(structure_seq)]
        
        matches = 0
        for i, (a, b) in enumerate(zip(structure_seq, ref_segment)):
            if a == b:
                matches += 1
        
        score = matches / len(structure_seq)
        if score > best_score:
            best_score = score
            best_start = start_pos
            best_match = ref_segment
    
    print(f"Best alignment found at position {best_start+1}")
    print(f"Match score: {best_score:.3f}")
    
    # Show detailed comparison
    print("\n4. Detailed sequence comparison:")
    print("Position | Structure | Reference | Match")
    print("-" * 40)
    
    mismatches = []
    for i in range(min(50, len(structure_seq))):  # Show first 50 positions
        struct_aa = structure_seq[i]
        ref_aa = best_match[i]
        match = "OK" if struct_aa == ref_aa else "NO"
        print(f"{i+1:8d} | {struct_aa:9s} | {ref_aa:9s} | {match}")
        
        if struct_aa != ref_aa:
            mismatches.append((i+1, struct_aa, ref_aa))
    
    print(f"\nFound {len(mismatches)} mismatches in first 50 positions")
    if len(mismatches) > 0:
        print("First few mismatches:")
        for pos, struct_aa, ref_aa in mismatches[:10]:
            print(f"  Position {pos}: {struct_aa} vs {ref_aa}")
    
    # Check if this is a partial sequence
    print(f"\n5. Sequence coverage analysis:")
    print(f"Reference sequence: {len(reference_seq)} amino acids")
    print(f"Structure sequence: {len(structure_seq)} amino acids")
    print(f"Coverage: {len(structure_seq)/len(reference_seq)*100:.1f}%")
    
    # Check for common patterns
    print(f"\n6. Pattern analysis:")
    print(f"Structure starts with: {structure_seq[:10]}")
    print(f"Structure ends with: {structure_seq[-10:]}")
    
    # Look for the specific overlapping sequence mentioned
    overlap_seq = "ITVSLAQGYERTSSPKPRFKSYAYTQAAYVTTSDPTRSPFPSQ"
    print(f"\n7. Looking for overlap sequence: {overlap_seq}")
    
    if overlap_seq in structure_seq:
        pos = structure_seq.find(overlap_seq)
        print(f"Found overlap sequence at position {pos+1} in structure")
    else:
        print("Overlap sequence not found in structure")
    
    if overlap_seq in reference_seq:
        pos = reference_seq.find(overlap_seq)
        print(f"Found overlap sequence at position {pos+1} in reference")
    else:
        print("Overlap sequence not found in reference")


def three_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_map.get(three_letter.upper(), 'X')


def main():
    """Main function."""
    analyze_sequences()


if __name__ == "__main__":
    main()
