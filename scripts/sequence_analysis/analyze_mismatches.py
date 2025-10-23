"""
Analyze the specific mismatches between structure and reference sequences.
"""

import os


def three_to_one_letter(three_letter: str) -> str:
    """Convert three-letter amino acid code to one-letter code."""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_map.get(three_letter.upper(), 'X')


def analyze_mismatches():
    """Analyze mismatches between structure and reference sequences."""
    
    print("=" * 60)
    print("MISMATCH ANALYSIS")
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
    
    # Extract structure sequence
    print("\n2. Extracting structure sequence...")
    structure_seq = ""
    residues = {}
    
    with open("dystrophin_final.pdb", 'r') as f:
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
    
    # Compare sequences
    print("\n3. Comparing sequences...")
    mismatches = []
    matches = 0
    
    for i in range(len(structure_seq)):
        struct_aa = structure_seq[i]
        ref_aa = reference_seq[i]
        
        if struct_aa == ref_aa:
            matches += 1
        else:
            mismatches.append((i+1, struct_aa, ref_aa))
    
    match_percentage = (matches / len(structure_seq)) * 100
    print(f"Matches: {matches}/{len(structure_seq)} ({match_percentage:.1f}%)")
    print(f"Mismatches: {len(mismatches)}")
    
    # Show first 20 mismatches
    print("\n4. First 20 mismatches:")
    print("Position | Structure | Reference | Difference")
    print("-" * 45)
    
    for i, (pos, struct_aa, ref_aa) in enumerate(mismatches[:20]):
        print(f"{pos:8d} | {struct_aa:9s} | {ref_aa:9s} | {struct_aa} -> {ref_aa}")
    
    # Analyze mismatch patterns
    print(f"\n5. Mismatch pattern analysis:")
    
    # Count different types of mismatches
    mismatch_types = {}
    for pos, struct_aa, ref_aa in mismatches:
        key = f"{struct_aa}->{ref_aa}"
        mismatch_types[key] = mismatch_types.get(key, 0) + 1
    
    print("Most common mismatches:")
    sorted_mismatches = sorted(mismatch_types.items(), key=lambda x: x[1], reverse=True)
    for mismatch, count in sorted_mismatches[:10]:
        print(f"  {mismatch}: {count} occurrences")
    
    # Check if mismatches are clustered
    print(f"\n6. Mismatch clustering analysis:")
    mismatch_positions = [pos for pos, _, _ in mismatches]
    
    # Find clusters of mismatches
    clusters = []
    current_cluster = [mismatch_positions[0]]
    
    for i in range(1, len(mismatch_positions)):
        if mismatch_positions[i] - mismatch_positions[i-1] <= 5:  # Within 5 positions
            current_cluster.append(mismatch_positions[i])
        else:
            if len(current_cluster) > 1:
                clusters.append(current_cluster)
            current_cluster = [mismatch_positions[i]]
    
    if len(current_cluster) > 1:
        clusters.append(current_cluster)
    
    print(f"Found {len(clusters)} clusters of mismatches:")
    for i, cluster in enumerate(clusters[:5]):  # Show first 5 clusters
        print(f"  Cluster {i+1}: positions {cluster[0]}-{cluster[-1]} ({len(cluster)} mismatches)")
    
    # Check for specific patterns
    print(f"\n7. Pattern analysis:")
    
    # Look for the overlap sequence
    overlap_seq = "ITVSLAQGYERTSSPKPRFKSYAYTQAAYVTTSDPTRSPFPSQ"
    print(f"Looking for overlap sequence: {overlap_seq}")
    
    if overlap_seq in structure_seq:
        pos = structure_seq.find(overlap_seq)
        print(f"Found overlap sequence at position {pos+1} in structure")
        
        # Check if it matches in reference
        ref_segment = reference_seq[pos:pos+len(overlap_seq)]
        if overlap_seq == ref_segment:
            print("Overlap sequence matches perfectly in reference")
        else:
            print("Overlap sequence has mismatches in reference")
            print(f"Structure: {overlap_seq}")
            print(f"Reference: {ref_segment}")
    else:
        print("Overlap sequence not found in structure")
    
    # Check for common protein motifs
    print(f"\n8. Checking for common motifs:")
    motifs = [
        "GXGXXG",  # Common in kinases
        "HEXXH",   # Common in metalloproteases
        "CXXC",    # Common in disulfide bonds
    ]
    
    for motif in motifs:
        if motif in structure_seq:
            pos = structure_seq.find(motif)
            print(f"Found motif {motif} at position {pos+1}")
        else:
            print(f"Motif {motif} not found")


def main():
    """Main function."""
    analyze_mismatches()


if __name__ == "__main__":
    main()
