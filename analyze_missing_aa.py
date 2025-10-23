#!/usr/bin/env python3
"""
Analyze the missing amino acids in the dystroseq analysis
"""

def analyze_missing_amino_acids():
    # Read the full protein sequence
    with open('dystrophin_analysis_updated/protein.fasta', 'r') as f:
        protein_lines = f.readlines()
        full_protein = ''.join([line.strip() for line in protein_lines[1:]])
    
    # Read the dystroseq analysis
    with open('dystrophin_dystroseq_analysis.fasta', 'r') as f:
        lines = f.readlines()
        # Find the sequence part (skip header and comments)
        sequence_lines = []
        in_sequence = False
        for line in lines:
            if line.startswith('>') or line.startswith('#'):
                continue
            if line.strip() and not line.startswith('#'):
                in_sequence = True
            if in_sequence and line.strip():
                sequence_lines.append(line.strip())
        
        dystroseq_seq = ''.join(sequence_lines)
    
    print("MISSING AMINO ACIDS ANALYSIS")
    print("=" * 60)
    
    # Find missing positions
    missing_positions = []
    for i, (orig_aa, dystro_aa) in enumerate(zip(full_protein, dystroseq_seq)):
        if dystro_aa == '_':
            missing_positions.append(i + 1)  # Convert to 1-based
    
    print(f"Total missing amino acids: {len(missing_positions)}")
    print(f"Missing positions: {missing_positions}")
    print()
    
    # Show context for each missing amino acid
    print("CONTEXT OF MISSING AMINO ACIDS:")
    print("-" * 60)
    
    for pos in missing_positions:
        start = max(0, pos - 6)  # Show 5 amino acids before
        end = min(len(full_protein), pos + 5)  # Show 5 amino acids after
        
        context_orig = full_protein[start:end]
        context_dystro = dystroseq_seq[start:end]
        
        # Highlight the missing position
        orig_display = context_orig[:pos-start] + f"[{context_orig[pos-start]}]" + context_orig[pos-start+1:]
        dystro_display = context_dystro[:pos-start] + f"[{context_dystro[pos-start]}]" + context_dystro[pos-start+1:]
        
        print(f"Position {pos:4d}: {orig_display}")
        print(f"          {dystro_display}")
        print()
    
    # Check if missing amino acids are in specific regions
    print("ANALYSIS BY PROTEIN REGIONS:")
    print("-" * 60)
    
    # Define protein regions (approximate)
    regions = [
        (1, 500, "N-terminal domain"),
        (501, 1000, "Actin-binding domain 1"),
        (1001, 2000, "Central rod domain (part 1)"),
        (2001, 3000, "Central rod domain (part 2)"),
        (3001, 3500, "Central rod domain (part 3)"),
        (3501, 3685, "C-terminal domain")
    ]
    
    for start, end, name in regions:
        region_missing = [pos for pos in missing_positions if start <= pos <= end]
        if region_missing:
            print(f"{name}: {len(region_missing)} missing at positions {region_missing}")
        else:
            print(f"{name}: No missing amino acids")
    
    print()
    print("SUMMARY:")
    print(f"- Total protein length: {len(full_protein)} amino acids")
    print(f"- Assigned by dystroseq: {len(full_protein) - len(missing_positions)} amino acids")
    print(f"- Missing amino acids: {len(missing_positions)} amino acids")
    print(f"- Coverage: {(len(full_protein) - len(missing_positions))/len(full_protein)*100:.1f}%")

if __name__ == '__main__':
    analyze_missing_amino_acids()

