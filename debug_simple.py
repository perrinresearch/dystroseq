#!/usr/bin/env python3
"""
Simple debug to understand the coordinate calculation issue
"""

def debug_simple():
    # Read the exon report
    with open('dystrophin_analysis_updated/exon_report.txt', 'r') as f:
        lines = f.readlines()
    
    # Read the full protein sequence
    with open('dystrophin_analysis_updated/protein.fasta', 'r') as f:
        protein_lines = f.readlines()
        full_protein = ''.join([line.strip() for line in protein_lines[1:]])
    
    print("SIMPLE DEBUG - CHECKING COORDINATE CALCULATION")
    print("=" * 60)
    
    # Look at the first few exons to understand the pattern
    for i, line in enumerate(lines[1:6]):  # First 5 exons
        parts = line.strip().split('\t')
        if len(parts) >= 9:
            exon_num = int(parts[0])
            start_idx = int(parts[6]) if parts[6] != '' else 0
            end_idx = int(parts[7]) if parts[7] != '' else 0
            protein_seq = parts[8] if parts[8] != '' else ''
            
            print(f"Exon {exon_num}: {start_idx}-{end_idx} ({end_idx-start_idx+1} aa)")
            print(f"  Expected: {full_protein[start_idx-1:end_idx]}")
            print(f"  Actual:   {protein_seq}")
            print(f"  Match:    {full_protein[start_idx-1:end_idx] == protein_seq}")
            print()
    
    # Check the specific gaps
    print("CHECKING SPECIFIC GAPS:")
    print("-" * 40)
    
    gaps = [
        (12, 13, 494),
        (18, 19, 764), 
        (44, 45, 2146),
        (46, 47, 2254),
        (55, 56, 2739),
        (57, 58, 2849),
        (69, 70, 3362),
        (75, 76, 3599)
    ]
    
    for prev_exon, next_exon, missing_pos in gaps:
        print(f"Gap between exon {prev_exon} and {next_exon}:")
        print(f"  Missing position: {missing_pos}")
        print(f"  Missing amino acid: '{full_protein[missing_pos-1]}'")
        print(f"  Context: ...{full_protein[missing_pos-3:missing_pos+2]}...")
        print()

if __name__ == '__main__':
    debug_simple()

