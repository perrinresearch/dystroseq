#!/usr/bin/env python3
"""
Debug codon assignment to identify where amino acids are being dropped
"""

def debug_codon_assignment():
    # Read the exon report
    with open('dystrophin_analysis_updated/exon_report.txt', 'r') as f:
        lines = f.readlines()
    
    # Read the full protein sequence
    with open('dystrophin_analysis_updated/protein.fasta', 'r') as f:
        protein_lines = f.readlines()
        full_protein = ''.join([line.strip() for line in protein_lines[1:]])
    
    print("DEBUGGING CODON ASSIGNMENT")
    print("=" * 80)
    
    # Parse exon data
    exons = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) >= 9:
            exon_num = int(parts[0])
            start_idx = parts[6] if parts[6] != '' else None
            end_idx = parts[7] if parts[7] != '' else None
            protein_seq = parts[8] if parts[8] != '' else ''
            
            if start_idx and end_idx and protein_seq != 'N/A':
                exons.append({
                    'exon': exon_num,
                    'start': int(start_idx),
                    'end': int(end_idx),
                    'seq': protein_seq,
                    'length': int(end_idx) - int(start_idx) + 1
                })
    
    print(f"Found {len(exons)} exons with protein assignments")
    print()
    
    # Check for gaps between exons
    print("CHECKING FOR GAPS BETWEEN EXONS:")
    print("-" * 50)
    
    total_assigned = 0
    for i, exon in enumerate(exons):
        total_assigned += exon['length']
        print(f"Exon {exon['exon']:2d}: positions {exon['start']:4d}-{exon['end']:4d} ({exon['length']:2d} aa) - {exon['seq'][:20]}...")
        
        # Check gap to next exon
        if i < len(exons) - 1:
            next_exon = exons[i + 1]
            gap = next_exon['start'] - exon['end'] - 1
            if gap > 0:
                print(f"  *** GAP: {gap} amino acids between exon {exon['exon']} and {next_exon['exon']} ***")
                print(f"      Missing positions: {exon['end'] + 1} to {next_exon['start'] - 1}")
                print(f"      Missing sequence: {full_protein[exon['end']:next_exon['start']-1]}")
    
    print()
    print(f"Total assigned amino acids: {total_assigned}")
    print(f"Full protein length: {len(full_protein)}")
    print(f"Missing: {len(full_protein) - total_assigned}")
    
    # Check specific missing positions
    print("\nCHECKING SPECIFIC MISSING POSITIONS:")
    print("-" * 50)
    
    missing_positions = [494, 764, 2146, 2254, 2739, 2849, 3362, 3599]
    
    for pos in missing_positions:
        print(f"\nPosition {pos}: '{full_protein[pos-1]}'")
        
        # Find which exons are nearby
        nearby_exons = []
        for exon in exons:
            if abs(exon['start'] - pos) <= 5 or abs(exon['end'] - pos) <= 5:
                nearby_exons.append(exon)
        
        if nearby_exons:
            print(f"  Nearby exons:")
            for exon in nearby_exons:
                print(f"    Exon {exon['exon']}: {exon['start']}-{exon['end']}")
        else:
            print(f"  No nearby exons found")
    
    # Check if there are any overlapping exons
    print("\nCHECKING FOR OVERLAPPING EXONS:")
    print("-" * 50)
    
    for i in range(len(exons) - 1):
        current = exons[i]
        next_exon = exons[i + 1]
        
        if current['end'] >= next_exon['start']:
            print(f"OVERLAP: Exon {current['exon']} ({current['start']}-{current['end']}) overlaps with Exon {next_exon['exon']} ({next_exon['start']}-{next_exon['end']})")
    
    # Check codon boundaries
    print("\nCHECKING CODON BOUNDARIES:")
    print("-" * 50)
    
    for i in range(len(exons) - 1):
        current = exons[i]
        next_exon = exons[i + 1]
        
        # Check if there's a codon split between exons
        current_end_aa = current['end']
        next_start_aa = next_exon['start']
        
        if next_start_aa - current_end_aa == 1:
            print(f"Consecutive exons: {current['exon']} ends at {current_end_aa}, {next_exon['exon']} starts at {next_start_aa}")
            print(f"  This should be continuous - no gap expected")
        elif next_start_aa - current_end_aa > 1:
            gap = next_start_aa - current_end_aa - 1
            print(f"Gap between {current['exon']} and {next_exon['exon']}: {gap} amino acids")
            print(f"  Missing: {full_protein[current_end_aa:next_start_aa-1]}")

if __name__ == '__main__':
    debug_codon_assignment()

