#!/usr/bin/env python3
"""
Analyze dystrophin exon protein coordinates from dystroseq output
"""

def analyze_dystrophin_exons():
    # Read the exon report
    with open('dystrophin_analysis_fixed/exon_report.txt', 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[0].strip().split('\t')
    print('Dystrophin Exon Protein Coordinates:')
    print('=' * 100)
    print(f"{'Exon':<6} {'Start_AA':<8} {'End_AA':<8} {'AA_Count':<8} {'Protein_Sequence':<30}")
    print('-' * 100)
    
    total_aa = 0
    exon_count = 0
    
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) >= 9:
            exon_num = parts[0]
            start_idx = parts[6] if parts[6] != '' else 'N/A'
            end_idx = parts[7] if parts[7] != '' else 'N/A'
            protein_seq = parts[8] if parts[8] != '' else 'N/A'
            
            if start_idx != 'N/A' and end_idx != 'N/A':
                aa_count = int(end_idx) - int(start_idx) + 1
                total_aa += aa_count
                exon_count += 1
                print(f"{exon_num:<6} {start_idx:<8} {end_idx:<8} {aa_count:<8} {protein_seq[:30]:<30}")
            else:
                print(f"{exon_num:<6} {'N/A':<8} {'N/A':<8} {'0':<8} {'N/A':<30}")
    
    print('-' * 100)
    print(f'Total amino acids from {exon_count} exons: {total_aa}')
    print('=' * 100)
    
    # Read the full protein sequence
    with open('dystrophin_analysis_fixed/protein.fasta', 'r') as f:
        lines = f.readlines()
        protein_seq = ''.join([line.strip() for line in lines[1:]])  # Skip header
    
    print(f'\nFull protein sequence length: {len(protein_seq)} amino acids')
    print(f'Difference: {len(protein_seq) - total_aa} amino acids')
    
    # Show first 50 amino acids
    print(f'\nFirst 50 amino acids: {protein_seq[:50]}')
    print(f'Last 50 amino acids: {protein_seq[-50:]}')
    
    # Known dystrophin protein length (from literature)
    known_length = 3685  # Dystrophin Dp427m isoform
    print(f'\nKnown dystrophin Dp427m length: {known_length} amino acids')
    print(f'Our calculation vs known: {len(protein_seq)} vs {known_length} (diff: {len(protein_seq) - known_length})')

if __name__ == '__main__':
    analyze_dystrophin_exons()
