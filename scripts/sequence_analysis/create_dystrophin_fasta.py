#!/usr/bin/env python3
"""
Create a FASTA file for dystrophin amino acids with underscores for missing positions
"""

def create_dystrophin_fasta():
    # Read the updated exon report
    with open('dystrophin_analysis_fixed/exon_report.txt', 'r') as f:
        lines = f.readlines()
    
    # Read the full protein sequence
    with open('dystrophin_analysis_fixed/protein.fasta', 'r') as f:
        protein_lines = f.readlines()
        full_protein = ''.join([line.strip() for line in protein_lines[1:]])  # Skip header
    
    # Create a mapping of amino acid positions to exons
    aa_to_exon = {}
    exon_sequences = {}
    
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) >= 9:
            exon_num = int(parts[0])
            start_idx = parts[6] if parts[6] != '' else None
            end_idx = parts[7] if parts[7] != '' else None
            protein_seq = parts[8] if parts[8] != '' else ''
            
            if start_idx and end_idx and protein_seq != 'N/A':
                start_pos = int(start_idx) - 1  # Convert to 0-based indexing
                end_pos = int(end_idx) - 1
                
                # Map amino acids to this exon
                for i in range(start_pos, end_pos + 1):
                    if i < len(full_protein):
                        aa_to_exon[i] = exon_num
                
                # Store the exon sequence
                exon_sequences[exon_num] = protein_seq
    
    # Create the dystroseq-determined sequence with underscores for missing positions
    dystroseq_sequence = []
    for i, aa in enumerate(full_protein):
        if i in aa_to_exon:
            dystroseq_sequence.append(aa)
        else:
            dystroseq_sequence.append('_')
    
    dystroseq_seq = ''.join(dystroseq_sequence)
    
    # Count statistics
    total_aa = len(full_protein)
    assigned_aa = sum(1 for aa in dystroseq_seq if aa != '_')
    missing_aa = total_aa - assigned_aa
    
    # Create FASTA content
    fasta_content = f">Dystrophin_Dp427m_dystroseq_analysis\n"
    fasta_content += f"# Total amino acids: {total_aa}\n"
    fasta_content += f"# Assigned by dystroseq: {assigned_aa}\n"
    fasta_content += f"# Missing (underscores): {missing_aa}\n"
    fasta_content += f"# Coverage: {assigned_aa/total_aa*100:.1f}%\n"
    fasta_content += f"# Underscores indicate amino acids not assigned to any exon\n\n"
    
    # Add the sequence in 60-character lines
    for i in range(0, len(dystroseq_seq), 60):
        fasta_content += dystroseq_seq[i:i+60] + '\n'
    
    # Write to file
    with open('dystrophin_dystroseq_analysis.fasta', 'w') as f:
        f.write(fasta_content)
    
    print("Dystrophin FASTA file created: dystrophin_dystroseq_analysis.fasta")
    print(f"Total amino acids: {total_aa}")
    print(f"Assigned by dystroseq: {assigned_aa}")
    print(f"Missing (underscores): {missing_aa}")
    print(f"Coverage: {assigned_aa/total_aa*100:.1f}%")
    
    # Show some examples of missing regions
    print("\nMissing amino acid regions (showing first 10):")
    missing_regions = []
    in_missing = False
    start_pos = 0
    
    for i, aa in enumerate(dystroseq_seq):
        if aa == '_' and not in_missing:
            start_pos = i + 1  # Convert to 1-based
            in_missing = True
        elif aa != '_' and in_missing:
            missing_regions.append((start_pos, i))
            in_missing = False
    
    # Handle case where sequence ends with missing amino acids
    if in_missing:
        missing_regions.append((start_pos, len(dystroseq_seq)))
    
    for i, (start, end) in enumerate(missing_regions[:10]):
        print(f"  Region {i+1}: positions {start}-{end} ({end-start+1} amino acids)")
    
    if len(missing_regions) > 10:
        print(f"  ... and {len(missing_regions)-10} more regions")
    
    # Show first 100 characters of the sequence
    print(f"\nFirst 100 amino acids of dystroseq analysis:")
    print(dystroseq_seq[:100])
    
    # Show last 100 characters
    print(f"\nLast 100 amino acids of dystroseq analysis:")
    print(dystroseq_seq[-100:])

if __name__ == '__main__':
    create_dystrophin_fasta()
