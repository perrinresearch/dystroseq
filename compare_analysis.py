#!/usr/bin/env python3
"""
Compare the original and updated dystrophin analysis results
"""

def compare_analyses():
    print("COMPARISON: Original vs Updated Dystrophin Analysis")
    print("=" * 80)
    
    # Read original analysis
    with open('dystrophin_analysis/exon_report.txt', 'r') as f:
        original_lines = f.readlines()
    
    # Read updated analysis  
    with open('dystrophin_analysis_updated/exon_report.txt', 'r') as f:
        updated_lines = f.readlines()
    
    print(f"{'Metric':<30} {'Original':<15} {'Updated':<15} {'Difference':<15}")
    print("-" * 80)
    
    # Calculate totals for both
    orig_total = 0
    upd_total = 0
    
    for i, (orig_line, upd_line) in enumerate(zip(original_lines[1:], updated_lines[1:])):
        orig_parts = orig_line.strip().split('\t')
        upd_parts = upd_line.strip().split('\t')
        
        if len(orig_parts) >= 9 and len(upd_parts) >= 9:
            orig_start = orig_parts[6] if orig_parts[6] != '' else '0'
            orig_end = orig_parts[7] if orig_parts[7] != '' else '0'
            upd_start = upd_parts[6] if upd_parts[6] != '' else '0'
            upd_end = upd_parts[7] if upd_parts[7] != '' else '0'
            
            if orig_start != 'N/A' and orig_end != 'N/A':
                orig_count = int(orig_end) - int(orig_start) + 1
                orig_total += orig_count
            else:
                orig_count = 0
                
            if upd_start != 'N/A' and upd_end != 'N/A':
                upd_count = int(upd_end) - int(upd_start) + 1
                upd_total += upd_count
            else:
                upd_count = 0
            
            if orig_count != upd_count:
                print(f"Exon {i+1:<25} {orig_count:<15} {upd_count:<15} {upd_count - orig_count:<15}")
    
    print("-" * 80)
    print(f"{'TOTAL AMINO ACIDS':<30} {orig_total:<15} {upd_total:<15} {upd_total - orig_total:<15}")
    print("=" * 80)
    
    # Read protein sequences
    with open('dystrophin_analysis/protein.fasta', 'r') as f:
        orig_protein = ''.join([line.strip() for line in f.readlines()[1:]])
    
    with open('dystrophin_analysis_updated/protein.fasta', 'r') as f:
        upd_protein = ''.join([line.strip() for line in f.readlines()[1:]])
    
    print(f"\nProtein sequence lengths:")
    print(f"Original: {len(orig_protein)} amino acids")
    print(f"Updated:  {len(upd_protein)} amino acids")
    print(f"Difference: {len(upd_protein) - len(orig_protein)} amino acids")
    
    print(f"\nCoverage improvement:")
    print(f"Original coverage: {orig_total}/{len(orig_protein)} = {orig_total/len(orig_protein)*100:.1f}%")
    print(f"Updated coverage:  {upd_total}/{len(upd_protein)} = {upd_total/len(upd_protein)*100:.1f}%")
    print(f"Improvement: {upd_total/len(upd_protein)*100 - orig_total/len(orig_protein)*100:.1f} percentage points")

if __name__ == '__main__':
    compare_analyses()

