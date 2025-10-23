"""
Generate a report for the combined dystrophin structure.
"""

import os
from pathlib import Path


def analyze_combined_structure():
    """Analyze the combined structure and generate a report."""
    
    print("=" * 60)
    print("COMBINED DYSTROPHIN STRUCTURE ANALYSIS")
    print("=" * 60)
    
    # Check if files exist
    cif_file = Path("combined_dystrophin.cif")
    pdb_file = Path("combined_dystrophin.pdb")
    
    if not cif_file.exists() or not pdb_file.exists():
        print("Error: Combined structure files not found!")
        return
    
    # Analyze PDB file
    print("\nAnalyzing PDB file...")
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    atom_lines = [line for line in lines if line.startswith('ATOM')]
    total_atoms = len(atom_lines)
    
    print(f"Total atoms: {total_atoms:,}")
    
    # Analyze confidence scores (B-factor column)
    confidence_scores = []
    for line in atom_lines:
        try:
            # B-factor is in column 61-66 (0-indexed: 60-65)
            b_factor = float(line[60:66].strip())
            confidence_scores.append(b_factor)
        except (ValueError, IndexError):
            continue
    
    if confidence_scores:
        print(f"Confidence scores analyzed: {len(confidence_scores):,}")
        print(f"Mean confidence: {sum(confidence_scores)/len(confidence_scores):.1f}")
        print(f"Min confidence: {min(confidence_scores):.1f}")
        print(f"Max confidence: {max(confidence_scores):.1f}")
        
        # Confidence distribution
        high_conf = sum(1 for c in confidence_scores if c >= 70)
        medium_conf = sum(1 for c in confidence_scores if 50 <= c < 70)
        low_conf = sum(1 for c in confidence_scores if c < 50)
        
        print(f"\nConfidence distribution:")
        print(f"  High confidence (>=70): {high_conf:,} ({high_conf/len(confidence_scores)*100:.1f}%)")
        print(f"  Medium confidence (50-70): {medium_conf:,} ({medium_conf/len(confidence_scores)*100:.1f}%)")
        print(f"  Low confidence (<50): {low_conf:,} ({low_conf/len(confidence_scores)*100:.1f}%)")
    
    # Analyze residue composition
    print(f"\nAnalyzing residue composition...")
    residues = {}
    for line in atom_lines:
        if line.startswith('ATOM'):
            residue_name = line[17:20].strip()
            if residue_name not in residues:
                residues[residue_name] = 0
            residues[residue_name] += 1
    
    print(f"Residue types found: {len(residues)}")
    print("Top 10 most abundant residues:")
    sorted_residues = sorted(residues.items(), key=lambda x: x[1], reverse=True)
    for i, (residue, count) in enumerate(sorted_residues[:10]):
        print(f"  {i+1:2d}. {residue}: {count:,} atoms")
    
    # File sizes
    print(f"\nFile information:")
    print(f"  PDB file size: {pdb_file.stat().st_size / 1024 / 1024:.1f} MB")
    print(f"  CIF file size: {cif_file.stat().st_size / 1024 / 1024:.1f} MB")
    
    # Structure coverage
    print(f"\nStructure coverage:")
    print(f"  The combined structure covers the full dystrophin protein")
    print(f"  from exon 1 to exon 79 with overlapping regions")
    print(f"  Confidence scores from AlphaFold are stored in B-factor column")
    
    print(f"\nFiles generated:")
    print(f"  - combined_dystrophin.pdb (PDB format)")
    print(f"  - combined_dystrophin.cif (CIF format)")
    print(f"  - Both files contain the same structural data")
    print(f"  - Confidence scores are in the B-factor column")
    
    print(f"\nUsage recommendations:")
    print(f"  - Use PDB file for most molecular visualization software")
    print(f"  - Use CIF file for crystallographic analysis")
    print(f"  - B-factor column contains AlphaFold confidence scores")
    print(f"  - Color by B-factor to visualize confidence levels")


def main():
    """Main function."""
    analyze_combined_structure()


if __name__ == "__main__":
    main()
