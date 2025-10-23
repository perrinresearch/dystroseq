"""
Utilities for analyzing AlphaFold confidence scores from DystroSeq structures.
"""

import json
import os
import glob
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np


def parse_alphafold_confidence_data(json_file: str) -> Dict:
    """
    Parse AlphaFold confidence data from a JSON file.
    
    Args:
        json_file: Path to the AlphaFold full_data JSON file
        
    Returns:
        Dictionary containing confidence data
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    return {
        'atom_plddts': data.get('atom_plddts', []),
        'pae': data.get('pae', []),
        'ptm': data.get('ptm', None),
        'iptm': data.get('iptm', None)
    }


def get_confidence_scores_for_structure(structure_dir: str) -> List[Dict]:
    """
    Extract confidence scores from all models in a structure directory.
    
    Args:
        structure_dir: Path to the structure directory (e.g., dystroseq_exons_1_9_protein)
        
    Returns:
        List of confidence data dictionaries for each model
    """
    confidence_data = []
    
    # Find all full_data JSON files
    pattern = os.path.join(structure_dir, "fold_*_full_data_*.json")
    json_files = glob.glob(pattern)
    
    for json_file in sorted(json_files):
        try:
            data = parse_alphafold_confidence_data(json_file)
            confidence_data.append(data)
        except Exception as e:
            print(f"Warning: Could not parse {json_file}: {e}")
    
    return confidence_data


def extract_exon_range_from_dirname(dirname: str) -> Tuple[int, int]:
    """
    Extract exon range from directory name.
    
    Args:
        dirname: Directory name like 'dystroseq_exons_1_9_protein'
        
    Returns:
        Tuple of (start_exon, end_exon)
    """
    # Extract the exon range from the directory name
    parts = dirname.split('_')
    if len(parts) >= 4 and parts[1] == 'exons':
        try:
            start_exon = int(parts[2])
            end_exon = int(parts[3])
            return start_exon, end_exon
        except ValueError:
            pass
    
    # Fallback: try to extract numbers from the string
    import re
    numbers = re.findall(r'\d+', dirname)
    if len(numbers) >= 2:
        return int(numbers[0]), int(numbers[1])
    
    return None, None


def calculate_average_confidence_per_residue(confidence_data_list: List[Dict]) -> List[float]:
    """
    Calculate average confidence (pLDDT) per residue across all models.
    
    Args:
        confidence_data_list: List of confidence data dictionaries
        
    Returns:
        List of average confidence scores per residue
    """
    if not confidence_data_list:
        return []
    
    # Get the length from the first model
    first_model_length = len(confidence_data_list[0]['atom_plddts'])
    
    # Initialize arrays for averaging
    total_confidence = np.zeros(first_model_length)
    model_count = 0
    
    for data in confidence_data_list:
        plddts = data['atom_plddts']
        if len(plddts) == first_model_length:
            total_confidence += np.array(plddts)
            model_count += 1
    
    if model_count == 0:
        return []
    
    # Calculate average
    average_confidence = total_confidence / model_count
    return average_confidence.tolist()


def analyze_all_structures(base_dir: str) -> pd.DataFrame:
    """
    Analyze confidence scores across all AlphaFold structures.
    
    Args:
        base_dir: Base directory containing all structure folders
        
    Returns:
        DataFrame with confidence analysis results
    """
    results = []
    
    # Find all structure directories
    structure_dirs = [d for d in os.listdir(base_dir) 
                     if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('dystroseq_exons_')]
    
    for structure_dir in sorted(structure_dirs):
        print(f"Processing {structure_dir}...")
        
        # Extract exon range
        start_exon, end_exon = extract_exon_range_from_dirname(structure_dir)
        if start_exon is None or end_exon is None:
            print(f"Warning: Could not extract exon range from {structure_dir}")
            continue
        
        # Get confidence data
        full_path = os.path.join(base_dir, structure_dir)
        confidence_data_list = get_confidence_scores_for_structure(full_path)
        
        if not confidence_data_list:
            print(f"Warning: No confidence data found for {structure_dir}")
            continue
        
        # Calculate average confidence per residue
        avg_confidence = calculate_average_confidence_per_residue(confidence_data_list)
        
        if avg_confidence:
            # Calculate statistics
            mean_confidence = np.mean(avg_confidence)
            std_confidence = np.std(avg_confidence)
            min_confidence = np.min(avg_confidence)
            max_confidence = np.max(avg_confidence)
            
            # Count residues by confidence level
            high_conf = sum(1 for c in avg_confidence if c >= 70)
            medium_conf = sum(1 for c in avg_confidence if 50 <= c < 70)
            low_conf = sum(1 for c in avg_confidence if c < 50)
            
            results.append({
                'structure': structure_dir,
                'start_exon': start_exon,
                'end_exon': end_exon,
                'num_models': len(confidence_data_list),
                'num_residues': len(avg_confidence),
                'mean_confidence': mean_confidence,
                'std_confidence': std_confidence,
                'min_confidence': min_confidence,
                'max_confidence': max_confidence,
                'high_confidence_residues': high_conf,
                'medium_confidence_residues': medium_conf,
                'low_confidence_residues': low_conf,
                'high_confidence_percent': (high_conf / len(avg_confidence)) * 100,
                'medium_confidence_percent': (medium_conf / len(avg_confidence)) * 100,
                'low_confidence_percent': (low_conf / len(avg_confidence)) * 100
            })
            
            print(f"  Exons {start_exon}-{end_exon}: {len(avg_confidence)} residues, "
                  f"mean confidence: {mean_confidence:.1f}, "
                  f"high confidence: {high_conf} ({high_conf/len(avg_confidence)*100:.1f}%)")
    
    return pd.DataFrame(results)


def create_residue_confidence_table(base_dir: str, output_file: str) -> None:
    """
    Create a detailed table of confidence scores for each residue across all structures.
    
    Args:
        base_dir: Base directory containing all structure folders
        output_file: Output CSV file path
    """
    all_residue_data = []
    
    # Find all structure directories
    structure_dirs = [d for d in os.listdir(base_dir) 
                     if os.path.isdir(os.path.join(base_dir, d)) and d.startswith('dystroseq_exons_')]
    
    for structure_dir in sorted(structure_dirs):
        print(f"Processing {structure_dir}...")
        
        # Extract exon range
        start_exon, end_exon = extract_exon_range_from_dirname(structure_dir)
        if start_exon is None or end_exon is None:
            continue
        
        # Get confidence data
        full_path = os.path.join(base_dir, structure_dir)
        confidence_data_list = get_confidence_scores_for_structure(full_path)
        
        if not confidence_data_list:
            continue
        
        # Calculate average confidence per residue
        avg_confidence = calculate_average_confidence_per_residue(confidence_data_list)
        
        # Add residue-level data
        for i, confidence in enumerate(avg_confidence):
            all_residue_data.append({
                'structure': structure_dir,
                'start_exon': start_exon,
                'end_exon': end_exon,
                'residue_index': i + 1,  # 1-based indexing
                'confidence': confidence,
                'confidence_category': 'high' if confidence >= 70 else 'medium' if confidence >= 50 else 'low'
            })
    
    # Create DataFrame and save
    df = pd.DataFrame(all_residue_data)
    df.to_csv(output_file, index=False)
    print(f"Residue confidence table saved to {output_file}")


def main():
    """Main function to analyze AlphaFold confidence scores."""
    base_dir = "../../../Data/dystrophin/folds_2025_10_22_19_55"
    
    print("Analyzing AlphaFold confidence scores...")
    
    # Analyze all structures
    results_df = analyze_all_structures(base_dir)
    
    if not results_df.empty:
        # Save summary results
        summary_file = "alphafold_confidence_summary.csv"
        results_df.to_csv(summary_file, index=False)
        print(f"\nSummary results saved to {summary_file}")
        
        # Print summary
        print("\nSummary of results:")
        print(results_df[['structure', 'start_exon', 'end_exon', 'num_residues', 
                          'mean_confidence', 'high_confidence_percent']].to_string(index=False))
        
        # Create detailed residue table
        residue_file = "alphafold_residue_confidence.csv"
        create_residue_confidence_table(base_dir, residue_file)
        
        # Overall statistics
        print(f"\nOverall statistics:")
        print(f"Total structures analyzed: {len(results_df)}")
        print(f"Total residues: {results_df['num_residues'].sum()}")
        print(f"Overall mean confidence: {results_df['mean_confidence'].mean():.1f}")
        print(f"Overall high confidence percentage: {results_df['high_confidence_percent'].mean():.1f}%")
    
    else:
        print("No structures found or processed.")


if __name__ == "__main__":
    main()
