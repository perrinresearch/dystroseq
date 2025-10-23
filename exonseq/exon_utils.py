"""
Utilities for working with exon reports and generating FASTA files from exon ranges.
"""

import json
import os
from typing import List, Dict, Any, Union, Optional
from pathlib import Path


def parse_exon_report(file_path: Union[str, Path]) -> List[Dict[str, Any]]:
    """
    Parse an exon report file (either .txt or .json) and return a list of exon dictionaries.
    
    Args:
        file_path: Path to the exon report file
        
    Returns:
        List of dictionaries containing exon information
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"Exon report file not found: {file_path}")
    
    if file_path.suffix.lower() == '.json':
        with open(file_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    # Parse tab-separated text file
    exons = []
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    if not lines:
        return exons
    
    # Parse header
    header = lines[0].strip().split('\t')
    
    # Parse data rows
    for line in lines[1:]:
        if line.strip():  # Skip empty lines
            values = line.strip().split('\t')
            if len(values) == len(header):
                exon_dict = {}
                for i, col in enumerate(header):
                    value = values[i]
                    # Convert numeric fields
                    if col in ['exon_number', 'genomic_start', 'genomic_end', 'strand', 
                              'protein_start_index', 'protein_end_index']:
                        try:
                            exon_dict[col] = int(value) if value else None
                        except ValueError:
                            exon_dict[col] = value
                    else:
                        exon_dict[col] = value
                exons.append(exon_dict)
    
    return exons


def convert_exon_report_to_json(input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """
    Convert an exon report text file to JSON format.
    
    Args:
        input_path: Path to the input exon report text file
        output_path: Path to the output JSON file
    """
    exons = parse_exon_report(input_path)
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(exons, f, indent=2, ensure_ascii=False)


def generate_fasta_from_exon_range(
    exon_data: List[Dict[str, Any]], 
    start_exon: int, 
    end_exon: int, 
    output_path: Union[str, Path],
    sequence_type: str = 'dna',
    header_prefix: str = 'exon_range'
) -> None:
    """
    Generate a FASTA file from a range of exons.
    
    Args:
        exon_data: List of exon dictionaries from parse_exon_report
        start_exon: Starting exon number (1-based)
        end_exon: Ending exon number (1-based, inclusive)
        output_path: Path to the output FASTA file
        sequence_type: Type of sequence to extract ('dna', 'rna', 'protein')
        header_prefix: Prefix for the FASTA header
    """
    if sequence_type not in ['dna', 'rna', 'protein']:
        raise ValueError("sequence_type must be 'dna', 'rna', or 'protein'")
    
    # Filter exons in the specified range
    selected_exons = []
    for exon in exon_data:
        exon_num = exon.get('exon_number')
        if exon_num and start_exon <= exon_num <= end_exon:
            selected_exons.append(exon)
    
    if not selected_exons:
        raise ValueError(f"No exons found in range {start_exon}-{end_exon}")
    
    # Sort by exon number
    selected_exons.sort(key=lambda x: x.get('exon_number', 0))
    
    # Concatenate sequences
    sequences = []
    for exon in selected_exons:
        if sequence_type == 'protein':
            seq = exon.get('protein_seq', '')
        elif sequence_type == 'rna':
            seq = exon.get('rna', '')
        else:  # dna
            seq = exon.get('dna', '')
        sequences.append(seq)
    
    concatenated_seq = ''.join(sequences)
    
    # Create FASTA header
    header = f">{header_prefix}_exons_{start_exon}-{end_exon}_{sequence_type}"
    
    # Write FASTA file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        # Write sequence in 60-character lines
        for i in range(0, len(concatenated_seq), 60):
            f.write(concatenated_seq[i:i+60] + '\n')


def generate_fasta_from_exon_list(
    exon_data: List[Dict[str, Any]], 
    exon_numbers: List[int], 
    output_path: Union[str, Path],
    sequence_type: str = 'dna',
    header_prefix: str = 'exon_list'
) -> None:
    """
    Generate a FASTA file from a specific list of exons.
    
    Args:
        exon_data: List of exon dictionaries from parse_exon_report
        exon_numbers: List of exon numbers to include (1-based)
        output_path: Path to the output FASTA file
        sequence_type: Type of sequence to extract ('dna', 'rna', 'protein')
        header_prefix: Prefix for the FASTA header
    """
    if sequence_type not in ['dna', 'rna', 'protein']:
        raise ValueError("sequence_type must be 'dna', 'rna', or 'protein'")
    
    # Create a mapping of exon numbers to exons
    exon_map = {exon.get('exon_number'): exon for exon in exon_data}
    
    # Get selected exons in the specified order
    selected_exons = []
    for exon_num in exon_numbers:
        if exon_num in exon_map:
            selected_exons.append(exon_map[exon_num])
        else:
            print(f"Warning: Exon {exon_num} not found in data")
    
    if not selected_exons:
        raise ValueError("No valid exons found in the specified list")
    
    # Concatenate sequences
    sequences = []
    for exon in selected_exons:
        if sequence_type == 'protein':
            seq = exon.get('protein_seq', '')
        elif sequence_type == 'rna':
            seq = exon.get('rna', '')
        else:  # dna
            seq = exon.get('dna', '')
        sequences.append(seq)
    
    concatenated_seq = ''.join(sequences)
    
    # Create FASTA header
    exon_list_str = ','.join(map(str, exon_numbers))
    header = f">{header_prefix}_exons_{exon_list_str}_{sequence_type}"
    
    # Write FASTA file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        # Write sequence in 60-character lines
        for i in range(0, len(concatenated_seq), 60):
            f.write(concatenated_seq[i:i+60] + '\n')


def get_exon_info(exon_data: List[Dict[str, Any]], exon_number: int) -> Optional[Dict[str, Any]]:
    """
    Get information for a specific exon.
    
    Args:
        exon_data: List of exon dictionaries from parse_exon_report
        exon_number: Exon number to look up (1-based)
        
    Returns:
        Dictionary with exon information or None if not found
    """
    for exon in exon_data:
        if exon.get('exon_number') == exon_number:
            return exon
    return None


def list_exons(exon_data: List[Dict[str, Any]]) -> List[int]:
    """
    Get a list of available exon numbers.
    
    Args:
        exon_data: List of exon dictionaries from parse_exon_report
        
    Returns:
        List of exon numbers in ascending order
    """
    exon_numbers = [exon.get('exon_number') for exon in exon_data if exon.get('exon_number')]
    return sorted(exon_numbers)


def generate_formatted_protein_fasta(
    exon_data: List[Dict[str, Any]], 
    output_path: Union[str, Path],
    line_length: int = 120,
    header_prefix: str = 'formatted_protein'
) -> None:
    """
    Generate a formatted protein FASTA file with specified line length and spaces between exons.
    
    Args:
        exon_data: List of exon dictionaries from parse_exon_report
        output_path: Path to the output FASTA file
        line_length: Number of amino acids per line (default: 120)
        header_prefix: Prefix for the FASTA header
    """
    # Sort exons by exon number
    sorted_exons = sorted(exon_data, key=lambda x: x.get('exon_number', 0))
    
    # Collect protein sequences with exon separators
    protein_sequences = []
    for exon in sorted_exons:
        protein_seq = exon.get('protein_seq', '')
        if protein_seq:
            protein_sequences.append(protein_seq)
    
    # Join with spaces between exons
    full_sequence = ' '.join(protein_sequences)
    
    # Create FASTA header
    header = f">{header_prefix}_formatted_{line_length}aa_per_line"
    
    # Write FASTA file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(header + '\n')
        
        # Write sequence with specified line length
        current_pos = 0
        while current_pos < len(full_sequence):
            # Find the end of the current line
            line_end = min(current_pos + line_length, len(full_sequence))
            line = full_sequence[current_pos:line_end]
            f.write(line + '\n')
            current_pos = line_end
