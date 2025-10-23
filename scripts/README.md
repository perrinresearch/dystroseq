# exonSeq Scripts

This directory contains organized Python scripts for various analysis tasks related to the exonSeq package.

## Directory Structure

### `structure_analysis/`
Scripts for protein structure analysis and manipulation:
- `alignment_summary.py` - Summary of structural alignments
- `create_clean_structure.py` - Create cleaned protein structures
- `create_final_structure.py` - Generate final combined structures
- `final_structure_summary.py` - Summary of final structure analysis
- `fix_pdb_columns.py` - Fix PDB file column formatting
- `fix_pdb_format.py` - Fix PDB file format issues
- `fix_structure_format.py` - Fix general structure format issues
- `renumbering_summary.py` - Summary of residue renumbering
- `structure_combination_report.py` - Report on structure combination process
- `test_structure_parsing.py` - Test structure file parsing

### `sequence_analysis/`
Scripts for sequence analysis and comparison:
- `analyze_dystrophin.py` - Analyze dystrophin sequences
- `analyze_mismatches.py` - Analyze sequence mismatches
- `analyze_missing_aa.py` - Analyze missing amino acids
- `compare_analysis.py` - Compare different analyses
- `create_dystrophin_fasta.py` - Create dystrophin FASTA files

### `utilities/`
Debug and utility scripts:
- `debug_codon_assignment.py` - Debug codon assignment issues
- `debug_sequence_issue.py` - Debug sequence-related issues
- `debug_simple.py` - Simple debugging utilities

## Usage

These scripts are designed to work with the exonSeq package and can be run independently for specific analysis tasks. Each script should be run from the project root directory to ensure proper import paths.

## Dependencies

All scripts require the exonSeq package and its dependencies as specified in `requirements.txt`.
