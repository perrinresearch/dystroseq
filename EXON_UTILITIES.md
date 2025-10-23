# exonSeq Exon Utilities

This document describes the new exon utilities added to exonSeq for working with exon reports and generating FASTA files from exon ranges.

## Overview

The exon utilities provide functionality to:
1. Convert exon report text files to JSON format
2. Generate FASTA files from ranges of exons
3. Generate FASTA files from specific lists of exons
4. List available exons in a report
5. Get detailed information about specific exons

## Usage

All exon utilities are accessed through the main DystroSeq CLI using the `exon` command:

```bash
python -m exonseq exon <command> [arguments]
```

## Commands

### 1. Convert Exon Report to JSON

Convert a tab-separated exon report file to JSON format:

```bash
python -m exonseq exon convert <input_file> <output_file>
```

**Example:**
```bash
python -m exonseq exon convert dystrophin_analysis/exon_report.txt exon_report.json
```

### 2. Generate FASTA from Exon Range

Generate a FASTA file containing sequences from a range of exons:

```bash
python -m exonseq exon fasta-range <exon_report> <start_exon> <end_exon> <output_file> [--type <sequence_type>] [--prefix <header_prefix>]
```

**Parameters:**
- `exon_report`: Path to exon report file (.txt or .json)
- `start_exon`: Starting exon number (1-based)
- `end_exon`: Ending exon number (1-based, inclusive)
- `output_file`: Path to output FASTA file
- `--type`: Sequence type - `dna`, `rna`, or `protein` (default: `dna`)
- `--prefix`: Header prefix for FASTA file (default: `exon_range`)

**Examples:**
```bash
# Generate DNA sequence for exons 1-5
python -m exonseq exon fasta-range exon_report.txt 1 5 exons_1-5_dna.fasta

# Generate protein sequence for exons 10-20
python -m exonseq exon fasta-range exon_report.txt 10 20 exons_10-20_protein.fasta --type protein

# Generate RNA sequence with custom prefix
python -m exonseq exon fasta-range exon_report.txt 1 3 exons_1-3_rna.fasta --type rna --prefix custom
```

### 3. Generate FASTA from Exon List

Generate a FASTA file containing sequences from a specific list of exons:

```bash
python -m exonseq exon fasta-list <exon_report> <exon_numbers...> <output_file> [--type <sequence_type>] [--prefix <header_prefix>]
```

**Parameters:**
- `exon_report`: Path to exon report file (.txt or .json)
- `exon_numbers`: Space-separated list of exon numbers
- `output_file`: Path to output FASTA file
- `--type`: Sequence type - `dna`, `rna`, or `protein` (default: `dna`)
- `--prefix`: Header prefix for FASTA file (default: `exon_list`)

**Examples:**
```bash
# Generate DNA sequence for exons 1, 3, 5, 7
python -m exonseq exon fasta-list exon_report.txt 1 3 5 7 exons_selected_dna.fasta

# Generate protein sequence for specific exons
python -m exonseq exon fasta-list exon_report.txt 10 15 20 exons_selected_protein.fasta --type protein
```

### 4. List Available Exons

Display all available exon numbers in a report:

```bash
python -m exonseq exon list <exon_report>
```

**Example:**
```bash
python -m exonseq exon list exon_report.txt
```

**Output:**
```
Available exons: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...]
Total exons: 79
Range: 1-79
```

### 5. Get Exon Information

Display detailed information about a specific exon:

```bash
python -m exonseq exon info <exon_report> <exon_number>
```

### 6. Generate Formatted Protein FASTA

Generate a formatted protein FASTA file with spaces between exons for better readability:

```bash
python -m exonseq exon formatted-protein <exon_report> <output_file> [--line-length <length>] [--prefix <prefix>]
```

**Parameters:**
- `exon_report`: Path to exon report file (.txt or .json)
- `output_file`: Path to output FASTA file
- `--line-length`: Number of amino acids per line (default: 120)
- `--prefix`: Header prefix for FASTA file (default: formatted_protein)

**Examples:**
```bash
# Generate formatted protein with default 120 amino acids per line
python -m exonseq exon formatted-protein exon_report.txt formatted_protein.fasta

# Generate with 60 amino acids per line
python -m exonseq exon formatted-protein exon_report.txt formatted_protein_60aa.fasta --line-length 60

# Generate with custom header prefix
python -m exonseq exon formatted-protein exon_report.txt dystrophin_formatted.fasta --prefix dystrophin
```

**Example:**
```bash
python -m exonseq exon info exon_report.txt 1
```

**Output:**
```
Exon 1 information:
  exon_number: 1
  genomic_start: 33211282
  genomic_end: 33211549
  strand: -1
  dna: ATCAGTTACTGTGTTGACTCACTCAGTGTTGGGATCACTCACTTTCCCCCTACAGGACTCAGATCTGGGAGGCAATTACCTTCGGAGAAAAACGAATAGGAAAAACTGAAGTGTTACTTTTTTTAAAGCTGCTGAAGTTTGTTGGTTTCTCATTGTTTTTAAGCCTACTGGAGCAATAAAGTTTGAAGAACTTTTACCAGGTTTTTTTTATCGCTGCCTTGATATACACTTTTCAAAATGCTTTGGTGGGAAGAAGTAGAGGACTGTT
  rna: AUCAGUUACUGUGUUGACUCACUCAGUGUUGGGAUCACUCACUUUCCCCCUACAGGACUCAGAUCUGGGAGGCAAUUACCUUCGGAGAAAAACGAAUAGGAAAAACUGAAGUGUUACUUUUUUUAAAGCUGCUGAAGUUUGUUGGUUUCUCAUUGUUUUUAAGCCUACUGGAGCAAUAAAGUUUGAAGAACUUUUACCAGGUUUUUUUUAUCGCUGCCUUGAUAUACACUUUUCAAAAUGCUUUGGUGGGAAGAAGUAGAGGACUGUU
  protein_start_index: 1
  protein_end_index: 10
  protein_seq: MLWWEEVEDC
```

## File Formats

### Input Files

The utilities accept both tab-separated text files (`.txt`) and JSON files (`.json`) as input. The text format is the standard exon report format generated by DystroSeq with the following columns:

- `exon_number`: Exon number (1-based)
- `genomic_start`: Genomic start position
- `genomic_end`: Genomic end position
- `strand`: Strand orientation (-1 or 1)
- `dna`: DNA sequence
- `rna`: RNA sequence
- `protein_start_index`: Protein start position
- `protein_end_index`: Protein end position
- `protein_seq`: Protein sequence

### Output Files

- **JSON files**: Structured data format for programmatic access
- **FASTA files**: Standard sequence format with 60-character line wrapping

## Examples

### Complete Workflow

```bash
# 1. Convert exon report to JSON
python -m exonseq exon convert dystrophin_analysis/exon_report.txt exon_report.json

# 2. List available exons
python -m exonseq exon list exon_report.json

# 3. Generate DNA sequence for exons 1-10
python -m exonseq exon fasta-range exon_report.json 1 10 exons_1-10_dna.fasta

# 4. Generate protein sequence for specific exons
python -m exonseq exon fasta-list exon_report.json 1 5 10 15 20 exons_selected_protein.fasta --type protein

# 5. Get information about exon 1
python -m exonseq exon info exon_report.json 1
```

### Use Cases

1. **Extract specific exon regions**: Use `fasta-range` to get sequences from contiguous exon ranges
2. **Create custom exon combinations**: Use `fasta-list` to combine non-contiguous exons
3. **Convert data formats**: Use `convert` to create JSON files for programmatic access
4. **Explore exon data**: Use `list` and `info` to understand available exons and their properties

## Integration with Existing DystroSeq

These utilities work seamlessly with the existing DystroSeq workflow:

1. Generate exon reports using the main DystroSeq export command
2. Use the exon utilities to process and extract specific sequences
3. Combine with variant analysis and other DystroSeq features

The utilities are designed to be flexible and work with any exon report generated by DystroSeq or compatible tools.
