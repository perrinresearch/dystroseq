"""
CLI interface for exon utilities.
"""

import argparse
import sys
from pathlib import Path
from typing import List

from .exon_utils import (
    parse_exon_report, convert_exon_report_to_json, 
    generate_fasta_from_exon_range, generate_fasta_from_exon_list,
    get_exon_info, list_exons, generate_formatted_protein_fasta
)
from .report_generator import generate_reproducibility_report


def parse_exon_args() -> argparse.Namespace:
    """Parse command line arguments for exon utilities."""
    parser = argparse.ArgumentParser(description="exonSeq utilities for working with exon reports")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Convert command
    convert_parser = subparsers.add_parser('convert', help='Convert exon report to JSON')
    convert_parser.add_argument('input', help='Input exon report file (.txt)')
    convert_parser.add_argument('output', help='Output JSON file')
    
    # Generate FASTA from range command
    fasta_range_parser = subparsers.add_parser('fasta-range', help='Generate FASTA from exon range')
    fasta_range_parser.add_argument('exon_report', help='Exon report file (.txt or .json)')
    fasta_range_parser.add_argument('start', type=int, help='Starting exon number')
    fasta_range_parser.add_argument('end', type=int, help='Ending exon number')
    fasta_range_parser.add_argument('output', help='Output FASTA file')
    fasta_range_parser.add_argument('--type', choices=['dna', 'rna', 'protein'], 
                                  default='dna', help='Sequence type (default: dna)')
    fasta_range_parser.add_argument('--prefix', default='exon_range', 
                                  help='Header prefix (default: exon_range)')
    
    # Generate FASTA from list command
    fasta_list_parser = subparsers.add_parser('fasta-list', help='Generate FASTA from exon list')
    fasta_list_parser.add_argument('exon_report', help='Exon report file (.txt or .json)')
    fasta_list_parser.add_argument('exons', nargs='+', type=int, help='Exon numbers to include')
    fasta_list_parser.add_argument('output', help='Output FASTA file')
    fasta_list_parser.add_argument('--type', choices=['dna', 'rna', 'protein'], 
                                  default='dna', help='Sequence type (default: dna)')
    fasta_list_parser.add_argument('--prefix', default='exon_list', 
                                  help='Header prefix (default: exon_list)')
    
    # List exons command
    list_parser = subparsers.add_parser('list', help='List available exons')
    list_parser.add_argument('exon_report', help='Exon report file (.txt or .json)')
    
    # Info command
    info_parser = subparsers.add_parser('info', help='Get information about a specific exon')
    info_parser.add_argument('exon_report', help='Exon report file (.txt or .json)')
    info_parser.add_argument('exon_number', type=int, help='Exon number to look up')
    
    # Formatted protein command
    formatted_parser = subparsers.add_parser('formatted-protein', help='Generate formatted protein FASTA with spaces between exons')
    formatted_parser.add_argument('exon_report', help='Exon report file (.txt or .json)')
    formatted_parser.add_argument('output', help='Output FASTA file')
    formatted_parser.add_argument('--line-length', type=int, default=120, help='Amino acids per line (default: 120)')
    formatted_parser.add_argument('--prefix', default='formatted_protein', help='Header prefix (default: formatted_protein)')
    
    return parser.parse_args()


def main() -> None:
    """Main entry point for exon utilities CLI."""
    args = parse_exon_args()
    
    if not args.command:
        print("Error: No command specified. Use --help for usage information.")
        sys.exit(1)
    
    try:
        if args.command == 'convert':
            convert_exon_report_to_json(args.input, args.output)
            print(f"Converted {args.input} to {args.output}")
            
        elif args.command == 'fasta-range':
            exon_data = parse_exon_report(args.exon_report)
            generate_fasta_from_exon_range(
                exon_data, args.start, args.end, args.output, 
                args.type, args.prefix
            )
            print(f"Generated FASTA file: {args.output}")
            print(f"Exons {args.start}-{args.end}, type: {args.type}")
            
        elif args.command == 'fasta-list':
            exon_data = parse_exon_report(args.exon_report)
            generate_fasta_from_exon_list(
                exon_data, args.exons, args.output, 
                args.type, args.prefix
            )
            print(f"Generated FASTA file: {args.output}")
            print(f"Exons {args.exons}, type: {args.type}")
            
        elif args.command == 'list':
            exon_data = parse_exon_report(args.exon_report)
            available_exons = list_exons(exon_data)
            print(f"Available exons: {available_exons}")
            print(f"Total exons: {len(available_exons)}")
            if available_exons:
                print(f"Range: {min(available_exons)}-{max(available_exons)}")
            
        elif args.command == 'info':
            exon_data = parse_exon_report(args.exon_report)
            exon_info = get_exon_info(exon_data, args.exon_number)
            if exon_info:
                print(f"Exon {args.exon_number} information:")
                for key, value in exon_info.items():
                    print(f"  {key}: {value}")
            else:
                print(f"Exon {args.exon_number} not found")
                sys.exit(1)
                
        elif args.command == 'formatted-protein':
            exon_data = parse_exon_report(args.exon_report)
            generate_formatted_protein_fasta(
                exon_data, args.output, args.line_length, args.prefix
            )
            print(f"Generated formatted protein FASTA: {args.output}")
            print(f"Line length: {args.line_length} amino acids")
            print("Spaces added between exons for readability")
        
        # Generate reproducibility report for exon utilities
        try:
            command_args = {
                "exon_command": args.command,
                "exon_args": sys.argv[1:]  # All arguments passed to exon utilities
            }
            # Use the output directory from the last command or current directory
            output_dir = Path(".")
            if hasattr(args, 'output') and args.output:
                output_dir = Path(args.output).parent
            elif hasattr(args, 'exon_report') and args.exon_report:
                output_dir = Path(args.exon_report).parent
            
            report_path = generate_reproducibility_report(command_args, output_dir, "exon")
            print(f"\nReproducibility report generated: {report_path}")
        except Exception as report_error:
            print(f"Warning: Could not generate reproducibility report: {report_error}")
                
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
