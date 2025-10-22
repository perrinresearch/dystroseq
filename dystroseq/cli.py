import argparse
import os
import sys
from typing import Optional

from .ensembl_client import EnsemblClient
from .seq_utils import (
	to_rna, translate_dna_cds, pick_transcript, transcript_exons_in_order, 
	build_cdna_exon_ranges, exon_protein_segments, parse_variant_notation,
	get_variant_affected_exons, extract_variant_sequence, create_variant_report,
	generate_variant_protein_sequence
)


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Export DNA, RNA, and protein sequences with exon mapping from Ensembl")
	p.add_argument("--species", required=True, help="Ensembl species name, e.g., homo_sapiens")
	p.add_argument("--gene", required=True, help="Gene symbol, e.g., DMD")
	p.add_argument("--reference", required=False, default=None, help="Reference genome coord system version, e.g., GRCh38")
	p.add_argument("--coords", required=False, default=None, help="Genomic region CHR:START-END:STRAND for DNA export")
	p.add_argument("--transcript", required=False, default=None, help="Optional Ensembl transcript ID override")
	p.add_argument("--variant", required=False, default=None, help="Variant notation (e.g., chrX:g.123456-789012del, c.123-456_789+012dup)")
	p.add_argument("--out", required=True, help="Output directory")
	return p.parse_args()


def write_text(path: str, text: str) -> None:
	with open(path, "w", encoding="utf-8") as f:
		f.write(text)


def write_fasta(path: str, header: str, sequence: str, width: int = 60) -> None:
	lines = [f">{header}"]
	for i in range(0, len(sequence), width):
		lines.append(sequence[i:i + width])
	write_text(path, "\n".join(lines) + "\n")


def main() -> None:
	args = parse_args()
	os.makedirs(args.out, exist_ok=True)

	client = EnsemblClient()
	gene_obj = client.lookup_gene_by_symbol(args.species, args.gene)
	gene_id = gene_obj["id"]

	if args.transcript:
		tr = client.lookup_id(args.transcript, expand=True)
	else:
		trs = client.gene_transcripts(gene_id)
		tr = pick_transcript(trs)

	transcript_id = tr["id"]

	# Fetch sequences
	cdna_seq = client.sequence_id(transcript_id, "cdna")
	cds_seq = client.sequence_id(transcript_id, "cds")
	protein_seq = client.sequence_id(transcript_id, "protein")

	# Region DNA if provided
	if args.coords:
		region_dna = client.sequence_region(args.species, args.coords, coord_system_version=args.reference)
		write_fasta(os.path.join(args.out, "dna_region.fasta"), f"{args.species}|{args.coords}", region_dna)

	# Write transcript sequences
	write_fasta(os.path.join(args.out, "transcript_cdna.fasta"), f"{transcript_id}|cDNA", cdna_seq)
	write_fasta(os.path.join(args.out, "transcript_cds.fasta"), f"{transcript_id}|CDS", cds_seq)
	write_fasta(os.path.join(args.out, "protein.fasta"), f"{transcript_id}|protein", protein_seq)

	# Per-exon mapping
	exons_ordered = transcript_exons_in_order(tr)
	
	# Add exon rank information to preserve original exon numbering
	# The original exons list from Ensembl is in the correct order
	original_exons = tr.get("Exon", [])
	for i, exon in enumerate(exons_ordered):
		# Find the original index of this exon in the Ensembl data
		original_index = None
		for j, orig_exon in enumerate(original_exons):
			if orig_exon["start"] == exon["start"] and orig_exon["end"] == exon["end"]:
				original_index = j + 1  # Convert to 1-based
				break
		exon["rank"] = original_index if original_index else i + 1
	
	exon_cdna_ranges = build_cdna_exon_ranges(exons_ordered)
	protein_segments = exon_protein_segments(cdna_seq, cds_seq, exon_cdna_ranges, protein_seq)

	# For exon DNA and RNA sequences, fetch directly per exon region
	exon_lines = []
	for i, exon in enumerate(exons_ordered, start=1):
		chrom = gene_obj.get("seq_region_name") or tr.get("seq_region_name")
		strand = exon.get("strand", tr.get("strand", 1))
		region = f"{chrom}:{exon['start']}-{exon['end']}:{strand}"
		dna = client.sequence_region(args.species, region, coord_system_version=args.reference)
		rna = to_rna(dna)
		prot_info = protein_segments[i - 1]
		exon_lines.append(
			"\t".join([
				str(i),
				str(exon["start"]),
				str(exon["end"]),
				str(strand),
				dna,
				rna,
				str(prot_info.get("protein_start_index") or ""),
				str(prot_info.get("protein_end_index") or ""),
				prot_info.get("protein_seq") or "",
			])
		)

	headers = [
		"exon_number",
		"genomic_start",
		"genomic_end",
		"strand",
		"dna",
		"rna",
		"protein_start_index",
		"protein_end_index",
		"protein_seq",
	]
	write_text(os.path.join(args.out, "exon_report.txt"), "\t".join(headers) + "\n" + "\n".join(exon_lines) + "\n")

	# Handle variant analysis if provided
	if args.variant:
		try:
			variant = parse_variant_notation(args.variant)
			print(f"Parsed variant: {variant}")
			
			# If chromosome not specified, use the gene's chromosome
			if variant.chrom == "unknown":
				variant.chrom = gene_obj.get("seq_region_name") or tr.get("seq_region_name")
				print(f"Using gene chromosome: {variant.chrom}")
			
			# Get affected exons
			affected_exons = get_variant_affected_exons(variant, exons_ordered)
			print(f"Found {len(affected_exons)} affected exons")
			
			# Extract variant sequence
			variant_seq = extract_variant_sequence(client, args.species, variant, args.reference)
			
			# Get protein segments for affected exons
			affected_protein_segments = []
			for exon in affected_exons:
				exon_num = exons_ordered.index(exon) + 1
				if exon_num <= len(protein_segments):
					affected_protein_segments.append(protein_segments[exon_num - 1])
			
			# Create variant report
			variant_report = create_variant_report(variant, affected_exons, variant_seq, affected_protein_segments, exons_ordered)
			write_text(os.path.join(args.out, "variant_report.txt"), variant_report)
			
			# Write variant sequence to FASTA
			write_fasta(os.path.join(args.out, "variant_sequence.fasta"), 
					   f"{variant.variant_type}_{variant.chrom}_{variant.start}_{variant.end}", 
					   variant_seq)
			
			# Generate variant protein sequence
			variant_protein_seq = generate_variant_protein_sequence(variant, protein_segments, affected_exons, exons_ordered)
			
			# Write variant protein sequence to FASTA
			write_fasta(os.path.join(args.out, "variant_protein.fasta"), 
					   f"{variant.variant_type}_{variant.chrom}_{variant.start}_{variant.end}_protein", 
					   variant_protein_seq)
			
			print(f"Variant analysis complete. Check {args.out}/variant_report.txt, {args.out}/variant_sequence.fasta, and {args.out}/variant_protein.fasta")
			
		except Exception as e:
			print(f"Error processing variant: {e}")
			sys.exit(1)  # Exit with error code for variant processing failures


if __name__ == "__main__":
	main()
