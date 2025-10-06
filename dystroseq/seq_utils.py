import re
from typing import Dict, List, Tuple, Optional, Union

from Bio.Seq import Seq


def to_rna(dna: str) -> str:
	return dna.upper().replace("T", "U")


def translate_dna_cds(dna_cds: str) -> str:
	return str(Seq(dna_cds).translate(to_stop=False))


def pick_transcript(transcripts: List[Dict], prefer_canonical: bool = True) -> Dict:
	if not transcripts:
		raise ValueError("No transcripts found for gene")
	if prefer_canonical:
		for t in transcripts:
			if t.get("is_canonical") or t.get("is_canonical") == 1 or t.get("canonical_transcript"):
				return t
	# fallback: longest CDS length
	def cds_len(t: Dict) -> int:
		try:
			translation = t.get("Translation") or t.get("translation") or {}
			length = int(translation.get("length")) if translation and translation.get("length") else 0
			return length
		except Exception:
			return 0
	return sorted(transcripts, key=cds_len, reverse=True)[0]


def transcript_exons_in_order(transcript: Dict) -> List[Dict]:
	exons = transcript.get("Exon") or transcript.get("exons") or []
	if not exons:
		raise ValueError("Transcript lacks exon annotations (need expand=1 lookup)")
	strand = transcript.get("strand", 1)
	if strand == 1:
		return sorted(exons, key=lambda e: (e["start"], e["end"]))
	return sorted(exons, key=lambda e: (-(e["end"]), -(e["start"])))


def build_cdna_exon_ranges(exons_in_order: List[Dict]) -> List[Tuple[int, int, Dict]]:
	ranges: List[Tuple[int, int, Dict]] = []
	pos = 0
	for exon in exons_in_order:
		length_nt = int(exon["end"]) - int(exon["start"]) + 1
		start = pos
		end = pos + length_nt
		ranges.append((start, end, exon))
		pos = end
	return ranges


def _normalize(seq: str) -> str:
	# Remove FASTA header if present (lines starting with >)
	lines = seq.split('\n')
	seq_lines = [line for line in lines if not line.startswith('>')]
	seq_only = ''.join(seq_lines)
	
	# Remove common Ensembl sequence prefixes that appear in the sequence
	prefixes_to_remove = [
		"VERSIONQUERYENSTDESCNULLMOLECULEDNAIDENSTSEQ",
		"QUERYENSTSEQ",
		"MOLECULEDNASEQ",
		"SEQ"
	]
	
	cleaned = seq_only.upper()
	for prefix in prefixes_to_remove:
		if cleaned.startswith(prefix):
			cleaned = cleaned[len(prefix):]
			break
	
	return "".join(ch for ch in cleaned if ch.isalpha())


def find_cds_in_cdna(cdna_seq: str, cds_seq: str) -> Tuple[int, int]:
	cdna_norm = _normalize(cdna_seq)
	cds_norm = _normalize(cds_seq)
	
	# Find the start codon (ATG) in both sequences
	cds_atg_pos = cds_norm.find("ATG")
	if cds_atg_pos == -1:
		raise ValueError("No start codon (ATG) found in CDS sequence")
	
	# Extract CDS from start codon onwards
	cds_from_start = cds_norm[cds_atg_pos:]
	
	# Find this CDS sequence in the cDNA
	idx = cdna_norm.find(cds_from_start)
	if idx == -1:
		raise ValueError(f"CDS not found within cDNA sequence; transcript mismatch?\nCDNA length: {len(cdna_norm)}, CDS length: {len(cds_norm)}")
	
	return idx, idx + len(cds_from_start)


def exon_protein_segments(
	cdna_seq: str,
	cds_seq: str,
	exon_cdna_ranges: List[Tuple[int, int, Dict]],
	protein_seq: str,
) -> List[Dict]:
	cds_start, cds_end = find_cds_in_cdna(cdna_seq, cds_seq)
	results: List[Dict] = []
	for i, (ex_start, ex_end, exon) in enumerate(exon_cdna_ranges, start=1):
		inter_start = max(ex_start, cds_start)
		inter_end = min(ex_end, cds_end)
		nt_len = max(0, inter_end - inter_start)
		protein_subseq = ""
		protein_start_index: Optional[int] = None
		protein_end_index: Optional[int] = None
		if nt_len >= 3:
			first_codon_index = (inter_start - cds_start) // 3
			full_codons = nt_len // 3
			protein_start_index = first_codon_index + 1
			protein_end_index = first_codon_index + full_codons
			protein_subseq = _normalize(protein_seq)[first_codon_index:first_codon_index + full_codons]
		results.append({
			"exon_number": i,
			"exon_genomic_start": int(exon["start"]),
			"exon_genomic_end": int(exon["end"]),
			"strand": int(exon.get("strand", 1)),
			"protein_start_index": protein_start_index,
			"protein_end_index": protein_end_index,
			"protein_seq": protein_subseq,
		})
	return results


# Variant notation parsing and handling
class VariantNotation:
	"""Represents a genetic variant notation"""
	def __init__(self, variant_type: str, chrom: str, start: int, end: int, strand: int = 1):
		self.variant_type = variant_type  # "del", "dup", "ins", "inv"
		self.chrom = chrom
		self.start = start
		self.end = end
		self.strand = strand
	
	def __repr__(self):
		return f"{self.variant_type}({self.chrom}:{self.start}-{self.end})"


def parse_variant_notation(notation: str) -> VariantNotation:
	"""
	Parse common genetic variant notations:
	- c.123-456_789+012del (deletion)
	- c.123-456_789+012dup (duplication)
	- g.123456_789012del (genomic deletion)
	- g.123456_789012dup (genomic duplication)
	- chrX:g.123456-789012del
	- chrX:g.123456-789012dup
	"""
	notation = notation.strip().upper()
	
	# Handle genomic coordinates with chromosome prefix
	genomic_match = re.match(r'CHR?(\w+):G\.(\d+)-(\d+)(DEL|DUP)', notation)
	if genomic_match:
		chrom, start, end, var_type = genomic_match.groups()
		return VariantNotation(var_type.lower(), chrom, int(start), int(end))
	
	# Handle genomic coordinates without chromosome prefix
	genomic_simple = re.match(r'G\.(\d+)_(\d+)(DEL|DUP)', notation)
	if genomic_simple:
		start, end, var_type = genomic_simple.groups()
		return VariantNotation(var_type.lower(), "unknown", int(start), int(end))
	
	# Handle cDNA coordinates (c. notation)
	cdna_match = re.match(r'C\.(\d+)([+-]\d+)?_(\d+)([+-]\d+)?(DEL|DUP)', notation)
	if cdna_match:
		start_base, start_offset, end_base, end_offset, var_type = cdna_match.groups()
		# For cDNA coordinates, we'll need to map to genomic coordinates later
		# For now, return a placeholder that indicates cDNA coordinates
		start_offset = int(start_offset) if start_offset else 0
		end_offset = int(end_offset) if end_offset else 0
		return VariantNotation(
			var_type.lower(), 
			"cDNA", 
			int(start_base) + start_offset, 
			int(end_base) + end_offset
		)
	
	# Handle simple genomic ranges
	range_match = re.match(r'(\d+)-(\d+)(DEL|DUP)', notation)
	if range_match:
		start, end, var_type = range_match.groups()
		return VariantNotation(var_type.lower(), "unknown", int(start), int(end))
	
	raise ValueError(f"Unable to parse variant notation: {notation}")


def get_variant_affected_exons(variant: VariantNotation, exons: List[Dict]) -> List[Dict]:
	"""Return exons that are affected by the variant"""
	affected = []
	
	for exon in exons:
		exon_start = int(exon["start"])
		exon_end = int(exon["end"])
		
		# Check if variant overlaps with exon
		if (variant.start <= exon_end and variant.end >= exon_start):
			affected.append(exon)
	
	return affected


def extract_variant_sequence(client, species: str, variant: VariantNotation, coord_system_version: Optional[str] = None) -> str:
	"""Extract sequence for a variant region"""
	if variant.chrom == "cDNA":
		raise ValueError("cDNA coordinates not yet supported for direct sequence extraction")
	
	region = f"{variant.chrom}:{variant.start}-{variant.end}:{variant.strand}"
	return client.sequence_region(species, region, coord_system_version)


def create_variant_report(variant: VariantNotation, affected_exons: List[Dict], 
                         variant_seq: str, protein_segments: List[Dict]) -> str:
	"""Create a detailed report for the variant"""
	lines = []
	lines.append(f"Variant: {variant.variant_type.upper()}({variant.chrom}:{variant.start}-{variant.end})")
	lines.append(f"Variant sequence length: {len(variant_seq)} bp")
	lines.append(f"Number of affected exons: {len(affected_exons)}")
	lines.append("")
	
	lines.append("Affected exons:")
	for i, exon in enumerate(affected_exons, 1):
		lines.append(f"  Exon {i}: {exon['start']}-{exon['end']} (strand {exon.get('strand', 1)})")
	
	lines.append("")
	lines.append("Variant sequence:")
	# Format sequence in 60-character lines
	for i in range(0, len(variant_seq), 60):
		lines.append(variant_seq[i:i+60])
	
	lines.append("")
	lines.append("Protein segments affected:")
	for segment in protein_segments:
		if segment.get("protein_seq"):
			lines.append(f"  Exon {segment['exon_number']}: {segment['protein_seq']}")
	
	return "\n".join(lines)
