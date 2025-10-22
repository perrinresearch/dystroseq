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
	
	# Create a mapping of amino acid positions to exons
	aa_to_exon = {}
	
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
			
			# Map amino acids to this exon
			for aa_pos in range(first_codon_index, first_codon_index + full_codons):
				aa_to_exon[aa_pos] = i
		
		# Handle boundary cases - check if this exon contains part of a codon
		# that spans to the next exon
		if i < len(exon_cdna_ranges):
			next_ex_start, next_ex_end, _ = exon_cdna_ranges[i]  # i is 1-based, so exon_cdna_ranges[i] is the next exon
			next_inter_start = max(next_ex_start, cds_start)
			
			# Check if there's a partial codon at the end of this exon
			remaining_nt = (inter_end - inter_start) % 3
			if remaining_nt > 0 and next_inter_start <= cds_end:
				# This exon has partial codon that continues to next exon
				partial_codon_start = inter_end - remaining_nt
				partial_codon_end = min(inter_end + (3 - remaining_nt), next_inter_start)
				
				# Calculate which amino acid this partial codon contributes to
				partial_codon_aa_index = (partial_codon_start - cds_start) // 3
				if partial_codon_aa_index < len(_normalize(protein_seq)):
					# Assign this amino acid to the exon with more nucleotides
					exon_contribution = remaining_nt
					next_exon_contribution = 3 - remaining_nt
					
					if exon_contribution >= next_exon_contribution:
						aa_to_exon[partial_codon_aa_index] = i
		
		results.append({
			"exon_number": i,
			"exon_genomic_start": int(exon["start"]),
			"exon_genomic_end": int(exon["end"]),
			"strand": int(exon.get("strand", 1)),
			"protein_start_index": protein_start_index,
			"protein_end_index": protein_end_index,
			"protein_seq": protein_subseq,
		})
	
	# Now update the results to include all amino acids assigned to each exon
	for i, result in enumerate(results, start=1):
		exon_amino_acids = []
		for aa_pos, assigned_exon in aa_to_exon.items():
			if assigned_exon == i and aa_pos < len(_normalize(protein_seq)):
				exon_amino_acids.append((aa_pos, _normalize(protein_seq)[aa_pos]))
		
		if exon_amino_acids:
			# Sort by position and create the complete protein sequence for this exon
			exon_amino_acids.sort()
			complete_protein_seq = ''.join([aa for _, aa in exon_amino_acids])
			
			# Update the result with the complete sequence
			result["protein_seq"] = complete_protein_seq
			if exon_amino_acids:
				result["protein_start_index"] = exon_amino_acids[0][0] + 1
				result["protein_end_index"] = exon_amino_acids[-1][0] + 1
	
	# Fill in missing amino acids by assigning them to the nearest exon
	# This handles the case where there are gaps between exons
	all_aa_positions = set(range(len(_normalize(protein_seq))))
	assigned_positions = set(aa_to_exon.keys())
	missing_positions = all_aa_positions - assigned_positions
	
	for missing_pos in missing_positions:
		# Find the nearest exon
		min_distance = float('inf')
		nearest_exon = 1
		
		for i, (ex_start, ex_end, exon) in enumerate(exon_cdna_ranges, start=1):
			inter_start = max(ex_start, cds_start)
			inter_end = min(ex_end, cds_end)
			
			if inter_start <= cds_end and inter_end >= cds_start:
				# Calculate distance to this exon
				exon_start_aa = (inter_start - cds_start) // 3
				exon_end_aa = (inter_end - cds_start) // 3
				
				if exon_start_aa <= missing_pos <= exon_end_aa:
					# Missing position is within this exon's range
					nearest_exon = i
					break
				else:
					# Calculate distance to exon boundaries
					dist_to_start = abs(missing_pos - exon_start_aa)
					dist_to_end = abs(missing_pos - exon_end_aa)
					min_dist_to_exon = min(dist_to_start, dist_to_end)
					
					if min_dist_to_exon < min_distance:
						min_distance = min_dist_to_exon
						nearest_exon = i
		
		# Assign the missing amino acid to the nearest exon
		aa_to_exon[missing_pos] = nearest_exon
	
	# Update results with the filled assignments
	for i, result in enumerate(results, start=1):
		exon_amino_acids = []
		for aa_pos, assigned_exon in aa_to_exon.items():
			if assigned_exon == i and aa_pos < len(_normalize(protein_seq)):
				exon_amino_acids.append((aa_pos, _normalize(protein_seq)[aa_pos]))
		
		if exon_amino_acids:
			# Sort by position and create the complete protein sequence for this exon
			exon_amino_acids.sort()
			complete_protein_seq = ''.join([aa for _, aa in exon_amino_acids])
			
			# Update the result with the complete sequence
			result["protein_seq"] = complete_protein_seq
			if exon_amino_acids:
				result["protein_start_index"] = exon_amino_acids[0][0] + 1
				result["protein_end_index"] = exon_amino_acids[-1][0] + 1
	
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
	"""Return exons that are affected by the variant with detailed overlap information"""
	affected = []
	
	for exon in exons:
		exon_start = int(exon["start"])
		exon_end = int(exon["end"])
		
		# For deletions, duplications, and insertions, check if exon is within the variant region
		if variant.variant_type in ["del", "dup", "ins"]:
			# Check different types of overlap for deletions
			if variant.variant_type == "del":
				# Exon is completely deleted if it's entirely within the deletion
				if (exon_start >= variant.start and exon_end <= variant.end):
					exon["deletion_type"] = "complete"
					affected.append(exon)
				# Exon is partially deleted if it overlaps with the deletion boundaries
				elif (variant.start <= exon_end and variant.end >= exon_start):
					exon["deletion_type"] = "partial"
					affected.append(exon)
			else:
				# For duplications and insertions, use simple overlap
				if (variant.start <= exon_end and variant.end >= exon_start):
					exon["deletion_type"] = "overlap"
					affected.append(exon)
		else:
			# For inversions and other variant types, use simple overlap
			if (variant.start <= exon_end and variant.end >= exon_start):
				exon["deletion_type"] = "overlap"
				affected.append(exon)
	
	return affected


def extract_variant_sequence(client, species: str, variant: VariantNotation, coord_system_version: Optional[str] = None) -> str:
	"""Extract sequence for a variant region"""
	if variant.chrom == "cDNA":
		raise ValueError("cDNA coordinates not yet supported for direct sequence extraction")
	
	region = f"{variant.chrom}:{variant.start}-{variant.end}:{variant.strand}"
	return client.sequence_region(species, region, coord_system_version)


def generate_variant_protein_sequence(variant: VariantNotation, protein_segments: List[Dict], 
                                     affected_exons: List[Dict], all_exons: List[Dict]) -> str:
	"""Generate the protein sequence for a variant"""
	if variant.variant_type == "del":
		# For deletions, handle different types of exon deletions
		variant_protein_parts = []
		
		for i, segment in enumerate(protein_segments):
			# Check if this exon is affected by the variant
			exon_num = i + 1  # protein_segments is 1-indexed
			exon_index = exon_num - 1  # Convert to 0-indexed for all_exons
			
			if exon_index < len(all_exons):
				exon = all_exons[exon_index]
				
				# Check if this exon is in the affected exons list
				affected_exon = None
				for aff_exon in affected_exons:
					if aff_exon["start"] == exon["start"] and aff_exon["end"] == exon["end"]:
						affected_exon = aff_exon
						break
				
				if affected_exon is None:
					# This exon is not affected, include its full protein sequence
					if segment.get("protein_seq"):
						variant_protein_parts.append(segment["protein_seq"])
				else:
					# This exon is affected by the deletion
					deletion_type = affected_exon.get("deletion_type", "complete")
					
					if deletion_type == "complete":
						# Exon is completely deleted, skip it entirely
						continue
					elif deletion_type == "partial":
						# Exon is partially deleted, need to calculate which part remains
						remaining_protein = _calculate_partial_exon_protein(
							variant, exon, segment, protein_segments, i
						)
						if remaining_protein:
							variant_protein_parts.append(remaining_protein)
		
		return "".join(variant_protein_parts)
	
	elif variant.variant_type == "dup":
		# For duplications, we would need to duplicate the affected exons
		# This is more complex and would require additional logic
		# For now, return the original protein sequence
		return "".join([segment.get("protein_seq", "") for segment in protein_segments])
	
	elif variant.variant_type == "ins":
		# For insertions, we would need to insert new sequence
		# This is complex and would require additional logic
		# For now, return the original protein sequence
		return "".join([segment.get("protein_seq", "") for segment in protein_segments])
	
	else:
		# For other variant types, return original sequence
		return "".join([segment.get("protein_seq", "") for segment in protein_segments])


def _calculate_partial_exon_protein(variant: VariantNotation, exon: Dict, segment: Dict, 
                                  protein_segments: List[Dict], segment_index: int) -> str:
	"""Calculate the remaining protein sequence for a partially deleted exon."""
	exon_start = int(exon["start"])
	exon_end = int(exon["end"])
	exon_strand = int(exon.get("strand", 1))
	
	# Get the protein sequence for this exon
	protein_seq = segment.get("protein_seq", "")
	if not protein_seq:
		return ""
	
	# For partial deletions, we need to determine which part of the exon is deleted
	# This is complex because we need to map genomic coordinates to protein coordinates
	
	# For now, implement a simplified approach:
	# If the deletion affects the beginning of the exon, remove the corresponding protein
	# If the deletion affects the end of the exon, remove the corresponding protein
	# If the deletion is in the middle, we would need more complex logic
	
	# Calculate how much of the exon is deleted
	deletion_start_in_exon = max(variant.start, exon_start) - exon_start
	deletion_end_in_exon = min(variant.end, exon_end) - exon_start
	
	# Calculate the fraction of the exon that is deleted
	exon_length = exon_end - exon_start + 1
	deleted_fraction = (deletion_end_in_exon - deletion_start_in_exon) / exon_length
	
	# For simplicity, if more than 50% of the exon is deleted, remove the entire exon
	# If less than 50% is deleted, keep the entire exon
	# This is a simplified approach - a more sophisticated implementation would
	# calculate the exact protein residues to remove based on codon boundaries
	
	if deleted_fraction > 0.5:
		# More than half the exon is deleted, remove it entirely
		return ""
	else:
		# Less than half is deleted, keep the entire exon
		# In a more sophisticated implementation, we would calculate
		# which specific amino acids to remove
		return protein_seq


def create_variant_report(variant: VariantNotation, affected_exons: List[Dict], 
                         variant_seq: str, protein_segments: List[Dict], 
                         all_exons: List[Dict] = None) -> str:
	"""Create a detailed report for the variant"""
	lines = []
	lines.append(f"Variant: {variant.variant_type.upper()}({variant.chrom}:{variant.start}-{variant.end})")
	lines.append(f"Variant sequence length: {len(variant_seq)} bp")
	lines.append(f"Number of affected exons: {len(affected_exons)}")
	lines.append("")
	
	lines.append("Affected exons:")
	for exon in affected_exons:
		# Use the rank information that was added to the exon data
		if "rank" in exon:
			exon_number = exon["rank"]
		elif all_exons:
			try:
				# Fallback to list index + 1
				exon_number = all_exons.index(exon) + 1
			except ValueError:
				exon_number = "Unknown"
		else:
			exon_number = "Unknown"
		
		# Add deletion type information if available
		deletion_type = exon.get("deletion_type", "unknown")
		if deletion_type != "unknown":
			lines.append(f"  Exon {exon_number}: {exon['start']}-{exon['end']} (strand {exon.get('strand', 1)}, {deletion_type} deletion)")
		else:
			lines.append(f"  Exon {exon_number}: {exon['start']}-{exon['end']} (strand {exon.get('strand', 1)})")
	
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
