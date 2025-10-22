"""Unit tests for dystroseq.seq_utils module."""

import pytest
from dystroseq.seq_utils import (
    to_rna, translate_dna_cds, pick_transcript, transcript_exons_in_order,
    build_cdna_exon_ranges, exon_protein_segments, parse_variant_notation,
    get_variant_affected_exons, extract_variant_sequence, create_variant_report,
    generate_variant_protein_sequence, VariantNotation
)


class TestBasicFunctions:
    """Test basic sequence utility functions."""
    
    def test_to_rna(self):
        """Test DNA to RNA conversion."""
        assert to_rna("ATCG") == "AUCG"
        assert to_rna("atcg") == "AUCG"
        assert to_rna("TTAAGGCC") == "UUAAGGCC"
    
    def test_translate_dna_cds(self):
        """Test DNA CDS translation."""
        # Test standard translation
        assert translate_dna_cds("ATG") == "M"
        assert translate_dna_cds("ATGTAA") == "M*"
        assert translate_dna_cds("ATGGCC") == "MA"
        
        # Test with stop codon (ensure length is multiple of 3)
        cds = "ATGAAAGAGAAGATGTTCAAAAGAAAACATTCACAAAATGGGTAAATGCACAATTTTCTAAG"  # 62 chars, pad to 63
        cds_padded = cds + "A"  # Make it 63 chars (multiple of 3)
        # Test that we can translate a proper CDS
        result = translate_dna_cds(cds_padded)
        assert len(result) == 21  # 63/3 = 21 amino acids
        assert result.startswith("M")  # Should start with methionine


class TestTranscriptFunctions:
    """Test transcript-related functions."""
    
    def test_pick_transcript_canonical(self):
        """Test picking canonical transcript."""
        transcripts = [
            {"id": "t1", "is_canonical": False},
            {"id": "t2", "is_canonical": True},
            {"id": "t3", "is_canonical": False}
        ]
        result = pick_transcript(transcripts, prefer_canonical=True)
        assert result["id"] == "t2"
    
    def test_pick_transcript_fallback(self):
        """Test fallback to longest CDS when no canonical."""
        transcripts = [
            {"id": "t1", "Translation": {"length": 100}},
            {"id": "t2", "Translation": {"length": 200}},
            {"id": "t3", "Translation": {"length": 150}}
        ]
        result = pick_transcript(transcripts, prefer_canonical=True)
        assert result["id"] == "t2"
    
    def test_transcript_exons_in_order_forward_strand(self):
        """Test exon ordering for forward strand."""
        transcript = {
            "strand": 1,
            "Exon": [
                {"start": 100, "end": 200},
                {"start": 300, "end": 400},
                {"start": 150, "end": 250}
            ]
        }
        result = transcript_exons_in_order(transcript)
        assert result[0]["start"] == 100
        assert result[1]["start"] == 150
        assert result[2]["start"] == 300
    
    def test_transcript_exons_in_order_reverse_strand(self):
        """Test exon ordering for reverse strand."""
        transcript = {
            "strand": -1,
            "Exon": [
                {"start": 100, "end": 200},
                {"start": 300, "end": 400},
                {"start": 150, "end": 250}
            ]
        }
        result = transcript_exons_in_order(transcript)
        assert result[0]["start"] == 300
        assert result[1]["start"] == 150
        assert result[2]["start"] == 100


class TestVariantParsing:
    """Test variant notation parsing."""
    
    def test_parse_genomic_variant_with_chrom(self):
        """Test parsing genomic variant with chromosome."""
        result = parse_variant_notation("chrX:g.31805775-31932165del")
        assert result.variant_type == "del"
        assert result.chrom == "X"
        assert result.start == 31805775
        assert result.end == 31932165
    
    def test_parse_genomic_variant_simple(self):
        """Test parsing simple genomic variant."""
        result = parse_variant_notation("g.123456_789012del")
        assert result.variant_type == "del"
        assert result.chrom == "unknown"
        assert result.start == 123456
        assert result.end == 789012
    
    def test_parse_duplication_variant(self):
        """Test parsing duplication variant."""
        result = parse_variant_notation("chrX:g.123456-789012dup")
        assert result.variant_type == "dup"
        assert result.chrom == "X"
        assert result.start == 123456
        assert result.end == 789012
    
    def test_parse_cDNA_variant(self):
        """Test parsing cDNA variant."""
        result = parse_variant_notation("c.123-456_789+012del")
        assert result.variant_type == "del"
        assert result.chrom == "cDNA"
        assert result.start == 123 - 456
        assert result.end == 789 + 12
    
    def test_parse_invalid_variant(self):
        """Test parsing invalid variant notation."""
        with pytest.raises(ValueError, match="Unable to parse variant notation"):
            parse_variant_notation("invalid_notation")


class TestVariantAffectedExons:
    """Test variant affected exon identification."""
    
    def test_get_variant_affected_exons_deletion(self):
        """Test identifying affected exons for deletion."""
        variant = VariantNotation("del", "X", 31805775, 31932165)
        exons = [
            {"start": 31932080, "end": 31932227, "strand": -1},  # Should be affected (overlaps end)
            {"start": 31929596, "end": 31929745, "strand": -1},  # Should be affected (within)
            {"start": 31875188, "end": 31875373, "strand": -1},  # Should be affected (within)
            {"start": 31836718, "end": 31836819, "strand": -1},  # Should be affected (within)
            {"start": 31819975, "end": 31820083, "strand": -1},  # Should be affected (overlaps start)
            {"start": 31773960, "end": 31774192, "strand": -1},  # Should NOT be affected (outside)
        ]
        
        affected = get_variant_affected_exons(variant, exons)
        assert len(affected) == 5
        assert affected[0]["start"] == 31932080
        assert affected[4]["start"] == 31819975
        
        # Check deletion types
        assert affected[0]["deletion_type"] == "partial"  # Exon 46: overlaps end
        assert affected[1]["deletion_type"] == "complete"  # Exon 47: within deletion
        assert affected[2]["deletion_type"] == "complete"  # Exon 48: within deletion
        assert affected[3]["deletion_type"] == "complete"  # Exon 49: within deletion
        assert affected[4]["deletion_type"] == "complete"  # Exon 50: within deletion
    
    def test_get_variant_affected_exons_duplication(self):
        """Test identifying affected exons for duplication."""
        variant = VariantNotation("dup", "X", 100, 200)
        exons = [
            {"start": 50, "end": 80, "strand": 1},   # Outside
            {"start": 90, "end": 110, "strand": 1},  # Overlaps
            {"start": 150, "end": 180, "strand": 1}, # Within
            {"start": 190, "end": 210, "strand": 1}, # Overlaps
            {"start": 250, "end": 280, "strand": 1}, # Outside
        ]
        
        affected = get_variant_affected_exons(variant, exons)
        assert len(affected) == 3
        assert affected[0]["start"] == 90
        assert affected[1]["start"] == 150
        assert affected[2]["start"] == 190
    
    def test_get_variant_affected_exons_partial_deletion(self):
        """Test identifying affected exons for partial deletion."""
        variant = VariantNotation("del", "X", 31820000, 31830000)
        exons = [
            {"start": 31819975, "end": 31820083, "strand": -1},  # Should be affected (partial)
            {"start": 31836718, "end": 31836819, "strand": -1},  # Should NOT be affected (outside)
        ]
        
        affected = get_variant_affected_exons(variant, exons)
        assert len(affected) == 1
        assert affected[0]["start"] == 31819975
        assert affected[0]["deletion_type"] == "partial"
    
    def test_get_variant_affected_exons_intronic_deletion(self):
        """Test identifying affected exons for intronic deletion."""
        variant = VariantNotation("del", "X", 31840000, 31850000)
        exons = [
            {"start": 31819975, "end": 31820083, "strand": -1},  # Should NOT be affected (outside)
            {"start": 31836718, "end": 31836819, "strand": -1},  # Should NOT be affected (outside)
        ]
        
        affected = get_variant_affected_exons(variant, exons)
        assert len(affected) == 0


class TestVariantProteinGeneration:
    """Test variant protein sequence generation."""
    
    def test_generate_variant_protein_sequence_deletion(self):
        """Test generating protein sequence for deletion variant."""
        variant = VariantNotation("del", "X", 31805775, 31932165)
        
        # Create affected exons (exons 46-50)
        affected_exons = [
            {"start": 31932080, "end": 31932227, "strand": -1},
            {"start": 31929596, "end": 31929745, "strand": -1},
            {"start": 31875188, "end": 31875373, "strand": -1},
            {"start": 31836718, "end": 31836819, "strand": -1},
            {"start": 31819975, "end": 31820083, "strand": -1},
        ]
        
        # Create all exons (including unaffected ones)
        all_exons = affected_exons + [
            {"start": 31773960, "end": 31774192, "strand": -1},  # Exon 51 - not affected
        ]
        
        # Create protein segments
        protein_segments = [
            {"exon_number": 46, "protein_seq": "ELQDGIGQRQTVVRTLNATGEEIIQQSSKTDASILQEKLGSLNLRWQEVCKQLSDRKK"},
            {"exon_number": 47, "protein_seq": "RLEEQKNILSEFQRDLNEFVLWLEEADNIASIPLEPGKEQQLKEKLEQV"},
            {"exon_number": 48, "protein_seq": "LLVEELPLRQGILKQLNETGGPVLVSAPISPEEQDKLENKLKQTNLQWIK"},
            {"exon_number": 49, "protein_seq": "VSRALPEKQGEIEAQIKDLGQLEKKLEDLEEQLNHLLLWLSPIRNQLEIYNQPNQEGPFDVK"},
            {"exon_number": 50, "protein_seq": "ETEIAVQAKQPDVEEILSKGQHLYKEKPATQPVK"},
            {"exon_number": 51, "protein_seq": "RKLEDLSSEWKAVNRLLQELRAKQPDLAPGLTTIGA"},
        ]
        
        result = generate_variant_protein_sequence(variant, protein_segments, affected_exons, all_exons)
        
        # Should only include exon 51's protein sequence (the unaffected exon)
        expected = "RKLEDLSSEWKAVNRLLQELRAKQPDLAPGLTTIGA"
        assert result == expected
    
    def test_generate_variant_protein_sequence_duplication(self):
        """Test generating protein sequence for duplication variant."""
        variant = VariantNotation("dup", "X", 100, 200)
        affected_exons = []
        all_exons = []
        protein_segments = [
            {"exon_number": 1, "protein_seq": "ABC"},
            {"exon_number": 2, "protein_seq": "DEF"},
        ]
        
        result = generate_variant_protein_sequence(variant, protein_segments, affected_exons, all_exons)
        # For duplications, should return original sequence (placeholder implementation)
        expected = "ABCDEF"
        assert result == expected
    
    def test_generate_variant_protein_sequence_partial_deletion(self):
        """Test generating protein sequence for partial deletion variant."""
        variant = VariantNotation("del", "X", 31820000, 31830000)
        
        # Create affected exons (exon 50 with partial deletion)
        affected_exons = [
            {"start": 31819975, "end": 31820083, "strand": -1, "deletion_type": "partial"}
        ]
        
        # Create all exons
        all_exons = [
            {"start": 31819975, "end": 31820083, "strand": -1},  # Exon 50 - partially affected
            {"start": 31836718, "end": 31836819, "strand": -1},  # Exon 49 - not affected
        ]
        
        # Create protein segments
        protein_segments = [
            {"exon_number": 50, "protein_seq": "RKLEDLSSEWKAVNRLLQELRAKQPDLAPGLTTIGA"},
            {"exon_number": 49, "protein_seq": "ETEIAVQAKQPDVEEILSKGQHLYKEKPATQPVK"},
        ]
        
        result = generate_variant_protein_sequence(variant, protein_segments, affected_exons, all_exons)
        
        # Should include exon 49 (unaffected) and potentially part of exon 50
        # The exact result depends on the partial deletion calculation logic
        assert len(result) > 0
        assert "ETEIAVQAKQPDVEEILSKGQHLYKEKPATQPVK" in result  # Exon 49 should be included


class TestVariantReport:
    """Test variant report generation."""
    
    def test_create_variant_report(self):
        """Test creating variant report."""
        variant = VariantNotation("del", "X", 31805775, 31932165)
        affected_exons = [
            {"start": 31932080, "end": 31932227, "strand": -1},
            {"start": 31929596, "end": 31929745, "strand": -1},
        ]
        variant_seq = "CACTGCTGATGATAAGAAATATTTTTTCCCCTGGAGGAAGCCTTAAGAAGATTAGGGAAA"
        protein_segments = [
            {"exon_number": 46, "protein_seq": "ELQDGIGQRQTVVRTLNATGEEIIQQSSKTDASILQEKLGSLNLRWQEVCKQLSDRKK"},
            {"exon_number": 47, "protein_seq": "RLEEQKNILSEFQRDLNEFVLWLEEADNIASIPLEPGKEQQLKEKLEQV"},
        ]
        all_exons = [
            {"start": 31932080, "end": 31932227, "strand": -1},
            {"start": 31929596, "end": 31929745, "strand": -1},
        ]
        
        report = create_variant_report(variant, affected_exons, variant_seq, protein_segments, all_exons)
        
        assert "Variant: DEL(X:31805775-31932165)" in report
        assert "Number of affected exons: 2" in report
        assert "Exon 1: 31932080-31932227" in report
        assert "Exon 2: 31929596-31929745" in report
        assert "CACTGCTGATGATAAGAAATATTTTTTCCCCTGGAGGAAGCCTTAAGAAGATTAGGGAAA" in report
        assert "ELQDGIGQRQTVVRTLNATGEEIIQQSSKTDASILQEKLGSLNLRWQEVCKQLSDRKK" in report
