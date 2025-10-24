"""Integration tests for exonSeq CLI functionality."""

import pytest
import os
import subprocess
import tempfile
from pathlib import Path


class TestCLIIntegration:
    """Integration tests for CLI functionality."""
    
    def test_basic_gene_export(self, temp_dir):
        """Test basic gene export without variants."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should succeed
        assert result.returncode == 0
        
        # Check that required files were created
        expected_files = [
            "transcript_cdna.fasta",
            "transcript_cds.fasta", 
            "protein.fasta",
            "exon_report.txt"
        ]
        
        for filename in expected_files:
            filepath = os.path.join(temp_dir, filename)
            assert os.path.exists(filepath), f"Expected file {filename} not found"
            assert os.path.getsize(filepath) > 0, f"File {filename} is empty"
    
    def test_variant_analysis_deletion(self, temp_dir):
        """Test variant analysis with deletion."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--variant", "chrX:g.31805775-31932165del",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should succeed
        assert result.returncode == 0
        
        # Check that variant-specific files were created
        expected_files = [
            "variant_report.txt",
            "variant_sequence.fasta",
            "variant_protein.fasta"
        ]
        
        for filename in expected_files:
            filepath = os.path.join(temp_dir, filename)
            assert os.path.exists(filepath), f"Expected file {filename} not found"
            assert os.path.getsize(filepath) > 0, f"File {filename} is empty"
        
        # Check variant report content
        with open(os.path.join(temp_dir, "variant_report.txt"), 'r') as f:
            report_content = f.read()
            
        assert "Variant: DEL(X:31805775-31932165)" in report_content
        assert "Number of affected exons: 5" in report_content
        assert "Exon 46:" in report_content
        assert "Exon 47:" in report_content
        assert "Exon 48:" in report_content
        assert "Exon 49:" in report_content
        assert "Exon 50:" in report_content
        
        # Check variant sequence file
        with open(os.path.join(temp_dir, "variant_sequence.fasta"), 'r') as f:
            variant_seq_content = f.read()
            
        assert variant_seq_content.startswith(">del_X_31805775_31932165")
        # Check that we have sequence content (may be formatted in 60-char lines)
        sequence_lines = variant_seq_content.split('\n')[1:]
        total_sequence_length = sum(len(line) for line in sequence_lines if line.strip())
        assert total_sequence_length > 100  # Should have substantial sequence
        
        # Check variant protein file
        with open(os.path.join(temp_dir, "variant_protein.fasta"), 'r') as f:
            variant_protein_content = f.read()
            
        assert variant_protein_content.startswith(">del_X_31805775_31932165_protein")
        protein_seq = ''.join(variant_protein_content.split('\n')[1:])
        assert len(protein_seq) > 100  # Should have substantial protein sequence
    
    def test_variant_analysis_duplication(self, temp_dir):
        """Test variant analysis with duplication."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--variant", "chrX:g.123456-789012dup",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should succeed
        assert result.returncode == 0
        
        # Check that variant-specific files were created
        expected_files = [
            "variant_report.txt",
            "variant_sequence.fasta",
            "variant_protein.fasta"
        ]
        
        for filename in expected_files:
            filepath = os.path.join(temp_dir, filename)
            assert os.path.exists(filepath), f"Expected file {filename} not found"
    
    def test_invalid_gene(self, temp_dir):
        """Test handling of invalid gene symbol."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "INVALID_GENE",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should fail
        assert result.returncode != 0
        assert "No valid lookup found" in result.stderr or "error" in result.stderr.lower()
    
    def test_invalid_variant_notation(self, temp_dir):
        """Test handling of invalid variant notation."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--variant", "invalid_variant_notation",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should fail during variant parsing
        assert result.returncode != 0
        assert "Unable to parse variant notation" in result.stdout
    
    def test_output_directory_creation(self, temp_dir):
        """Test that output directory is created if it doesn't exist."""
        new_temp_dir = os.path.join(temp_dir, "new_subdir")
        
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--out", new_temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should succeed and create the directory
        assert result.returncode == 0
        assert os.path.exists(new_temp_dir)
        assert os.path.exists(os.path.join(new_temp_dir, "protein.fasta"))
    
    def test_file_contents_validity(self, temp_dir):
        """Test that generated files contain valid content."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "homo_sapiens",
            "--gene", "DMD",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        
        # Check FASTA file format
        with open(os.path.join(temp_dir, "protein.fasta"), 'r') as f:
            protein_content = f.read()
            
        assert protein_content.startswith(">")
        lines = protein_content.strip().split('\n')
        assert len(lines) >= 2
        assert all(line.startswith('>') or all(c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*' for c in line) for line in lines)
        
        # Check exon report format
        with open(os.path.join(temp_dir, "exon_report.txt"), 'r') as f:
            exon_content = f.read()
            
        lines = exon_content.strip().split('\n')
        assert len(lines) > 1
        header = lines[0]
        assert "exon_number" in header
        assert "genomic_start" in header
        assert "genomic_end" in header
        assert "protein_seq" in header


class TestErrorHandling:
    """Test error handling scenarios."""
    
    def test_missing_required_arguments(self):
        """Test that missing required arguments cause appropriate errors."""
        # Missing --gene
        cmd = ["python", "-m", "exonseq", "export", "--species", "homo_sapiens"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode != 0
        
        # Missing --species
        cmd = ["python", "-m", "exonseq", "export", "--gene", "DMD"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode != 0
        
        # Missing --out
        cmd = ["python", "-m", "exonseq", "export", "--species", "homo_sapiens", "--gene", "DMD"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode != 0
    
    def test_invalid_species(self, temp_dir):
        """Test handling of invalid species."""
        cmd = [
            "python", "-m", "exonseq", "export",
            "--species", "invalid_species",
            "--gene", "DMD",
            "--out", temp_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode != 0


class TestFilePermissions:
    """Test file permission handling."""
    
    def test_readonly_output_directory(self, temp_dir):
        """Test handling of read-only output directory."""
        # Create a read-only subdirectory
        readonly_dir = os.path.join(temp_dir, "readonly")
        os.makedirs(readonly_dir)
        
        import platform
        if platform.system() == "Windows":
            pytest.skip("Read-only directory test not reliable on Windows")
            
        try:
            os.chmod(readonly_dir, 0o444)  # Read-only
            
            cmd = [
                "python", "-m", "exonseq", "export",
                "--species", "homo_sapiens",
                "--gene", "DMD",
                "--out", readonly_dir
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # Should fail due to permission error
            assert result.returncode != 0
        except (OSError, PermissionError):
            # Skip test on systems where chmod doesn't work as expected
            pytest.skip("Cannot create read-only directory on this system")
        
        # Restore permissions for cleanup
        os.chmod(readonly_dir, 0o755)
