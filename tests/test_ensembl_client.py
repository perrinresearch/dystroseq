"""Unit tests for dystroseq.ensembl_client module."""

import pytest
from unittest.mock import Mock, patch, MagicMock
from dystroseq.ensembl_client import EnsemblClient


class TestEnsemblClient:
    """Test EnsemblClient functionality."""
    
    @pytest.fixture
    def client(self):
        """Create a test client instance."""
        return EnsemblClient(base_url="https://test.ensembl.org", timeout_s=5, max_retries=2)
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_lookup_gene_by_symbol_success(self, mock_request, client):
        """Test successful gene lookup by symbol."""
        # Mock successful response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "id": "ENSG00000198947",
            "display_name": "DMD",
            "seq_region_name": "X"
        }
        mock_request.return_value = mock_response
        
        result = client.lookup_gene_by_symbol("homo_sapiens", "DMD")
        
        assert result["id"] == "ENSG00000198947"
        assert result["display_name"] == "DMD"
        assert result["seq_region_name"] == "X"
        
        # Verify the request was made correctly
        mock_request.assert_called_once()
        call_args = mock_request.call_args
        assert "homo_sapiens/DMD" in call_args[0][1]  # URL is the second positional argument
        assert call_args[1]["params"]["expand"] == 1
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_lookup_gene_by_symbol_not_found(self, mock_request, client):
        """Test gene lookup when gene is not found."""
        # Mock 404 response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_response.text = '{"error":"No valid lookup found for symbol INVALID"}'
        mock_request.return_value = mock_response
        
        with pytest.raises(RuntimeError, match="Ensembl API error 404"):
            client.lookup_gene_by_symbol("homo_sapiens", "INVALID")
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_sequence_region_success(self, mock_request, client):
        """Test successful sequence region retrieval."""
        # Mock successful response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "ATCGATCGATCG"
        mock_request.return_value = mock_response
        
        result = client.sequence_region("homo_sapiens", "X:100-200:1")
        
        assert result == "ATCGATCGATCG"
        mock_request.assert_called_once()
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_sequence_region_json_response(self, mock_request, client):
        """Test sequence region with JSON response."""
        # Mock JSON response (some APIs return JSON instead of plain text)
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '{"seq": "ATCGATCGATCG"}'
        mock_response.startswith.return_value = True
        mock_request.return_value = mock_response
        
        result = client.sequence_region("homo_sapiens", "X:100-200:1")
        
        assert result == "ATCGATCGATCG"
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_sequence_id_success(self, mock_request, client):
        """Test successful sequence ID retrieval."""
        # Mock successful response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "ATCGATCGATCG"
        mock_request.return_value = mock_response
        
        result = client.sequence_id("ENST00000357033", "cdna")
        
        assert result == "ATCGATCGATCG"
        
        # Verify the request was made correctly
        mock_request.assert_called_once()
        call_args = mock_request.call_args
        assert "ENST00000357033" in call_args[0][1]  # URL is the second positional argument
        assert call_args[1]["params"]["type"] == "cdna"
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_retry_on_rate_limit(self, mock_request, client):
        """Test retry mechanism on rate limiting."""
        # Mock rate limit response followed by success
        rate_limit_response = Mock()
        rate_limit_response.status_code = 429
        rate_limit_response.headers = {"Retry-After": "1"}
        
        success_response = Mock()
        success_response.status_code = 200
        success_response.json.return_value = {"id": "test"}
        
        mock_request.side_effect = [rate_limit_response, success_response]
        
        result = client.lookup_id("test_id")
        
        assert result["id"] == "test"
        assert mock_request.call_count == 2
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_max_retries_exceeded(self, mock_request, client):
        """Test behavior when max retries are exceeded."""
        # Mock rate limit response
        rate_limit_response = Mock()
        rate_limit_response.status_code = 429
        rate_limit_response.headers = {"Retry-After": "1"}
        
        mock_request.return_value = rate_limit_response
        
        with pytest.raises(RuntimeError, match="Exhausted retries"):
            client.lookup_id("test_id")
        
        # Should have tried max_retries times (client has max_retries=2)
        assert mock_request.call_count == client.max_retries
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_gene_transcripts(self, mock_request, client):
        """Test gene transcripts retrieval."""
        # Mock successful response with transcripts
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "Transcript": [
                {"id": "ENST00000357033", "is_canonical": True},
                {"id": "ENST00000470390", "is_canonical": False}
            ]
        }
        mock_request.return_value = mock_response
        
        result = client.gene_transcripts("ENSG00000198947")
        
        assert len(result) == 2
        assert result[0]["id"] == "ENST00000357033"
        assert result[1]["id"] == "ENST00000470390"
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_lookup_id_with_expand(self, mock_request, client):
        """Test lookup ID with expand parameter."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"id": "test", "display_name": "test"}
        mock_request.return_value = mock_response
        
        result = client.lookup_id("test_id", expand=True)
        
        assert result["id"] == "test"
        mock_request.assert_called_once()
        call_args = mock_request.call_args
        assert call_args[1]["params"]["expand"] == 1
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_lookup_id_without_expand(self, mock_request, client):
        """Test lookup ID without expand parameter."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"id": "test", "display_name": "test"}
        mock_request.return_value = mock_response
        
        result = client.lookup_id("test_id", expand=False)
        
        assert result["id"] == "test"
        mock_request.assert_called_once()
        call_args = mock_request.call_args
        assert call_args[1]["params"] is None
    
    def test_client_initialization(self):
        """Test client initialization with different parameters."""
        client = EnsemblClient(
            base_url="https://custom.ensembl.org",
            timeout_s=60,
            max_retries=10
        )
        
        assert client.base_url == "https://custom.ensembl.org"
        assert client.timeout_s == 60
        assert client.max_retries == 10
        
        # Test default values
        default_client = EnsemblClient()
        assert default_client.base_url == "https://rest.ensembl.org"
        assert default_client.timeout_s == 30
        assert default_client.max_retries == 5
    
    @patch('dystroseq.ensembl_client.requests.Session.request')
    def test_sequence_region_with_coord_system_version(self, mock_request, client):
        """Test sequence region with coordinate system version."""
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = "ATCGATCGATCG"
        mock_request.return_value = mock_response
        
        result = client.sequence_region("homo_sapiens", "X:100-200:1", "GRCh38")
        
        assert result == "ATCGATCGATCG"
        mock_request.assert_called_once()
        call_args = mock_request.call_args
        assert call_args[1]["params"]["coord_system_version"] == "GRCh38"
