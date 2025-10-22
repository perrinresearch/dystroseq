"""Test configuration and fixtures for dystroSeq tests."""

import pytest
import tempfile
import os
from typing import Dict, List, Any


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_exons():
    """Sample exon data for testing."""
    return [
        {
            "start": 31932080,
            "end": 31932227,
            "strand": -1,
            "id": "exon_46"
        },
        {
            "start": 31929596,
            "end": 31929745,
            "strand": -1,
            "id": "exon_47"
        },
        {
            "start": 31875188,
            "end": 31875373,
            "strand": -1,
            "id": "exon_48"
        },
        {
            "start": 31836718,
            "end": 31836819,
            "strand": -1,
            "id": "exon_49"
        },
        {
            "start": 31819975,
            "end": 31820083,
            "strand": -1,
            "id": "exon_50"
        },
        {
            "start": 31773960,
            "end": 31774192,
            "strand": -1,
            "id": "exon_51"
        }
    ]


@pytest.fixture
def sample_protein_segments():
    """Sample protein segments for testing."""
    return [
        {
            "exon_number": 46,
            "protein_start_index": 2147,
            "protein_end_index": 2204,
            "protein_seq": "ELQDGIGQRQTVVRTLNATGEEIIQQSSKTDASILQEKLGSLNLRWQEVCKQLSDRKK"
        },
        {
            "exon_number": 47,
            "protein_start_index": 2205,
            "protein_end_index": 2253,
            "protein_seq": "RLEEQKNILSEFQRDLNEFVLWLEEADNIASIPLEPGKEQQLKEKLEQV"
        },
        {
            "exon_number": 48,
            "protein_start_index": 2255,
            "protein_end_index": 2304,
            "protein_seq": "LLVEELPLRQGILKQLNETGGPVLVSAPISPEEQDKLENKLKQTNLQWIK"
        },
        {
            "exon_number": 49,
            "protein_start_index": 2305,
            "protein_end_index": 2366,
            "protein_seq": "VSRALPEKQGEIEAQIKDLGQLEKKLEDLEEQLNHLLLWLSPIRNQLEIYNQPNQEGPFDVK"
        },
        {
            "exon_number": 50,
            "protein_start_index": 2367,
            "protein_end_index": 2400,
            "protein_seq": "ETEIAVQAKQPDVEEILSKGQHLYKEKPATQPVK"
        },
        {
            "exon_number": 51,
            "protein_start_index": 2401,
            "protein_end_index": 2436,
            "protein_seq": "RKLEDLSSEWKAVNRLLQELRAKQPDLAPGLTTIGA"
        }
    ]


@pytest.fixture
def sample_variant_notations():
    """Sample variant notations for testing."""
    return {
        "genomic_del": "chrX:g.31805775-31932165del",
        "genomic_dup": "chrX:g.123456-789012dup",
        "simple_del": "g.123456_789012del",
        "cDNA_del": "c.123-456_789+012del",
        "invalid": "invalid_notation"
    }

