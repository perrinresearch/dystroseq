# exonSeq Test Suite

This directory contains comprehensive unit tests and integration tests for the exonSeq bioinformatics tool.

## Test Structure

```
tests/
├── __init__.py                 # Test package initialization
├── conftest.py                 # Pytest configuration and fixtures
├── test_seq_utils.py          # Unit tests for sequence utilities
├── test_ensembl_client.py     # Unit tests for Ensembl API client
├── test_integration.py        # Integration tests for CLI functionality
└── README.md                  # This file
```

## Test Categories

### Unit Tests

**test_seq_utils.py** - Tests for core sequence manipulation functions:
- Basic sequence operations (DNA→RNA conversion, translation)
- Transcript processing (exon ordering, protein segment mapping)
- Variant parsing and notation handling
- Exon identification for variants
- Variant protein sequence generation
- Report generation

**test_ensembl_client.py** - Tests for Ensembl API interactions:
- Gene lookup by symbol
- Sequence retrieval (DNA, cDNA, CDS, protein)
- Error handling and retry logic
- Rate limiting handling
- Response parsing

### Integration Tests

**test_integration.py** - Tests for complete CLI workflows:
- Basic gene export functionality
- Variant analysis with deletions
- Variant analysis with duplications
- Error handling for invalid inputs
- File output validation
- Directory creation and permissions

## Running Tests

### Prerequisites

Install test dependencies:
```bash
pip install -r requirements.txt
```

### Quick Test Commands

```bash
# Run all tests
make test

# Run only unit tests
make test-unit

# Run only integration tests
make test-integration

# Run tests with coverage
make test-coverage

# Run fast tests (no network)
make test-fast

# Run example functionality test
make test-example
```

### Manual Test Commands

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_seq_utils.py

# Run specific test class
pytest tests/test_seq_utils.py::TestVariantParsing

# Run specific test function
pytest tests/test_seq_utils.py::TestVariantParsing::test_parse_genomic_variant_with_chrom

# Run with verbose output
pytest -v

# Run with coverage
pytest --cov=dystroseq --cov-report=html

# Run tests matching a pattern
pytest -k "variant"

# Run tests excluding network tests
pytest -m "not network"
```

## Test Data and Fixtures

### Fixtures (conftest.py)

- `temp_dir`: Creates temporary directories for test outputs
- `sample_exons`: Sample exon data with genomic coordinates
- `sample_protein_segments`: Sample protein segment data
- `sample_variant_notations`: Various variant notation formats for testing

### Mock Data

Tests use mock data to avoid dependency on external services:
- Mock Ensembl API responses
- Sample genomic coordinates and sequences
- Test variant notations

## Key Test Scenarios

### Exon Identification Tests

Tests verify that the exon identification logic correctly:
1. Identifies exons affected by deletion variants
2. Handles overlapping exon boundaries
3. Properly numbers exons in reports
4. Handles different variant types (del, dup, ins)

### Variant Protein Generation Tests

Tests verify that variant protein sequences:
1. Remove affected exons for deletion variants
2. Maintain correct sequence ordering
3. Generate valid FASTA format output
4. Handle edge cases (no affected exons, all exons affected)

### Integration Tests

Tests verify complete workflows:
1. CLI argument parsing and validation
2. File generation and output validation
3. Error handling for invalid inputs
4. Network error handling
5. File permission handling

## Continuous Integration

The test suite is configured to run automatically on:
- Push to main/develop branches
- Pull requests
- Multiple Python versions (3.9, 3.10, 3.11, 3.12)

### CI Pipeline

1. **Unit Tests**: Fast, isolated tests that don't require network access
2. **Integration Tests**: Tests that may require network access (marked to skip if unavailable)
3. **Linting**: Code style and type checking
4. **Coverage**: Code coverage reporting
5. **Example Test**: End-to-end functionality verification

## Test Markers

Tests are marked for different categories:

- `@pytest.mark.unit`: Unit tests (fast, isolated)
- `@pytest.mark.integration`: Integration tests (may require network)
- `@pytest.mark.slow`: Slow-running tests
- `@pytest.mark.network`: Tests requiring network access

## Coverage Goals

- **Unit Tests**: >95% coverage for core functions
- **Integration Tests**: Cover main user workflows
- **Error Handling**: Test all error conditions

## Adding New Tests

### For New Features

1. Add unit tests for individual functions in `test_seq_utils.py`
2. Add integration tests for CLI functionality in `test_integration.py`
3. Add fixtures for test data in `conftest.py`
4. Update this README with new test scenarios

### Test Naming Conventions

- Test files: `test_*.py`
- Test classes: `Test*`
- Test functions: `test_*`
- Use descriptive names that explain what is being tested

### Example Test Structure

```python
def test_function_name_scenario(self):
    """Test description explaining what is being tested."""
    # Arrange - Set up test data
    input_data = "test_input"
    expected_output = "expected_result"
    
    # Act - Execute the function
    result = function_under_test(input_data)
    
    # Assert - Verify the result
    assert result == expected_output
```

## Troubleshooting

### Common Issues

1. **Network Tests Failing**: Integration tests may fail if Ensembl API is unavailable
   - Use `make test-fast` to skip network tests
   - Check internet connectivity

2. **Permission Errors**: Some tests create temporary files
   - Ensure write permissions in test directory
   - Run with appropriate user permissions

3. **Missing Dependencies**: Ensure all test dependencies are installed
   ```bash
   pip install -r requirements.txt
   ```

### Debug Mode

Run tests with verbose output and no capture:
```bash
pytest -v -s
```

### Test Isolation

Each test should be independent:
- Use fixtures for setup/teardown
- Clean up temporary files
- Mock external dependencies
- Reset global state between tests

