# exonSeq

A small CLI to fetch DNA, RNA (cDNA), and translated protein sequences for a gene/transcript from Ensembl, and export per-exon mappings.

## Installation

```bash
python -m venv .venv
. .venv/Scripts/activate  # on Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Usage

exonSeq provides two main commands: `export` for gene sequence analysis and `exon` for exon utilities.

### Basic gene analysis
```bash
python -m exonseq export \
  --species homo_sapiens \
  --gene DMD \
  --reference GRCh38 \
  --coords X:31100000-33300000:1 \
  --out out
```

### With variant analysis (GRCh38/hg38)
```bash
python -m exonseq export \
  --species homo_sapiens \
  --gene DMD \
  --reference GRCh38 \
  --variant "chrX:g.31100000-31200000del" \
  --out out
```

### With variant analysis (GRCh37/hg19)
```bash
python -m exonseq export \
  --species homo_sapiens \
  --gene DMD \
  --reference GRCh37 \
  --variant "chrX:g.31137345-31200000del" \
  --out out_hg19
```

### Exon utilities
```bash
# Convert exon report to JSON
python -m exonseq exon convert exon_report.txt exon_report.json

# Generate FASTA from exon range
python -m exonseq exon fasta-range exon_report.txt 1 9 output.fasta --type protein

# Generate FASTA from specific exons
python -m exonseq exon fasta-list exon_report.txt 1 3 5 output.fasta --type protein

# List all exons in report
python -m exonseq exon list exon_report.txt

# Get detailed exon information
python -m exonseq exon info exon_report.txt 10

# Generate formatted protein FASTA
python -m exonseq exon formatted-protein exon_report.txt formatted_protein.fasta
```

### Version information
```bash
python -m exonseq --version
```

### Parameters

#### Export command parameters:
- `--species`: Ensembl species name, e.g., `homo_sapiens`, `mus_musculus`.
- `--gene`: HGNC/MGI symbol. Used to pick a transcript (canonical preferred).
- `--reference`: Reference genome/coord system version for sequence endpoints, e.g., `GRCh38` (hg38), `GRCh37` (hg19), `GRCm39`.
- `--coords`: Genomic region `CHR:START-END:STRAND` for region DNA export. Strand must be `1` or `-1`.
- `--transcript`: Optional Ensembl transcript ID to override automatic selection.
- `--variant`: Variant notation for deletion/duplication analysis (see supported formats below).
- `--out`: Output directory (created if missing).
- `--no-report`: Skip generating reproducibility report.

#### Exon command parameters:
- `convert`: Convert exon report to JSON format
- `fasta-range`: Generate FASTA from exon range
- `fasta-list`: Generate FASTA from specific exon numbers
- `list`: List all exons in report
- `info`: Get detailed information for specific exon
- `formatted-protein`: Generate formatted protein FASTA with spaces between exons

### Supported Variant Notations
- `chrX:g.123456-789012del` - Genomic deletion with chromosome
- `chrX:g.123456-789012dup` - Genomic duplication with chromosome  
- `g.123456_789012del` - Genomic deletion without chromosome
- `g.123456_789012dup` - Genomic duplication without chromosome
- `c.123-456_789+012del` - cDNA deletion (experimental)
- `c.123-456_789+012dup` - cDNA duplication (experimental)
- `123456-789012del` - Simple genomic range deletion
- `123456-789012dup` - Simple genomic range duplication

### Outputs
- `dna_region.fasta` (region DNA)
- `transcript_cdna.fasta` (spliced cDNA)
- `transcript_cds.fasta` (CDS only)
- `protein.fasta` (translated protein)
- `exon_report.txt` (per-exon DNA/RNA/protein coordinates and residues)
- `variant_report.txt` (variant analysis report, if --variant provided)
- `variant_sequence.fasta` (variant sequence, if --variant provided)

## Notes
- Uses Ensembl REST API (`https://rest.ensembl.org`).
- Transcript selection prefers canonical; falls back to longest CDS.
