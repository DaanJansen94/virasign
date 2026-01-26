# RefSelector

A pipeline for selecting optimal viral reference sequences from databases (RVDB, RefSeq, or custom accession). Maps reads using minimap2, calculates coverage metrics, filters by identity thresholds, and generates comprehensive interactive HTML reports with filterable tables, charts, and heatmaps.

## Installation

### 1. Create Conda Environment with Required Tools

Create a new conda environment with all required tools:

```bash
conda create -n refselector python=3.9 -y
conda activate refselector
conda install -c bioconda -c conda-forge minimap2 seqtk curl samtools -y
```

**Note**: `seqtk` and `curl` are optional (Python fallbacks are available), but `minimap2` is required.

### 2. Install RefSelector

```bash
git clone https://github.com/DaanJansen94/refselector.git
cd refselector
pip install .
```

## Usage

### Basic Usage

```bash
refselector -i <input_dir> -o <output_dir> [options]
```

**Note**: The database argument (`-d/--database`) is optional and defaults to `RVDB` if not specified.

### Command-Line Options

#### Required Arguments
- `-i, --input`: Input directory containing FASTQ files
- `-o, --output`: Output directory for results (use `.` to create `Refselector_output/` in current directory)

#### Optional Arguments
- `-d, --database`: Database name, accession number, text file with accessions, or path to FASTA file (default: `RVDB`)
  - **Database names**: `RVDB`, `RefSeq`, or `RVDB,RefSeq` (automatically downloads to `Databases/` directory)
  - **Accession number**: Single NCBI accession (e.g., `OZ254622.1`) - downloads and uses as database
  - **Text file with accessions**: Path to text file with one accession per line - downloads all and merges into database
  - **File path**: Path to your own FASTA database file
  - **Default**: If not specified, defaults to `RVDB`
  - Examples:
    ```bash
    # Use default RVDB database (no -d needed)
    refselector -i input/ -o output/
    
    # Use RefSeq database
    refselector -i input/ -d RefSeq -o output/
    
    # Use both RVDB and RefSeq (separate outputs for each)
    refselector -i input/ -d RVDB,RefSeq -o output/
    
    # Use a single accession as database
    refselector -i input/ -d OZ254622.1 -o output/
    
    # Use text file with accessions as database
    refselector -i input/ -d accessions.txt -o output/
    # (accessions.txt contains one accession per line, e.g.:
    #  OZ254622.1
    #  NC_123456.1
    #  AY123456.1)
    
    # Use accession as database and merge with additional accessions
    refselector -i input/ -d OZ254622.1 -a NC_123456.1,AY123456.1 -o output/
    
    # Use text file with accessions and merge with additional accessions
    refselector -i input/ -d accessions.txt -a PX852146.1 -o output/
    
    # Use custom database file
    refselector -i input/ -d /path/to/database.fasta -o output/
    ```

- `-a, --accession`: NCBI accession number(s) to download and merge with database
  - **Comma-separated list**: `-a PX852146.1,NC_123456.1,AY123456.1`
  - **Text file** (one accession per line): `-a accessions.txt`
  - Examples:
    ```bash
    # Comma-separated list
    refselector -i input/ -d RVDB -o output/ -a PX852146.1,NC_123456.1
    
    # Text file
    refselector -i input/ -d RVDB -o output/ -a accessions.txt
    ```

- `--min_identity`: Minimum average identity percentage
  - Default: `95%` for RefSeq, `80%` for RVDB
  - Example: `--min_identity 90`
- `--min_mapped_reads`: Minimum number of mapped reads (default: `100`)
- `--coverage_depth_threshold`: Minimum coverage depth threshold (default: `1.0`)
- `--coverage_breadth_threshold`: Minimum coverage breadth threshold (default: `0.1`)
- `-t, --threads`: Number of threads for minimap2 (default: `1`)

### Examples

```bash
# Basic usage with default RVDB database
refselector -i samples/ -o results/

# Use RefSeq with custom identity threshold
refselector -i samples/ -d RefSeq -o results/ --min_identity 97

# Use both databases with custom accessions
refselector -i samples/ -d RVDB,RefSeq -o results/ -a PX852146.1,NC_123456.1 -t 16

# Use a single accession as the database
refselector -i samples/ -d OZ254622.1 -o results/

# Use text file with accessions as database
refselector -i samples/ -d my_accessions.txt -o results/
# (my_accessions.txt contains one accession per line)

# Use accession as database and merge with additional accessions
refselector -i samples/ -d OZ254622.1 -a NC_123456.1,AY123456.1 -o results/

# Use text file with accessions and merge with additional accessions
refselector -i samples/ -d my_accessions.txt -a PX852146.1 -o results/

# Use custom database with accessions from file
refselector -i samples/ -d /path/to/database.fasta -o results/ -a my_accessions.txt
```

## Output Structure

When using a single database:
```
output/
├── refselector.log
├── results_summary_RefSeq.html  (or results_summary_RVDB.html)
└── sample_name/
    ├── sample_name_final_selected_references.json
    ├── sample_name_unfiltered_all_references.json
    └── PP826286.1/  (per-reference folder)
        ├── PP826286.1.fasta
        ├── PP826286.1.bam
        └── PP826286.1_mapped_reads.fastq
```

When using multiple databases (e.g., `RVDB,RefSeq`):
```
output/
├── refselector.log
├── results_summary_RefSeq.html
├── results_summary_RVDB.html
└── sample_name/
    ├── RVDB/
    │   ├── sample_name_final_selected_references.json
    │   ├── sample_name_unfiltered_all_references.json
    │   └── PP826286.1/  (per-reference folder)
    │       ├── PP826286.1.fasta
    │       ├── PP826286.1.bam
    │       └── PP826286.1_mapped_reads.fastq
    └── RefSeq/
        ├── sample_name_final_selected_references.json
        ├── sample_name_unfiltered_all_references.json
        └── NC_123456.1/  (per-reference folder)
            ├── NC_123456.1.fasta
            ├── NC_123456.1.bam
            └── NC_123456.1_mapped_reads.fastq
```

## Features

- **Automatic Database Download**: Downloads RVDB and RefSeq databases automatically to `Databases/` directory
- **NCBI Accession Support**: Download and merge custom accessions from NCBI
- **Database-Specific Thresholds**: Automatically uses appropriate identity thresholds (95% for RefSeq, 80% for RVDB)
- **Segmented Virus Support**: Handles segmented viruses (e.g., Lassa, Influenza) by keeping best reference for each segment
- **Maps reads** to reference database using minimap2
- **Calculates coverage depth and breadth** for each reference
- **Filters and curates references** based on identity and coverage thresholds
- **Re-maps reads** to selected references for accurate read counts
- **Automatically indexes databases** for faster mapping
- **Downloads and uses NCBI taxonomy database** for organism identification
- **Deduplicates references** by organism (preferring RefSeq over GenBank)
- **Interactive HTML Reports**: Generates comprehensive, interactive HTML visualization reports with:
  - Filterable and sortable tables per sample
  - Interactive charts (bar plots) for mapped reads, identity, coverage depth, and breadth
  - Heatmap showing viral species across all samples
  - Download functionality (CSV for tables, PNG for charts and heatmap)
  - Pipeline configuration display
  - Separate HTML files for each database (e.g., `results_summary_RefSeq.html`, `results_summary_RVDB.html`)

## Output Files

- **`results_summary_*.html`**: Interactive HTML visualization reports (one per database)
  - Filterable and sortable tables with viral reference statistics
  - Interactive bar charts for mapped reads, identity, coverage depth, and breadth
  - Heatmap visualization showing viral species presence across samples
  - Download options: CSV export for tables, PNG export for charts and heatmap
  - Sample selector to view results per sample
- **`*_final_selected_references.json`**: Final curated references after re-mapping (with accurate stats)
- **`*_unfiltered_all_references.json`**: All references found during initial mapping (with original stats)
- **Per-reference folders**: Each curated reference gets its own folder with:
  - `{accession}.fasta`: Reference sequence
  - `{accession}.bam`: BAM alignment file (indexed)
  - `{accession}_mapped_reads.fastq`: Reads that mapped to that specific reference

## Notes

- Databases are downloaded to `Databases/` directory in the current working directory
- Downloaded databases are cached and reused on subsequent runs
- Accessions are temporarily downloaded, merged with the database, then cleaned up
- The tool automatically handles database indexing for faster subsequent runs
- For segmented viruses, the tool keeps the best reference for each segment (e.g., L and S segments for Lassa virus)
- All SAM files are automatically removed to save space (BAM files are kept)
