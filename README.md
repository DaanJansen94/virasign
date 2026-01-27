# Virasign

Virasign (Viral Read **ASSIGN**ment from nanopore sequencing) is a viral classification and reference selection tool for nanopore data. It maps long-read sequencing data (via minimap2) against viral databases (RVDB, RefSeq, or custom) and generates comprehensive interactive HTML reports with filterable tables, charts and heatmaps. For each identified virus, Virasign also provides the closest reference sequence, mapped reads in FASTQ format, and BAM files which can be used to easily generate a consensus genome and visualize data.

## Installation

### 1. Create Conda Environment with Required Tools

Create a new conda environment with all required tools:

```bash
conda create -n virasign python=3.9 -y
conda activate virasign
conda install -c bioconda -c conda-forge minimap2=2.24 seqtk=1.3 curl=7.88 samtools=1.17 -y
```

### 2. Install Virasign

```bash
git clone https://github.com/DaanJansen94/virasign.git
cd virasign
pip install .
```

## Usage

### Basic Usage

```bash
virasign -i input_dir -o output_dir [options]
```

**Note**: The database argument (`-d/--database`) is optional and defaults to `RVDB` if not specified.

### Command-Line Options

#### Required Arguments
- `-i, --input`: Input directory containing FASTQ files
- `-o, --output`: Output directory for results (use `.` to create `Virasign_output/` in current directory)

#### Optional Arguments
- `-d, --database`: Database name (`RefSeq`, `RVDB` (default), or accession number)
- `-a, --accession`: NCBI accession number(s) to download and merge with database
- `--min_identity`: Minimum average identity percentage
- `--min_mapped_reads`: Minimum number of mapped reads (default: `100`)
- `--coverage_depth_threshold`: Minimum coverage depth threshold (default: `1.0`)
- `--coverage_breadth_threshold`: Minimum coverage breadth threshold (default: `0.1`)
- `-t, --threads`: Number of threads for minimap2 (default: `1`)

### Examples

```bash
# Basic usage with default RVDB database
virasign -i input_dir -o output_dir

# Use both databases with custom accessions and 16 threads
virasign -i input_dir -d RVDB,RefSeq -o output_dir -a PX852146.1,NC_123456.1 -t 16

# Use a single accession as the database
virasign -i input_dir -d OZ254622.1 -o output_dir

# Use text file with accessions as database
virasign -i input_dir -d my_accessions.txt -o output_dir
# (my_accessions.txt contains one accession per line)
```

## Output Structure

When using a single database:
```
output/
├── virasign.log
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
├── virasign.log
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

## Output Files

- **results_summary_*.html**: Interactive HTML reports (one per database)
- ***_final_selected_references.json**: Final curated references
- **Per-reference folders**: Each curated reference gets its own folder with:
  - **{accession}.fasta**: Reference sequence
  - **{accession}.bam**: BAM alignment file (indexed)
  - **{accession}_mapped_reads.fastq**: Reads that mapped to that specific reference
