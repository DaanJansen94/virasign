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

## Output Files

- **results_summary_*.html**: Interactive HTML reports (one per database)
- ***_final_selected_references.json**: Final curated references
- **Per-reference folders**: Each curated reference gets its own folder with:
  - **{accession}.fasta**: Reference sequence
  - **{accession}.bam**: BAM alignment file (indexed with {accession}.bam.bai)
  - **{accession}_mapped_reads.fastq**: Reads that mapped to that specific reference

## HTML Output example

Example HTML output for one sample:

![HTML Output Example](html_example.png)

## Citation

If you use Lassaseq in your research, please cite:

```
Jansen, D., & Vercauteren, K. (2026). Virasign: A viral classification tool designed for nanopore sequencing data (v0.0.1). Zenodo. https://doi.org/10.5281/zenodo.18387008
```

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.
