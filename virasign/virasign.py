#!/usr/bin/env python3

import os
import subprocess
import sys
from pathlib import Path
import argparse
import logging
import shutil
import json
import multiprocessing
import re
import hashlib
import time
import urllib.request
import urllib.error
import urllib.parse
import xml.etree.ElementTree as ET
import gzip
import zipfile
import sqlite3

def setup_logging(output_dir):
    """Set up logging configuration."""
    log_file = Path(output_dir) / 'virasign.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),  # 'w' mode overwrites the file
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

logger = logging.getLogger(__name__)

def extract_accession_from_header(header: str) -> str:
    """
    Extract accession number from FASTA header.
    Accession is the part between the 2nd and 3rd pipe (|).
    Example: 'acc|REFSEQ|NC_001474.2|Dengue virus 2, complete genome' -> 'NC_001474.2'
    If header doesn't have the expected format, returns the first whitespace-delimited token.
    """
    if not header:
        return ""
    parts = header.split("|")
    if len(parts) >= 3:
        # Accession is between 2nd and 3rd pipe
        return parts[2].strip()
    # Fallback: return first token if no pipes
    return header.split()[0] if header.split() else ""

def normalize_species_name(species: str) -> str:
    """
    Normalize species name for matching.
    Removes 'type' and normalizes numbers to allow matching between formats like:
    'dengue virus type 2' and 'dengue virus 2'
    """
    if not species:
        return ""
    # Lowercase and strip
    normalized = species.lower().strip()
    # Remove 'type' word (e.g., "dengue virus type 2" -> "dengue virus 2")
    normalized = normalized.replace(" type ", " ")
    normalized = normalized.replace("type ", "")
    # Remove extra spaces
    normalized = " ".join(normalized.split())
    return normalized

def extract_species_from_header(header: str) -> str:
    """
    Extract species from FASTA header.
    For GENBANK format: 'acc|GENBANK|...|dengue virus type 2|VRL|03-SEP-2010'
    Species is the part before the second-to-last pipe (third from last).
    For REFSEQ format: 'acc|REFSEQ|NC_001474.2|Dengue virus 2, complete genome'
    Species is extracted from description field (between 3rd and 4th pipe), before comma.
    Returns normalized species name (lowercase, trimmed, 'type' removed).
    """
    if not header:
        return ""
    parts = header.split("|")
    
    # Check if it's GENBANK format with VRL|date at end
    if len(parts) >= 4 and parts[-1].strip() and parts[-2].strip() == "VRL":
        # GENBANK format: ...|species|VRL|date
        species = parts[-3].strip()
    # Check if it's REFSEQ format (has REFSEQ in second position)
    elif len(parts) >= 4 and parts[1].strip() == "REFSEQ":
        # REFSEQ format: acc|REFSEQ|accession|description
        description = parts[3].strip()
        # Extract species from description (before comma, or first 2-3 words)
        if "," in description:
            species = description.split(",")[0].strip()
        else:
            # Take first 3 words as species name
            words = description.split()
            species = " ".join(words[:3]) if len(words) >= 3 else description
    # Fallback: try second-to-last part
    elif len(parts) >= 2:
        species = parts[-2].strip()
    else:
        species = header.split()[0] if header.split() else ""
    
    return normalize_species_name(species)

def fetch_segment_from_ncbi_genbank(accession: str, max_retries: int = 3) -> str:
    """
    Fetch segment information from NCBI GenBank record using Entrez API.
    Parses the GenBank record to extract /segment="X" from FEATURES section.
    Returns segment identifier (e.g., 'S', 'L', '1', '2', etc.) or empty string if not found.
    """
    import time
    import urllib.request
    import urllib.error
    
    # Remove version from accession if present
    accession_base = accession.split('.')[0] if '.' in accession else accession
    
    for attempt in range(max_retries):
        try:
            time.sleep(0.1)  # Rate limiting
            
            # Use efetch to get GenBank record
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_base}&rettype=gb&retmode=text"
            with urllib.request.urlopen(fetch_url, timeout=10) as response:
                gb_record = response.read().decode('utf-8')
                
                # Parse GenBank record to find /segment="X" in FEATURES section
                import re
                # Look for /segment="..." pattern in FEATURES section
                segment_match = re.search(r'/segment="([^"]+)"', gb_record)
                if segment_match:
                    segment = segment_match.group(1).strip()
                    return segment
                
                return ""
        except urllib.error.HTTPError as e:
            if e.code == 404:
                # Accession not found
                return ""
            if attempt < max_retries - 1:
                time.sleep(1)
            else:
                logger.debug(f"Failed to fetch GenBank record for {accession}: {e}")
        except Exception as e:
            if attempt < max_retries - 1:
                logger.debug(f"Attempt {attempt + 1} failed for {accession}: {e}")
                time.sleep(1)
            else:
                logger.debug(f"Failed to fetch segment for {accession} after {max_retries} attempts: {e}")
    
    return ""

def get_segment_from_database(accession: str, database_path: Path) -> str:
    """
    Get segment information from SQLite database.
    Returns segment identifier or empty string if not found.
    """
    if not database_path or not database_path.exists():
        return ""
    
    try:
        import sqlite3
        conn = sqlite3.connect(str(database_path))
        cursor = conn.cursor()
        
        # Try versioned accession first
        cursor.execute('SELECT segment FROM accession_to_segment WHERE accession = ?', (accession,))
        result = cursor.fetchone()
        if result:
            conn.close()
            return result[0] if result[0] else ""
        
        # Try base accession (without version)
        accession_base = accession.split('.')[0] if '.' in accession else accession
        cursor.execute('SELECT segment FROM accession_to_segment WHERE accession = ?', (accession_base,))
        result = cursor.fetchone()
        if result:
            conn.close()
            return result[0] if result[0] else ""
        
        conn.close()
        return ""
    except Exception as e:
        logger.debug(f"Error reading segment database: {e}")
        return ""

def save_segment_to_database(accession: str, segment: str, database_path: Path):
    """
    Save segment information to SQLite database.
    """
    try:
        import sqlite3
        conn = sqlite3.connect(str(database_path))
        cursor = conn.cursor()
        
        # Create table if it doesn't exist
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS accession_to_segment (
                accession TEXT PRIMARY KEY,
                segment TEXT
            )
        ''')
        
        # Create index for faster lookups
        cursor.execute('''
            CREATE INDEX IF NOT EXISTS idx_accession_segment ON accession_to_segment(accession)
        ''')
        
        # Insert or replace
        cursor.execute('INSERT OR REPLACE INTO accession_to_segment (accession, segment) VALUES (?, ?)', 
                      (accession, segment))
        
        # Also save base accession (without version)
        accession_base = accession.split('.')[0] if '.' in accession else accession
        if accession_base != accession:
            cursor.execute('INSERT OR REPLACE INTO accession_to_segment (accession, segment) VALUES (?, ?)', 
                          (accession_base, segment))
        
        conn.commit()
        conn.close()
    except Exception as e:
        logger.debug(f"Error saving to segment database: {e}")

def fetch_segment_from_ncbi(accession: str, database_path: Path = None, max_retries: int = 3, silent: bool = False) -> str:
    """
    Fetch segment information for an accession.
    First checks database, then fetches from NCBI GenBank if not found.
    
    Args:
        accession: Accession number
        database_path: Path to segment database
        max_retries: Maximum retry attempts
        silent: If True, don't log debug messages (for batch processing)
    """
    # Step 1: Check database first (already done in deduplicate_by_organism, but keep for direct calls)
    if database_path:
        segment = get_segment_from_database(accession, database_path)
        if segment:
            return segment
    
    # Step 2: Fetch from NCBI GenBank record
    segment = fetch_segment_from_ncbi_genbank(accession, max_retries)
    
    # Step 3: Save to database if found
    if segment and database_path:
        save_segment_to_database(accession, segment, database_path)
    
    return segment

def extract_segment_from_description(description: str, accession: str = None, database_path: Path = None) -> str:
    """
    Extract segment identifier from description for segmented viruses.
    
    First tries to get segment from NCBI GenBank database (if accession provided),
    then falls back to parsing description text.
    
    Handles various formats found in RefSeq and GenBank databases:
    - Arenaviruses (Lassa, etc.): "segment L", "segment S", "L segment", "S segment"
    - Influenza: "segment 1", "segment 2", etc. (or PB2, PB1, PA, HA, NP, NA, M, NS)
    - Bunyaviruses: "segment L", "segment M", "segment S"
    - Reoviruses: "segment 1", "segment 2", etc.
    - Multi-segmented: "segment RNA 1", "segment RNA2", "segment S1", "segment L2", "segment M1"
    
    Examples from databases:
    - "Lassa virus strain Z148 segment S, complete sequence"
    - "Influenza A virus segment 1, complete sequence"
    - "High Plains wheat mosaic virus segment RNA 2, complete sequence"
    - "Mahlapitsi virus segment S1, complete sequence"
    - "segment L, complete genome"
    
    Args:
        description: FASTA header description
        accession: Accession number (optional, used to query NCBI database)
        database_path: Path to segment database (optional)
    
    Returns segment identifier (e.g., 'L', 'S', '1', '2', 'M', 'PB2', 'S1', 'L2', etc.) 
    or empty string if not segmented.
    """
    # Step 1: Try to get segment from NCBI GenBank database (most reliable)
    if accession:
        segment = fetch_segment_from_ncbi(accession, database_path)
        if segment:
            return segment
    
    # Step 2: Fall back to parsing description text (less reliable but works offline)
    if not description:
        return ""
    
    description_lower = description.lower()
    
    import re
    
    # Valid segment identifiers
    valid_influenza = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NS']
    
    # Pattern 1a: "segment RNA X" or "segment RNAX" (e.g., "segment RNA 2", "segment RNA2")
    # Extract just the number, ignore "RNA"
    match = re.search(r'segment\s+rna\s+(\d+)(?:[,\.\s]|$)', description_lower)
    if match:
        return match.group(1)
    match = re.search(r'segment\s+rna(\d+)(?:[,\.\s]|$)', description_lower)
    if match:
        return match.group(1)
    
    # Pattern 1b: "segment X" where X is letter+number (e.g., "segment S1", "segment L2", "segment M1")
    # These are sub-segments and should be treated as distinct segments
    match = re.search(r'segment\s+([lsm]\d+)(?:[,\.\s]|$)', description_lower)
    if match:
        return match.group(1).upper()
    
    # Pattern 1c: "segment X" where X is single letter (L, S, M) - allow word boundary, comma, period, or whitespace after
    # Note: description_lower is lowercase, so match 'l', 's', 'm'
    match = re.search(r'segment\s+([lsm])(?:\b|[,\.\s]|$)', description_lower)
    if match:
        segment = match.group(1).upper()
        return segment
    
    # Pattern 1d: "segment X" where X is a number (segment 1, segment 2, etc.)
    match = re.search(r'segment\s+(\d+)(?:[,\.\s]|$)', description_lower)
    if match:
        return match.group(1)
    
    # Pattern 1e: "segment X" where X is influenza segment name (PB2, PB1, PA, HA, NP, NA, NS)
    for seg in valid_influenza:
        match = re.search(r'segment\s+' + re.escape(seg.lower()) + r'(?:\b|[,\.\s]|$)', description_lower)
        if match:
            return seg
    
    # Pattern 2: "X segment" (where X is L, S, or M) - less common but exists
    match = re.search(r'\b([lsm])\s+segment', description_lower)
    if match:
        return match.group(1).upper()
    
    # Pattern 3: "segment X RNA" or "segment X genome" - handles variations
    match = re.search(r'segment\s+([lsm])\s+(?:rna|genome|sequence)', description_lower)
    if match:
        return match.group(1).upper()
    
    match = re.search(r'segment\s+(\d+)\s+(?:rna|genome|sequence)', description_lower)
    if match:
        return match.group(1)
    
    # Skip "segment DNA" - this is molecule type, not a segment identifier
    # (we already checked for numbers and letters above, so this won't match)
    
    return ""

def is_refseq(description: str) -> bool:
    """Check if a reference is from RefSeq based on description."""
    return "|REFSEQ|" in description or description.startswith("acc|REFSEQ|")

def get_taxonomy_cache_path(output_dir: Path) -> Path:
    """Get path to individual accession-to-organism JSON cache."""
    return output_dir / "taxonomy_cache" / "accession_to_organism.json"

def get_taxonomy_database_path(output_dir: Path) -> Path:
    """Get path to comprehensive taxonomy database (SQLite format, faster than JSON)."""
    return output_dir / "taxonomy_cache" / "accession_to_organism_database.db"

def get_segment_database_path(output_dir: Path) -> Path:
    """Get path to segment database (SQLite format, stores accession -> segment mappings)."""
    cache_dir = output_dir / "taxonomy_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / "accession_to_segment_database.db"

def get_taxonomy_database_json_path(output_dir: Path) -> Path:
    """Get path to comprehensive taxonomy database JSON (legacy/fallback)."""
    return output_dir / "taxonomy_cache" / "accession_to_organism_database.json"

def get_ncbi_files_dir(database_dir: Path) -> Path:
    """Get directory for downloaded NCBI taxonomy dump files."""
    return database_dir / "taxonomy_cache" / "ncbi_files"

def get_virasign_databases_dir() -> Path:
    """Get path to virasign Databases directory (where downloaded databases are stored)."""
    # Create Databases directory in current working directory (where virasign is run from)
    # This way databases are stored in the user's workspace, not in conda environment
    databases_dir = Path.cwd() / "Databases"
    databases_dir.mkdir(parents=True, exist_ok=True)
    return databases_dir

def download_rvdb_database(databases_dir: Path, force_download: bool = False, rvdb_version: str = "C-RVDBv31.0.fasta.gz", enable_clustering: bool = False, cluster_identity: float = 0.98) -> Path:
    """
    Download and prepare RVDB database with complete genomes only, optionally cluster at specified identity.
    
    Args:
        databases_dir: Directory where databases are stored
        force_download: Force re-download even if database exists
        rvdb_version: RVDB version filename (e.g., "C-RVDBv30.0.fasta.gz" or "C-RVDBvCurrent.fasta.gz")
        enable_clustering: Whether to cluster sequences (default: True)
        cluster_identity: Identity threshold for clustering, e.g., 0.98 for 98% (default: 0.98)
    
    Returns path to RVDB database (e.g., RVDBv30.0_complete.fasta)
    """
    rvdb_dir = databases_dir / "RVDB"
    rvdb_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract version number from filename (e.g., "C-RVDBv30.0.fasta.gz" -> "v30.0")
    version_match = rvdb_version.replace("C-RVDB", "").replace(".fasta.gz", "")
    version_suffix = version_match if version_match else "Current"
    
    output_fasta = rvdb_dir / f"RVDB{version_suffix}_complete.fasta"
    
    # Check if already exists
    if output_fasta.exists() and not force_download:
        logger.info(f"RVDB database already exists: {output_fasta}")
        
        # Clean up any intermediate files that might exist
        compressed_file = rvdb_dir / rvdb_version
        temp_fasta = rvdb_dir / rvdb_version.replace(".gz", "").replace(".fasta", ".fasta")
        cluster_prefix = rvdb_dir / f"RVDB{version_suffix}_clu95"
        if compressed_file.exists():
            compressed_file.unlink()
            logger.info(f"  Cleaned up: {compressed_file.name}")
        if temp_fasta.exists():
            temp_fasta.unlink()
            logger.info(f"  Cleaned up: {temp_fasta.name}")
        # Clean up clustering temp files
        if (rvdb_dir / "tmp_clu95").exists():
            shutil.rmtree(rvdb_dir / "tmp_clu95")
        for f in rvdb_dir.glob(f"RVDB{version_suffix}_clu95*"):
            if f != output_fasta:
                f.unlink()
        
        # Check if metadata exists, if not create it
        metadata_file = rvdb_dir / "database_metadata.json"
        if not metadata_file.exists():
            logger.info("Creating metadata file for existing database...")
            metadata = {
                "database_name": "RVDB",
                "version": version_suffix,
                "download_date": "Unknown (database existed before metadata tracking)",
                "source_url": f"https://rvdb.dbi.udel.edu/download/{rvdb_version}",
                "filtered": True,
                "filter_criteria": "complete genomes only" + (", clustered" if enable_clustering else ""),
                "output_file": str(output_fasta.name),
                "note": "Metadata created retroactively"
            }
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
        return output_fasta
    
    compressed_file = rvdb_dir / rvdb_version
    temp_fasta = rvdb_dir / rvdb_version.replace(".gz", "")
    
    # Download if needed
    if not compressed_file.exists() or force_download:
        logger.info(f"Downloading RVDB database version {rvdb_version} (this may take a while, ~few GB)...")
        url = f"https://rvdb.dbi.udel.edu/download/{rvdb_version}"
        try:
            def show_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) // total_size)
                    if block_num % 1000 == 0:
                        logger.info(f"  Download progress: {percent}%")
            
            urllib.request.urlretrieve(url, compressed_file, show_progress)
            logger.info(f"Successfully downloaded {compressed_file.name}")
        except Exception as e:
            logger.error(f"Failed to download RVDB: {e}")
            if compressed_file.exists():
                compressed_file.unlink()
            raise
    
    # Extract if needed
    if not temp_fasta.exists() or force_download:
        logger.info("Extracting RVDB database...")
        try:
            with gzip.open(compressed_file, 'rb') as f_in:
                with open(temp_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info("Extraction complete")
        except Exception as e:
            logger.error(f"Failed to extract RVDB: {e}")
            raise
    
    # Filter for complete genomes
    logger.info("Filtering for complete genomes...")
    complete_ids_file = rvdb_dir / "complete_ids.txt"
    
    try:
        # Extract IDs of complete genomes
        with open(temp_fasta, 'r') as f_in, open(complete_ids_file, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # Check if it's a complete genome
                    if any(term in line.lower() for term in [
                        'complete virus segment', 'complete virion genome',
                        'complete viral segment', 'complete viral genome',
                        'complete sequence', 'complete genome'
                    ]):
                        # Extract ID (everything after '>' up to first space or '|'
                        header = line[1:].strip()
                        # Get first 5 pipe-separated parts or first part
                        parts = header.split('|')[:5]
                        if len(parts) > 1:
                            f_out.write('|'.join(parts) + '\n')
                        else:
                            f_out.write(header.split()[0] + '\n')
        
        # Use seqtk to extract sequences (if available) or use Python
        logger.info("Extracting complete genome sequences...")
        try:
            # Try seqtk first (faster)
            result = subprocess.run(
                ['seqtk', 'subseq', str(temp_fasta), str(complete_ids_file)],
                stdout=open(output_fasta, 'w'),
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
            logger.info("Used seqtk for extraction")
        except (subprocess.CalledProcessError, FileNotFoundError):
            # Fallback to Python extraction (slower but works)
            logger.info("seqtk not found, using Python extraction (this may take longer)...")
            complete_ids = set()
            with open(complete_ids_file, 'r') as f:
                for line in f:
                    complete_ids.add(line.strip())
            
            with open(temp_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
                write_seq = False
                for line in f_in:
                    if line.startswith('>'):
                        header = line[1:].strip()
                        # Check if this header matches any complete ID
                        header_parts = header.split('|')[:5]
                        header_key = '|'.join(header_parts) if len(header_parts) > 1 else header.split()[0]
                        write_seq = header_key in complete_ids
                    
                    if write_seq:
                        f_out.write(line)
            
            logger.info("Python extraction complete")
        
        # Clean up temporary files
        if complete_ids_file.exists():
            complete_ids_file.unlink()
        
        complete_count = sum(1 for _ in open(output_fasta)) // 2  # Approximate (headers + sequences)
        logger.info(f"RVDB complete genomes database ready: {output_fasta} ({complete_count} sequences)")
        
        # Optional clustering step
        clustered_count = complete_count
        if enable_clustering:
            cluster_identity_pct = int(cluster_identity * 100)
            logger.info(f"Clustering complete genomes at {cluster_identity_pct}% identity (80% coverage of shortest sequence)...")
            mmseqs_available = shutil.which('mmseqs') is not None
            
            if mmseqs_available:
                cluster_prefix_name = f"RVDB{version_suffix}_clu{cluster_identity_pct}"
                cluster_prefix = rvdb_dir / cluster_prefix_name
                tmp_cluster_dir = rvdb_dir / f"tmp_clu{cluster_identity_pct}"
                
                try:
                    # Run mmseqs clustering
                    cluster_cmd = [
                        'mmseqs', 'easy-linclust',
                        str(output_fasta),
                        str(cluster_prefix),
                        str(tmp_cluster_dir),
                        '--min-seq-id', str(cluster_identity),
                        '-c', '0.8',
                        '--cov-mode', '0',
                        '--threads', '16'
                    ]
                
                    logger.info(f"Running: {' '.join(cluster_cmd)}")
                    result = subprocess.run(cluster_cmd, check=True, capture_output=True, text=True)
                    
                    # Rename clustered representative sequences to final output
                    clustered_fasta = rvdb_dir / f"{cluster_prefix_name}_rep_seq.fasta"
                    if clustered_fasta.exists():
                        # Replace the unclustered file with clustered version
                        output_fasta.unlink()
                        clustered_fasta.rename(output_fasta)
                        logger.info(f"Clustered database ready: {output_fasta}")
                        
                        # Count clustered sequences
                        clustered_count = sum(1 for _ in open(output_fasta)) // 2
                        logger.info(f"Clustering reduced {complete_count} sequences -> {clustered_count} representatives")
                    else:
                        logger.warning(f"Clustered output not found at {clustered_fasta}, keeping unclustered version")
                        clustered_count = complete_count
                    
                    # Clean up clustering intermediate files
                    logger.info("Cleaning up clustering intermediate files...")
                    if tmp_cluster_dir.exists():
                        shutil.rmtree(tmp_cluster_dir)
                        logger.info(f"  Removed: {tmp_cluster_dir.name}/")
                    # Remove other clustering output files (keep only the final FASTA)
                    for f in rvdb_dir.glob(f"{cluster_prefix_name}*"):
                        if f != output_fasta:
                            f.unlink()
                            logger.info(f"  Removed: {f.name}")
                    
                except (subprocess.CalledProcessError, FileNotFoundError) as e:
                    logger.warning(f"mmseqs clustering failed: {e}. Keeping unclustered database.")
                    clustered_count = complete_count
            else:
                logger.warning("mmseqs not found. Skipping clustering step. Install mmseqs2 for clustering.")
                clustered_count = complete_count
        else:
            logger.info("Clustering disabled. Using unclustered complete genomes database.")
        
        # Clean up intermediate files to save disk space
        logger.info("Cleaning up intermediate files...")
        if compressed_file.exists():
            compressed_file.unlink()
            logger.info(f"  Removed: {compressed_file.name}")
        if temp_fasta.exists():
            temp_fasta.unlink()
            logger.info(f"  Removed: {temp_fasta.name}")
        clustering_status = ""
        if enable_clustering:
            mmseqs_available = shutil.which('mmseqs') is not None
            if mmseqs_available:
                clustering_status = f" + clustered at {int(cluster_identity*100)}%"
        logger.info(f"  Kept: {output_fasta.name} (complete genomes{clustering_status})")
        
        # Save metadata file
        metadata_file = rvdb_dir / "database_metadata.json"
        filter_criteria = "complete genomes only"
        if enable_clustering:
            mmseqs_available = shutil.which('mmseqs') is not None
            if mmseqs_available:
                filter_criteria += f", clustered at {int(cluster_identity*100)}% identity (80% coverage)"
        
        metadata = {
            "database_name": "RVDB",
            "version": version_suffix,
            "download_date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "source_url": f"https://rvdb.dbi.udel.edu/download/{rvdb_version}",
            "filtered": True,
            "filter_criteria": filter_criteria,
            "output_file": str(output_fasta.name),
            "sequence_count": clustered_count,
            "clustered": enable_clustering and shutil.which('mmseqs') is not None,
            "cluster_identity": cluster_identity if enable_clustering else None,
            "note": "Intermediate files (compressed, uncurated, and clustering temp) were removed to save disk space"
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Saved database metadata to {metadata_file}")
        
    except Exception as e:
        logger.error(f"Failed to filter RVDB for complete genomes: {e}")
        raise
    
    return output_fasta

def download_refseq_database(databases_dir: Path, force_download: bool = False) -> Path:
    """
    Download and prepare RefSeq viral database.
    Returns path to viral_refseq_complete.fna
    """
    refseq_dir = databases_dir / "RefSeq"
    refseq_dir.mkdir(parents=True, exist_ok=True)
    
    output_fasta = refseq_dir / "viral_refseq_complete.fna"
    
    # Check if already exists
    if output_fasta.exists() and not force_download:
        logger.info(f"RefSeq database already exists: {output_fasta}")
        # Check if metadata exists, if not create it
        metadata_file = refseq_dir / "database_metadata.json"
        if not metadata_file.exists():
            logger.info("Creating metadata file for existing database...")
            from datetime import datetime
            metadata = {
                "database_name": "RefSeq Viral",
                "version": "v1.1",
                "download_date": "Unknown (database existed before metadata tracking)",
                "source_urls": [
                    "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz",
                    "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"
                ],
                "filtered": False,
                "output_file": str(output_fasta.name),
                "note": "Complete viral RefSeq genomes - metadata created retroactively"
            }
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2)
        return output_fasta
    
    logger.info("Downloading RefSeq viral database (this may take a while, ~few GB)...")
    
    # Use curl to download and pipe directly to gunzip (like: curl ... | gunzip > output)
    # This mimics: curl https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{1,2}.1.genomic.fna.gz | gunzip > viral_refseq_complete.fna
    try:
        # Check if curl and gunzip are available
        curl_available = shutil.which('curl') is not None
        gunzip_available = shutil.which('gunzip') is not None
        
        if curl_available and gunzip_available:
            # Use curl to download and pipe to gunzip (most efficient, like the user's command)
            logger.info("Using curl to download and extract directly (like: curl ... | gunzip > output)...")
            
            # Download part 1 (always exists)
            url1 = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
            logger.info(f"Downloading and extracting: {url1}")
            curl_proc1 = subprocess.Popen(
                ['curl', '--progress-bar', url1],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            gunzip_proc1 = subprocess.Popen(
                ['gunzip', '-c'],
                stdin=curl_proc1.stdout,
                stdout=open(output_fasta, 'wb'),
                stderr=subprocess.PIPE
            )
            curl_proc1.stdout.close()
            gunzip_proc1.communicate()
            exit_code1 = curl_proc1.wait()
            gunzip_exit1 = gunzip_proc1.wait()
            
            if exit_code1 != 0:
                raise subprocess.CalledProcessError(exit_code1, 'curl')
            
            # Try to download part 2 (may not exist)
            url2 = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"
            logger.info(f"Attempting to download: {url2}")
            curl_proc2 = subprocess.Popen(
                ['curl', '--progress-bar', '--fail', '--silent', url2],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            gunzip_proc2 = subprocess.Popen(
                ['gunzip', '-c'],
                stdin=curl_proc2.stdout,
                stdout=open(output_fasta, 'ab'),  # Append mode
                stderr=subprocess.PIPE
            )
            curl_proc2.stdout.close()
            gunzip_proc2.communicate()
            exit_code2 = curl_proc2.wait()
            
            if exit_code2 == 0:
                gunzip_proc2.wait()
                logger.info("Part 2 downloaded and appended")
            else:
                logger.info("Part 2 not available (only part 1 exists - this is normal)")
        
        else:
            # Fallback to Python urllib if curl/gunzip not available
            logger.info("Using Python to download RefSeq database (curl/gunzip not available)...")
            url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
            temp_file = refseq_dir / "viral.1.1.genomic.fna.gz"
            
            def show_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) // total_size)
                    if block_num % 1000 == 0:
                        logger.info(f"  Download progress: {percent}%")
            
            urllib.request.urlretrieve(url, temp_file, show_progress)
            logger.info("Extracting...")
            
            with gzip.open(temp_file, 'rb') as f_in:
                with open(output_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Try part 2 if it exists
            url2 = "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz"
            try:
                temp_file2 = refseq_dir / "viral.2.1.genomic.fna.gz"
                urllib.request.urlretrieve(url2, temp_file2, show_progress)
                logger.info("Extracting part 2...")
                with gzip.open(temp_file2, 'rb') as f_in:
                    with open(output_fasta, 'ab') as f_out:  # Append mode
                        shutil.copyfileobj(f_in, f_out)
                temp_file2.unlink()  # Clean up
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    logger.info("Part 2 not available (only part 1 exists - this is normal)")
                else:
                    raise
            
            temp_file.unlink()  # Clean up
        
        logger.info(f"RefSeq database ready: {output_fasta}")
        
        # Save metadata file
        metadata_file = refseq_dir / "database_metadata.json"
        # Try to get version info from file or use release date
        from datetime import datetime
        metadata = {
            "database_name": "RefSeq Viral",
            "version": "v1.1",  # RefSeq release version
            "download_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "source_urls": [
                "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
            ],
            "note": "Only part 1 exists (part 2 not available on server)",
            "filtered": False,
            "output_file": str(output_fasta.name),
            "note": "Complete viral RefSeq genomes"
        }
        
        # Try to count sequences
        try:
            seq_count = sum(1 for line in open(output_fasta) if line.startswith('>'))
            metadata["sequence_count"] = seq_count
        except:
            pass
        
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Saved database metadata to {metadata_file}")
        
    except Exception as e:
        logger.error(f"Failed to download RefSeq: {e}")
        # Clean up partial downloads
        for temp_file in temp_files:
            if temp_file.exists():
                temp_file.unlink()
        raise
    
    return output_fasta

def download_accession_from_ncbi(accession: str, output_dir: Path = None) -> Path:
    """
    Download a single accession from NCBI and save as FASTA file.
    Returns path to the downloaded FASTA file.
    """
    if output_dir is None:
        output_dir = get_virasign_databases_dir() / "Custom"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Remove version number for base accession (e.g., PX852146.1 -> PX852146)
    base_accession = accession.split('.')[0] if '.' in accession else accession
    
    output_fasta = output_dir / f"{base_accession}.fasta"
    
    # Check if already downloaded
    if output_fasta.exists():
        logger.info(f"Accession {accession} already downloaded: {output_fasta}")
        return output_fasta
    
    logger.info(f"Downloading accession {accession} from NCBI...")
    
    # NCBI Entrez API URL for FASTA download
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'nuccore',
        'id': base_accession,  # Use base accession (NCBI handles versioning)
        'rettype': 'fasta',
        'retmode': 'text'
    }
    
    try:
        # Build URL with parameters
        query_string = urllib.parse.urlencode(params)
        full_url = f"{url}?{query_string}"
        
        # Download with retry logic
        max_retries = 3
        for attempt in range(max_retries):
            try:
                with urllib.request.urlopen(full_url, timeout=30) as response:
                    fasta_data = response.read().decode('utf-8')
                    
                    if not fasta_data or fasta_data.strip().startswith('ERROR') or len(fasta_data) < 100:
                        raise ValueError(f"NCBI returned invalid data for accession {accession}")
                    
                    # Save to file
                    with open(output_fasta, 'w') as f:
                        f.write(fasta_data)
                    
                    logger.info(f"Successfully downloaded {accession} -> {output_fasta}")
                    return output_fasta
                    
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    logger.error(f"Accession {accession} not found in NCBI database (404)")
                    raise ValueError(f"Accession {accession} not found in NCBI")
                elif attempt < max_retries - 1:
                    logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}. Retrying...")
                    time.sleep(2)
                    continue
                else:
                    raise
            except Exception as e:
                if attempt < max_retries - 1:
                    logger.warning(f"Attempt {attempt + 1} failed for {accession}: {e}. Retrying...")
                    time.sleep(2)
                    continue
                else:
                    raise
                    
    except Exception as e:
        logger.error(f"Failed to download accession {accession}: {e}")
        raise

def merge_fasta_files(fasta_files: list, output_fasta: Path) -> Path:
    """
    Merge multiple FASTA files into a single file.
    Returns path to the merged FASTA file.
    """
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Merging {len(fasta_files)} FASTA file(s) into {output_fasta.name}...")
    
    with open(output_fasta, 'w') as outfile:
        for fasta_file in fasta_files:
            fasta_path = Path(fasta_file)
            if not fasta_path.exists():
                logger.warning(f"FASTA file not found: {fasta_path}. Skipping...")
                continue
            
            # Handle gzipped files
            if fasta_path.suffix == '.gz':
                with gzip.open(fasta_path, 'rt') as infile:
                    outfile.write(infile.read())
            else:
                with open(fasta_path, 'r') as infile:
                    outfile.write(infile.read())
    
    logger.info(f"Successfully merged FASTA files -> {output_fasta}")
    return output_fasta

def is_accession_number(text: str) -> bool:
    """
    Check if a string looks like an NCBI accession number.
    Accession formats: NC_123456.1, AY123456.1, PX852146.1, etc.
    Pattern: 1-2 letters, optional underscore, numbers, optional version number.
    """
    if not text or len(text) < 4:
        return False
    
    # Remove version number if present (e.g., OZ254622.1 -> OZ254622)
    base_text = text.split('.')[0] if '.' in text else text
    
    # Pattern: 1-2 letters, optional underscore, then digits
    import re
    pattern = r'^[A-Z]{1,2}_?\d+$'
    return bool(re.match(pattern, base_text.upper()))

def resolve_database_path(database_arg: str, accessions: list = None, enable_clustering: bool = False, cluster_identity: float = 0.98, rvdb_version: str = None) -> Path:
    """
    Resolve database argument to actual file path.
    Supports:
    - File paths: '/path/to/database.fasta' -> returns Path
    - Database names: 'RVDB', 'RefSeq', 'refseq' -> downloads and returns Path
    - Multiple databases: 'RVDB,RefSeq' -> combines and returns Path
    - Accession numbers: 'OZ254622.1' -> downloads and uses as database
    - Text file with accessions: '/path/to/accessions.txt' -> downloads all and merges
    
    If accessions are provided, they will be downloaded and merged with the database.
    """
    database_arg = database_arg.strip()
    databases_dir = get_virasign_databases_dir()
    
    # Check if database_arg is a single accession number
    if is_accession_number(database_arg):
        logger.info(f"Database argument is an accession number: {database_arg}")
        # Download the accession to Custom directory
        custom_dir = databases_dir / "Custom"
        custom_dir.mkdir(parents=True, exist_ok=True)
        
        # Download the main accession
        main_acc_fasta = download_accession_from_ncbi(database_arg, custom_dir)
        
        # If additional accessions provided via -a, merge them
        if accessions:
            logger.info(f"Downloading {len(accessions)} additional accession(s) to merge...")
            accession_fasta_files = [main_acc_fasta]
            for accession in accessions:
                try:
                    acc_fasta = download_accession_from_ncbi(accession.strip(), custom_dir)
                    accession_fasta_files.append(acc_fasta)
                except Exception as e:
                    logger.error(f"Failed to download accession {accession}: {e}")
                    raise
            
            # Merge all accessions
            base_accession = database_arg.split('.')[0] if '.' in database_arg else database_arg
            merged_fasta = custom_dir / f"{base_accession}_merged.fasta"
            result = merge_fasta_files(accession_fasta_files, merged_fasta)
            logger.info(f"Created merged database from {len(accession_fasta_files)} accession(s): {merged_fasta}")
            return result
        else:
            return main_acc_fasta
    
    # Check if database_arg is a text file with accessions (not a FASTA file)
    db_path = Path(database_arg)
    if db_path.exists() and db_path.is_file():
        # Check if it's a text file (not a FASTA file)
        # FASTA files typically have extensions: .fasta, .fa, .fna, .fas, .faa, .fq, .fastq
        fasta_extensions = {'.fasta', '.fa', '.fna', '.fas', '.faa', '.fq', '.fastq', '.gz'}
        if db_path.suffix.lower() not in fasta_extensions:
            # It's likely a text file with accessions
            logger.info(f"Reading accessions from file: {db_path}")
            accession_list = []
            with open(db_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):  # Skip empty lines and comments
                        accession_list.append(line)
            
            if not accession_list:
                raise ValueError(f"No accessions found in file: {db_path}")
            
            logger.info(f"Found {len(accession_list)} accession(s) in file")
            
            # Download all accessions to Custom directory
            custom_dir = databases_dir / "Custom"
            custom_dir.mkdir(parents=True, exist_ok=True)
            
            accession_fasta_files = []
            for accession in accession_list:
                try:
                    acc_fasta = download_accession_from_ncbi(accession.strip(), custom_dir)
                    accession_fasta_files.append(acc_fasta)
                except Exception as e:
                    logger.error(f"Failed to download accession {accession}: {e}")
                    raise
            
            # If additional accessions provided via -a, add them
            if accessions:
                logger.info(f"Downloading {len(accessions)} additional accession(s) from -a argument...")
                for accession in accessions:
                    try:
                        acc_fasta = download_accession_from_ncbi(accession.strip(), custom_dir)
                        accession_fasta_files.append(acc_fasta)
                    except Exception as e:
                        logger.error(f"Failed to download accession {accession}: {e}")
                        raise
            
            # Merge all accessions
            file_stem = db_path.stem
            merged_fasta = custom_dir / f"{file_stem}_database.fasta"
            result = merge_fasta_files(accession_fasta_files, merged_fasta)
            logger.info(f"Created merged database from {len(accession_fasta_files)} accession(s): {merged_fasta}")
            return result
    
    # Download accessions if provided (will determine target directory after we know which database)
    accession_fasta_files = []
    temp_dir = None
    if accessions:
        logger.info(f"Downloading {len(accessions)} accession(s) from NCBI...")
        # We'll download to a temp location first, then move to database directory after merging
        temp_dir = databases_dir / ".temp_accessions"
        temp_dir.mkdir(parents=True, exist_ok=True)
        for accession in accessions:
            try:
                acc_fasta = download_accession_from_ncbi(accession.strip(), temp_dir)
                accession_fasta_files.append(acc_fasta)
            except Exception as e:
                logger.error(f"Failed to download accession {accession}: {e}")
                raise
    
    # If it's a file path (contains '/' or exists as file), return as-is (but merge with accessions if provided)
    if '/' in database_arg or '\\' in database_arg or Path(database_arg).exists():
        base_database = Path(database_arg)
        
        # If accessions provided, merge with base database
        if accession_fasta_files:
            # Merge in the same directory as the base database
            db_dir = base_database.parent
            merged_fasta = db_dir / f"{base_database.stem}_with_accessions.fasta"
            all_fastas = [base_database] + accession_fasta_files
            result = merge_fasta_files(all_fastas, merged_fasta)
            
            # Clean up temporary accession files after merging
            if temp_dir and temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary accession download directory")
            
            return result
        else:
            return base_database
    
    # Parse database names (case-insensitive)
    db_names = [name.strip().lower() for name in database_arg.split(',')]
    
    fasta_files = []
    
    for db_name in db_names:
        if db_name == 'rvdb':
            # Convert version number to filename format (e.g., "30.0" -> "C-RVDBv30.0.fasta.gz")
            if rvdb_version:
                # Remove any leading/trailing whitespace and ensure it's in the right format
                version_clean = rvdb_version.strip()
                # If user provided just the number (e.g., "30.0"), convert to filename
                if not version_clean.startswith("C-RVDB"):
                    rvdb_filename = f"C-RVDBv{version_clean}.fasta.gz"
                else:
                    # User provided full filename
                    rvdb_filename = version_clean
                logger.info(f"Using RVDB version: {rvdb_filename}")
            else:
                # Default to v31.0
                rvdb_filename = "C-RVDBv31.0.fasta.gz"
                logger.info(f"Using default RVDB version: {rvdb_filename}")
            fasta_file = download_rvdb_database(databases_dir, rvdb_version=rvdb_filename, enable_clustering=enable_clustering, cluster_identity=cluster_identity)
            fasta_files.append(fasta_file)
        elif db_name == 'refseq':
            fasta_file = download_refseq_database(databases_dir)
            fasta_files.append(fasta_file)
        else:
            raise ValueError(f"Unknown database name: {db_name}. Supported: 'RVDB', 'RefSeq', or an accession number (e.g., 'OZ254622.1')")
    
    # If accessions provided, merge with each database directly in the database directory
    if accession_fasta_files:
        merged_files = []
        for db_fasta in fasta_files:
            # Merge accessions directly into the database directory (not Custom/)
            db_dir = db_fasta.parent
            db_name = db_fasta.stem
            # Create merged file in the same directory as the database
            merged_fasta = db_dir / f"{db_name}_with_accessions.fasta"
            all_fastas = [db_fasta] + accession_fasta_files
            merged_files.append(merge_fasta_files(all_fastas, merged_fasta))
        
        # Clean up temporary accession files after merging
        if temp_dir and temp_dir.exists():
            import shutil
            shutil.rmtree(temp_dir)
            logger.info(f"Cleaned up temporary accession download directory")
        
        fasta_files = merged_files
    
    # Return list of database files (no combining)
    # The caller will handle mapping against each separately
    if len(fasta_files) > 1:
        return fasta_files  # Return list when multiple databases
    else:
        return fasta_files[0]  # Return single file when one database

def load_taxonomy_cache(cache_path: Path) -> dict:
    """Load taxonomy cache from JSON file."""
    if cache_path.exists():
        try:
            with open(cache_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load taxonomy cache: {e}")
    return {}

def save_taxonomy_cache(cache_path: Path, cache: dict):
    """Save taxonomy cache to JSON file."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=2)

def download_ncbi_taxonomy_files(database_dir: Path, force_download: bool = False):
    """
    Download NCBI taxonomy dump files (one-time download).
    Downloads:
    - nucl_gb.accession2taxid.gz (GenBank nucleotide accession -> taxid mapping, includes RefSeq)
    - taxdmp.zip (contains names.dmp for taxid -> organism name mapping)
    """
    ncbi_dir = get_ncbi_files_dir(database_dir)
    ncbi_dir.mkdir(parents=True, exist_ok=True)
    
    # Use nucl_gb.accession2taxid.gz (GenBank nucleotide, includes RefSeq)
    # Note: nuccore.accession2taxid.gz doesn't exist - nucl_gb is the correct file
    accession2taxid_file = ncbi_dir / "nucl_gb.accession2taxid.gz"
    taxdmp_file = ncbi_dir / "taxdmp.zip"
    names_file = ncbi_dir / "names.dmp"
    
    # Download accession2taxid if needed
    if not accession2taxid_file.exists() or force_download:
        # Try primary URL (GenBank nucleotide, includes RefSeq)
        urls = [
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
        ]
        
        downloaded = False
        for url in urls:
            logger.info(f"Downloading {url} (this may take a while, ~few GB)...")
            logger.info(f"Saving to: {accession2taxid_file}")
            try:
                # Use urlretrieve with a callback to show progress
                def show_progress(block_num, block_size, total_size):
                    if total_size > 0:
                        percent = min(100, (block_num * block_size * 100) // total_size)
                        if block_num % 1000 == 0:  # Log every 1000 blocks to avoid spam
                            logger.info(f"  Download progress: {percent}%")
                
                urllib.request.urlretrieve(url, accession2taxid_file, show_progress)
                logger.info(f"Successfully downloaded {accession2taxid_file.name}")
                downloaded = True
                break
            except urllib.error.HTTPError as e:
                logger.warning(f"Failed to download from {url}: HTTP {e.code} - {e.reason}")
                if os.path.exists(accession2taxid_file):
                    os.remove(accession2taxid_file)  # Remove partial download
            except Exception as e:
                logger.warning(f"Failed to download from {url}: {e}")
                if os.path.exists(accession2taxid_file):
                    os.remove(accession2taxid_file)  # Remove partial download
        
        if not downloaded:
            logger.error(f"Failed to download accession2taxid from all attempted URLs")
            raise Exception("Could not download accession2taxid file from NCBI")
    
    # Download taxdmp.zip if needed
    if not taxdmp_file.exists() or force_download:
        url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
        logger.info(f"Downloading {url} (this may take a while, ~few hundred MB)...")
        logger.info(f"Saving to: {taxdmp_file}")
        try:
            # Use urlretrieve with a callback to show progress
            def show_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) // total_size)
                    if block_num % 100 == 0:  # Log every 100 blocks
                        logger.info(f"  Download progress: {percent}%")
            
            urllib.request.urlretrieve(url, taxdmp_file, show_progress)
            logger.info(f"Successfully downloaded {taxdmp_file.name}")
        except Exception as e:
            logger.error(f"Failed to download taxdmp: {e}")
            if os.path.exists(taxdmp_file):
                os.remove(taxdmp_file)  # Remove partial download
            raise
    
    # Extract names.dmp and nodes.dmp from taxdmp.zip if needed
    names_file = ncbi_dir / "names.dmp"
    nodes_file = ncbi_dir / "nodes.dmp"
    
    if not names_file.exists() or not nodes_file.exists() or force_download:
        logger.info(f"Extracting names.dmp and nodes.dmp from {taxdmp_file.name}...")
        try:
            with zipfile.ZipFile(taxdmp_file, 'r') as zip_ref:
                if not names_file.exists() or force_download:
                    zip_ref.extract("names.dmp", ncbi_dir)
                    logger.info(f"Extracted names.dmp")
                if not nodes_file.exists() or force_download:
                    zip_ref.extract("nodes.dmp", ncbi_dir)
                    logger.info(f"Extracted nodes.dmp")
        except Exception as e:
            logger.error(f"Failed to extract files from taxdmp.zip: {e}")
            raise

def get_virus_taxids(nodes_file: Path) -> set:
    """
    Find all taxids that are viruses (descendants of taxid 10239 "Viruses").
    Uses nodes.dmp to traverse the taxonomy tree.
    Returns set of virus taxids.
    """
    VIRUS_ROOT_TAXID = "10239"  # Viruses kingdom
    
    logger.info("Identifying all virus taxids from nodes.dmp...")
    
    # Build parent -> children mapping
    parent_to_children = {}
    taxid_to_parent = {}
    
    with open(nodes_file, 'r') as f:
        for line in f:
            # nodes.dmp format: taxid | parent taxid | rank | ... |
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.endswith('\t|'):
                line = line[:-2]
            parts = [p.strip() for p in line.split('\t|\t')]
            if len(parts) >= 2:
                taxid = parts[0]
                parent_taxid = parts[1]
                taxid_to_parent[taxid] = parent_taxid
                if parent_taxid not in parent_to_children:
                    parent_to_children[parent_taxid] = []
                parent_to_children[parent_taxid].append(taxid)
    
    # Find all descendants of VIRUS_ROOT_TAXID using BFS
    virus_taxids = set([VIRUS_ROOT_TAXID])
    queue = [VIRUS_ROOT_TAXID]
    
    while queue:
        current_taxid = queue.pop(0)
        if current_taxid in parent_to_children:
            for child_taxid in parent_to_children[current_taxid]:
                if child_taxid not in virus_taxids:
                    virus_taxids.add(child_taxid)
                    queue.append(child_taxid)
    
    logger.info(f"Found {len(virus_taxids):,} virus taxids (including all descendants of taxid {VIRUS_ROOT_TAXID})")
    return virus_taxids

def is_complete_genome_accession(accession: str) -> bool:
    """
    Check if accession likely represents a complete genome based on prefix.
    NC_ = RefSeq complete genomes/chromosomes
    AC_ = GenBank complete genomes
    Note: This is not perfect - some complete genomes use other prefixes.
    """
    accession_upper = accession.upper()
    # Common prefixes for complete genomes
    complete_prefixes = ['NC_', 'AC_']
    return any(accession_upper.startswith(prefix) for prefix in complete_prefixes)

def build_taxonomy_database_from_ncbi_files(database_dir: Path, accession2taxid_file: Path, names_file: Path, nodes_file: Path = None, viruses_only: bool = True, complete_genomes_only: bool = False, reference_accessions: set = None) -> dict:
    """
    Build comprehensive taxonomy database from NCBI dump files.
    If viruses_only=True, only includes virus accessions (much smaller file size).
    Returns dict mapping accession -> organism_name.
    """
    logger.info("Building taxonomy database from NCBI files...")
    if viruses_only:
        logger.info("Filtering to VIRUSES ONLY (to reduce file size)")
    if complete_genomes_only:
        logger.info("Filtering to COMPLETE GENOMES ONLY (NC_ and AC_ prefixes)")
        logger.warning("Note: Complete genome filtering is approximate - some complete genomes may use other prefixes")
    if reference_accessions:
        logger.info(f"Filtering to ACCESSIONS IN REFERENCE DATABASE ONLY ({len(reference_accessions):,} accessions)")
    
    # Step 0: Get virus taxids if filtering to viruses only
    virus_taxids = set()
    if viruses_only:
        if nodes_file and nodes_file.exists():
            virus_taxids = get_virus_taxids(nodes_file)
        else:
            logger.warning("nodes.dmp not found - cannot filter to viruses. Building full database.")
            viruses_only = False
    
    # Step 1: Parse names.dmp to create taxid -> organism_name mapping
    logger.info("Parsing names.dmp (taxid -> organism name)...")
    taxid_to_organism = {}
    taxid_to_name = {}  # Store all names for species lookup
    with open(names_file, 'r') as f:
        for line in f:
            # names.dmp format: taxid | name | unique name | name class |
            # Fields are separated by \t|\t and line ends with \t|
            line = line.rstrip('\n\r')
            if not line:
                continue
            # Remove trailing \t| if present
            if line.endswith('\t|'):
                line = line[:-2]
            # Split by tab-pipe-tab delimiter
            parts = [p.strip() for p in line.split('\t|\t')]
            
            if len(parts) >= 4:
                taxid = parts[0]
                name = parts[1]
                name_type = parts[3]
                # Only keep scientific names
                # If filtering to viruses, only keep virus taxids
                if name_type == "scientific name":
                    taxid_to_name[taxid] = name  # Store all for species lookup
                    if not viruses_only or taxid in virus_taxids:
                        taxid_to_organism[taxid] = name
    
    logger.info(f"Loaded {len(taxid_to_organism):,} taxid -> organism mappings")
    
    # Step 1.5: Build taxid -> species mapping using nodes.dmp
    logger.info("Building taxid -> species mapping from nodes.dmp...")
    taxid_to_species = {}
    if nodes_file and nodes_file.exists():
        # Parse nodes.dmp to build parent mapping and rank mapping
        taxid_to_parent = {}
        taxid_to_rank = {}
        with open(nodes_file, 'r') as f:
            for line in f:
                parts = [p.strip() for p in line.rstrip('\n\r').rstrip('\t|').split('\t|\t')]
                if len(parts) >= 3:
                    child_taxid = parts[0]
                    parent_taxid = parts[1]
                    rank = parts[2]
                    taxid_to_parent[child_taxid] = parent_taxid
                    taxid_to_rank[child_taxid] = rank
        
        # For each taxid, walk up the tree to find species-level taxid
        for taxid in taxid_to_organism.keys():
            current_taxid = taxid
            visited = set()
            
            while current_taxid and current_taxid not in visited:
                visited.add(current_taxid)
                rank = taxid_to_rank.get(current_taxid, "").lower()
                
                if rank == "species":
                    species_name = taxid_to_name.get(current_taxid, "")
                    if species_name:
                        taxid_to_species[taxid] = species_name
                    break
                
                # Move to parent
                current_taxid = taxid_to_parent.get(current_taxid)
                if not current_taxid or current_taxid == "1":  # Root
                    break
        
        logger.info(f"Built {len(taxid_to_species):,} taxid -> species mappings")
    
    # Step 2: Parse nucl_gb.accession2taxid.gz to create accession -> taxid mapping
    logger.info("Parsing nucl_gb.accession2taxid.gz (accession -> taxid -> organism/species)...")
    accession_to_organism = {}
    accession_to_species = {}  # New: store species mapping
    line_count = 0
    virus_count = 0
    skipped_count = 0
    with gzip.open(accession2taxid_file, 'rt') as f:
        next(f)  # Skip header line
        for line in f:
            line_count += 1
            if line_count % 1000000 == 0:
                logger.info(f"  Processed {line_count:,} lines... (kept {virus_count:,} virus accessions)")
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                accession = parts[0].strip()
                accession_versioned = parts[1].strip()
                taxid = parts[2].strip()
                
                # If filtering to viruses, skip non-virus taxids
                if viruses_only and taxid not in virus_taxids:
                    skipped_count += 1
                    continue
                
                # If filtering to complete genomes only, skip non-complete genome accessions
                if complete_genomes_only:
                    # Check both versioned and unversioned accessions
                    if not (is_complete_genome_accession(accession) or is_complete_genome_accession(accession_versioned)):
                        skipped_count += 1
                        continue
                
                # If filtering to reference database accessions only, skip accessions not in reference database
                if reference_accessions:
                    # Check both versioned and unversioned accessions
                    if not (accession in reference_accessions or accession_versioned in reference_accessions or 
                            ('.' in accession and accession.split('.')[0] in reference_accessions) or
                            ('.' in accession_versioned and accession_versioned.split('.')[0] in reference_accessions)):
                        skipped_count += 1
                        continue
                
                # Map both versioned and unversioned accessions
                if taxid in taxid_to_organism:
                    organism = taxid_to_organism[taxid]
                    accession_to_organism[accession] = organism
                    accession_to_organism[accession_versioned] = organism
                    
                    # Also store species if available
                    if taxid in taxid_to_species:
                        species = taxid_to_species[taxid]
                        accession_to_species[accession] = species
                        accession_to_species[accession_versioned] = species
                    
                    virus_count += 1
    
    if viruses_only:
        logger.info(f"Filtered to viruses: kept {virus_count:,} accessions, skipped {skipped_count:,} non-virus accessions")
    if complete_genomes_only:
        logger.info(f"Filtered to complete genomes: kept accessions with NC_ or AC_ prefixes")
    if reference_accessions:
        logger.info(f"Filtered to reference database: kept accessions present in reference database")
    logger.info(f"Built database with {len(accession_to_organism):,} accession -> organism mappings")
    logger.info(f"Built database with {len(accession_to_species):,} accession -> species mappings")
    return accession_to_organism, accession_to_species

def build_taxonomy_database_from_ncbi_download(database_dir: Path, database_fasta: Path = None, viruses_only: bool = True, complete_genomes_only: bool = False) -> dict:
    """
    Download NCBI taxonomy files and build comprehensive database.
    If viruses_only=True, only includes virus accessions (much smaller file size).
    If complete_genomes_only=True, only includes complete genomes (NC_ and AC_ prefixes).
    If database_fasta is provided, only includes accessions present in that FASTA file.
    Saves the final database to SQLite for fast lookups.
    """
    ncbi_dir = get_ncbi_files_dir(database_dir)
    # Use nucl_gb.accession2taxid.gz (GenBank nucleotide, includes RefSeq)
    accession2taxid_file = ncbi_dir / "nucl_gb.accession2taxid.gz"
    names_file = ncbi_dir / "names.dmp"
    nodes_file = ncbi_dir / "nodes.dmp"
    database_path = get_taxonomy_database_path(database_dir)
    
    # Extract accessions from reference database if provided
    reference_accessions = None
    if database_fasta and Path(database_fasta).exists():
        logger.info(f"Extracting accessions from reference database: {Path(database_fasta).name}...")
        reference_accessions = extract_accessions_from_fasta(Path(database_fasta))
        logger.info(f"Found {len(reference_accessions):,} accessions in reference database")
    
    # Download files if needed
    download_ncbi_taxonomy_files(database_dir, force_download=False)
    
    # Build database (filtered to viruses only by default, optionally complete genomes only, optionally reference database only)
    accession_to_organism, accession_to_species = build_taxonomy_database_from_ncbi_files(
        database_dir, accession2taxid_file, names_file, nodes_file, 
        viruses_only=viruses_only, complete_genomes_only=complete_genomes_only,
        reference_accessions=reference_accessions
    )
    
    # Save to SQLite (much faster than JSON for large databases)
    database_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Ensure we're using .db extension (not .json)
    if database_path.suffix != '.db':
        database_path = database_path.with_suffix('.db')
    
    logger.info(f"Saving comprehensive taxonomy database to {database_path.name} (SQLite format)...")
    logger.info(f"Full path: {database_path}")
    
    # Remove old database if it exists (both .db and .json)
    if database_path.exists():
        database_path.unlink()
    json_path = database_path.with_suffix('.json')
    if json_path.exists():
        logger.info(f"Removing old JSON database: {json_path.name}")
        json_path.unlink()
    
    # Create SQLite database with indexed table
    try:
        conn = sqlite3.connect(str(database_path))
        cursor = conn.cursor()
        
        # Create table with accession as primary key (automatically indexed)
        # Include both organism and species columns
        cursor.execute('''
            CREATE TABLE accession_to_organism (
                accession TEXT PRIMARY KEY,
                organism TEXT NOT NULL,
                species TEXT
            )
        ''')
        
        # Create index on accession for faster lookups (though PRIMARY KEY is already indexed)
        # Also create index for versioned accessions (without version number)
        cursor.execute('''
            CREATE INDEX idx_accession_base ON accession_to_organism(substr(accession, 1, instr(accession || '.', '.') - 1))
        ''')
        
        # Insert all mappings in batches for speed (include species if available)
        batch_size = 100000
        items = []
        for accession, organism in accession_to_organism.items():
            species = accession_to_species.get(accession, None)
            items.append((accession, organism, species))
        
        total_batches = (len(items) + batch_size - 1) // batch_size
        
        for batch_num in range(total_batches):
            batch = items[batch_num * batch_size:(batch_num + 1) * batch_size]
            cursor.executemany('INSERT OR REPLACE INTO accession_to_organism (accession, organism, species) VALUES (?, ?, ?)', batch)
            if (batch_num + 1) % 10 == 0:
                logger.info(f"  Inserted batch {batch_num + 1}/{total_batches} ({len(batch)} entries)...")
        
        conn.commit()
        conn.close()
        
        # Get file size
        file_size_mb = database_path.stat().st_size / (1024 * 1024)
        logger.info(f"Saved {len(accession_to_organism):,} accession -> organism mappings and {len(accession_to_species):,} accession -> species mappings to SQLite database ({file_size_mb:.1f} MB)")
    except Exception as e:
        logger.error(f"Failed to save to SQLite database: {e}")
        logger.error(f"Falling back to JSON format (slower)...")
        # Fallback to JSON if SQLite fails
        json_path = database_path.with_suffix('.json')
        # Store both organism and species in JSON format
        json_data = {}
        for accession, organism in accession_to_organism.items():
            species = accession_to_species.get(accession)
            json_data[accession] = {'organism': organism, 'species': species}
        with open(json_path, 'w') as f:
            json.dump(json_data, f)
        file_size_mb = json_path.stat().st_size / (1024 * 1024)
        logger.info(f"Saved {len(accession_to_organism):,} mappings to JSON database ({file_size_mb:.1f} MB)")
        raise  # Re-raise to indicate there was an issue
    
    return accession_to_organism, accession_to_species

def fetch_organism_from_ncbi(accession: str, cache_path: Path = None, database_path: Path = None, max_retries: int = 3) -> str:
    """
    Fetch organism name from NCBI for a given accession.
    First checks comprehensive database, then falls back to API calls.
    """
    # Remove version from accession if present (e.g., "NC_001474.2" -> "NC_001474")
    accession_base = accession.split('.')[0] if '.' in accession else accession
    
    # Step 1: Check comprehensive database first (SQLite or JSON)
    if database_path and database_path.exists():
        try:
            if database_path.suffix == '.db':
                # SQLite database
                import sqlite3
                conn = sqlite3.connect(str(database_path))
                cursor = conn.cursor()
                # Try both versioned and unversioned
                cursor.execute('SELECT organism FROM accession_to_organism WHERE accession = ?', (accession,))
                result = cursor.fetchone()
                if result:
                    conn.close()
                    return result[0]
                cursor.execute('SELECT organism FROM accession_to_organism WHERE accession = ?', (accession_base,))
                result = cursor.fetchone()
                if result:
                    conn.close()
                    return result[0]
                conn.close()
            else:
                # JSON database (legacy)
                with open(database_path, 'r') as f:
                    database = json.load(f)
                    # Try both versioned and unversioned
                    if accession in database:
                        return database[accession]
                    if accession_base in database:
                        return database[accession_base]
        except Exception as e:
            logger.debug(f"Could not read comprehensive database: {e}")
    
    # Step 2: Check individual cache
    if cache_path:
        cache = load_taxonomy_cache(cache_path)
        if accession in cache:
            return cache[accession]
        if accession_base in cache:
            return cache[accession_base]
    
    # Step 3: Fetch from NCBI API (with rate limiting)
    for attempt in range(max_retries):
        try:
            time.sleep(0.1)  # Rate limiting
            
            # Use esearch to get UID
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={accession_base}&retmode=json"
            with urllib.request.urlopen(search_url, timeout=30) as response:
                data = json.loads(response.read())
                if 'esearchresult' not in data or not data['esearchresult'].get('idlist'):
                    logger.debug(f"No UID found for {accession}")
                    return ""
                uid = data['esearchresult']['idlist'][0]
            
            # Use esummary to get organism
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={uid}&retmode=json"
            with urllib.request.urlopen(summary_url, timeout=30) as response:
                data = json.loads(response.read())
                if 'result' in data and uid in data['result']:
                    # Try lowercase 'organism' first (correct field name), fallback to 'Organism' for compatibility
                    organism = data['result'][uid].get('organism', '') or data['result'][uid].get('Organism', '')
                    if organism:
                        # Save to cache
                        if cache_path:
                            cache = load_taxonomy_cache(cache_path)
                            cache[accession] = organism
                            cache[accession_base] = organism
                            save_taxonomy_cache(cache_path, cache)
                        logger.debug(f"Fetched organism '{organism}' from NCBI API for {accession}")
                        return organism
            logger.debug(f"No organism found in esummary response for {accession}")
            return ""
        except Exception as e:
            if attempt < max_retries - 1:
                logger.debug(f"Attempt {attempt + 1} failed for {accession}: {e}")
                time.sleep(1)
            else:
                logger.warning(f"Failed to fetch organism for {accession} after {max_retries} attempts: {e}")
    
    return ""

def enrich_stats_with_organism(stats_list: list, cache_path: Path = None, database_path: Path = None, nodes_file: Path = None, names_file: Path = None) -> list:
    """Add organism names to stats list by fetching from NCBI."""
    if not stats_list:
        return stats_list
    
    logger.info(f"Fetching ORGANISM names for {len(stats_list)} references...")
    
    # Load database ONCE (SQLite is much faster than JSON for large databases)
    database = {}
    if database_path and database_path.exists():
        try:
            logger.info(f"Loading comprehensive taxonomy database from {database_path.name}...")
            start_time = time.time()
            
            # Try SQLite first (faster)
            if database_path.suffix == '.db':
                conn = sqlite3.connect(str(database_path))
                conn.row_factory = sqlite3.Row
                cursor = conn.cursor()
                
                # Get count first
                cursor.execute('SELECT COUNT(*) as count FROM accession_to_organism')
                count = cursor.fetchone()['count']
                
                # Load all into memory for fast lookups (SQLite is fast but in-memory dict is even faster)
                # Check if species column exists
                cursor.execute("PRAGMA table_info(accession_to_organism)")
                columns = [col[1] for col in cursor.fetchall()]
                has_species = 'species' in columns
                
                if has_species:
                    cursor.execute('SELECT accession, organism, species FROM accession_to_organism')
                    database = {row['accession']: {'organism': row['organism'], 'species': row['species']} for row in cursor}
                else:
                    cursor.execute('SELECT accession, organism FROM accession_to_organism')
                    database = {row['accession']: {'organism': row['organism'], 'species': None} for row in cursor}
                
                conn.close()
                load_time = time.time() - start_time
                logger.info(f"Loaded {len(database):,} accession -> organism mappings from SQLite database in {load_time:.2f}s")
            else:
                # Fallback to JSON (legacy support)
                with open(database_path, 'r') as f:
                    database = json.load(f)
                load_time = time.time() - start_time
                logger.info(f"Loaded {len(database):,} accession -> organism mappings from JSON database in {load_time:.2f}s")
        except Exception as e:
            logger.warning(f"Could not load comprehensive database: {e}")
            # Try JSON fallback if SQLite fails
            json_path = database_path.with_suffix('.json') if database_path.suffix == '.db' else database_path
            if json_path.exists():
                try:
                    logger.info(f"Trying JSON fallback: {json_path.name}...")
                    with open(json_path, 'r') as f:
                        database = json.load(f)
                    logger.info(f"Loaded {len(database):,} accession -> organism mappings from JSON fallback")
                except Exception as e2:
                    logger.warning(f"Could not load JSON fallback either: {e2}")
    
    # Load cache ONCE
    cache = {}
    if cache_path:
        cache = load_taxonomy_cache(cache_path)
        if cache:
            logger.info(f"Loaded {len(cache)} entries from individual cache")
    
    cache_hits = 0
    database_hits = 0
    api_calls = 0
    
    for i, stat in enumerate(stats_list):
        accession = stat.get("accession", "")
        if not accession:
            stat["organism"] = ""
            continue
        
        accession_base = accession.split('.')[0] if '.' in accession else accession
        organism = ""
        
        # Step 1: Check individual cache first (fastest, if available)
        # Note: Cache is only used if database doesn't have it, but we check cache first
        # to avoid loading database for already-cached entries (though database is usually already loaded)
        if not organism and cache:
            if accession in cache:
                organism = cache[accession]
                cache_hits += 1
            elif accession_base in cache:
                organism = cache[accession_base]
                cache_hits += 1
        
        # Step 2: Check comprehensive database (primary source, usually has everything)
        if not organism and database:
            db_entry = None
            if accession in database:
                db_entry = database[accession]
                database_hits += 1
            elif accession_base in database:
                db_entry = database[accession_base]
                database_hits += 1
            
            if db_entry:
                # Handle both old format (string) and new format (dict)
                if isinstance(db_entry, dict):
                    organism = db_entry.get('organism', '')
                    # Also get species if available
                    species = db_entry.get('species')
                    if species:
                        stat["viral_species"] = species
                else:
                    organism = db_entry
        
        # Step 3: Fetch from NCBI API only if not found in database or cache
        if not organism:
            organism = fetch_organism_from_ncbi_with_loaded_data(
                accession, database, cache, cache_path
            )
            if organism:
                api_calls += 1
        
        stat["organism"] = organism if organism else ""  # Ensure it's set, even if empty
        
        # Extract viral species - use from database if available, otherwise try offline lookup, then API, then fallback to parsing
        # For custom accessions, organism might be empty from FASTA header, so we need to fetch from NCBI
        if "viral_species" not in stat or not stat.get("viral_species"):
            if organism:
                viral_species = extract_viral_species(organism, accession=accession, database_path=database_path, cache_path=cache_path, nodes_file=nodes_file, names_file=names_file)
                stat["viral_species"] = viral_species
            elif accession:
                # If no organism found, try fetching species directly from NCBI using accession
                # This helps for custom accessions where organism wasn't in the FASTA header
                viral_species = fetch_species_from_ncbi(accession, database_path=database_path, cache_path=cache_path, nodes_file=nodes_file, names_file=names_file)
                if viral_species:
                    stat["viral_species"] = viral_species
                    # Also try to get organism from NCBI API if we still don't have it
                    if not organism:
                        organism = fetch_organism_from_ncbi_with_loaded_data(accession, database, cache, cache_path)
                        if organism:
                            stat["organism"] = organism
                            api_calls += 1
                else:
                    stat["viral_species"] = ""
            else:
                stat["viral_species"] = ""
        
        if (i + 1) % 10 == 0:
            logger.info(f"  Processed {i + 1}/{len(stats_list)} organisms... (database: {database_hits}, cache: {cache_hits}, API: {api_calls})")
    
    logger.info(f"Completed fetching ORGANISM names and extracting viral species (database: {database_hits}, cache: {cache_hits}, API: {api_calls})")
    return stats_list

def fetch_organism_from_ncbi_with_loaded_data(accession: str, database: dict = None, cache: dict = None, cache_path: Path = None, max_retries: int = 3) -> str:
    """
    Fetch organism name from NCBI API when not found in pre-loaded database/cache.
    This version accepts already-loaded database and cache dictionaries for efficiency.
    """
    accession_base = accession.split('.')[0] if '.' in accession else accession
    
    # Double-check database and cache (in case they weren't checked before)
    if database:
        if accession in database:
            return database[accession]
        if accession_base in database:
            return database[accession_base]
    
    if cache:
        if accession in cache:
            return cache[accession]
        if accession_base in cache:
            return cache[accession_base]
    
    # Fetch from NCBI API (with rate limiting)
    for attempt in range(max_retries):
        try:
            time.sleep(0.1)  # Rate limiting
            
            # Use esearch to get UID
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={accession_base}&retmode=json"
            with urllib.request.urlopen(search_url) as response:
                data = json.loads(response.read())
                if 'esearchresult' not in data or not data['esearchresult'].get('idlist'):
                    return ""
                uid = data['esearchresult']['idlist'][0]
            
            # Use esummary to get organism
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={uid}&retmode=json"
            with urllib.request.urlopen(summary_url, timeout=30) as response:
                data = json.loads(response.read())
                if 'result' in data and uid in data['result']:
                    # Try lowercase 'organism' first (correct field name), fallback to 'Organism' for compatibility
                    organism = data['result'][uid].get('organism', '') or data['result'][uid].get('Organism', '')
                    if organism:
                        # Save to cache
                        if cache_path:
                            cache_to_save = load_taxonomy_cache(cache_path)
                            cache_to_save[accession] = organism
                            cache_to_save[accession_base] = organism
                            save_taxonomy_cache(cache_path, cache_to_save)
                        return organism
            return ""
        except Exception as e:
            if attempt < max_retries - 1:
                logger.debug(f"Attempt {attempt + 1} failed for {accession}: {e}")
                time.sleep(1)
            else:
                logger.warning(f"Failed to fetch organism for {accession} after {max_retries} attempts: {e}")
    
    return ""

def get_species_from_taxid(taxid: str, nodes_file: Path = None, names_file: Path = None, taxid_to_species_cache: dict = None) -> str:
    """
    Get species name from taxid using local taxonomy files (offline).
    Walks up the taxonomy tree to find the species-level taxid.
    
    Args:
        taxid: Taxonomy ID
        nodes_file: Path to nodes.dmp file
        names_file: Path to names.dmp file
        taxid_to_species_cache: Optional cache dict (taxid -> species) to avoid re-parsing
    
    Returns:
        Species name or empty string if not found
    """
    if not taxid:
        return ""
    
    # Use cache if provided
    if taxid_to_species_cache is not None and taxid in taxid_to_species_cache:
        return taxid_to_species_cache[taxid]
    
    if not nodes_file or not nodes_file.exists() or not names_file or not names_file.exists():
        return ""
    
    # Parse nodes.dmp to build parent mapping: taxid -> parent_taxid, rank
    taxid_to_parent = {}
    taxid_to_rank = {}
    with open(nodes_file, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip('\n\r').rstrip('\t|').split('\t|\t')]
            if len(parts) >= 3:
                child_taxid = parts[0]
                parent_taxid = parts[1]
                rank = parts[2]
                taxid_to_parent[child_taxid] = parent_taxid
                taxid_to_rank[child_taxid] = rank
    
    # Parse names.dmp to get taxid -> name mapping
    taxid_to_name = {}
    with open(names_file, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip('\n\r').rstrip('\t|').split('\t|\t')]
            if len(parts) >= 4:
                tid = parts[0]
                name = parts[1]
                name_type = parts[3]
                if name_type == "scientific name":
                    taxid_to_name[tid] = name
    
    # Walk up the tree to find species-level taxid
    current_taxid = taxid
    visited = set()
    
    while current_taxid and current_taxid not in visited:
        visited.add(current_taxid)
        
        # Check if current taxid is at species rank
        rank = taxid_to_rank.get(current_taxid, "").lower()
        if rank == "species":
            species_name = taxid_to_name.get(current_taxid, "")
            if species_name:
                # Update cache if provided
                if taxid_to_species_cache is not None:
                    taxid_to_species_cache[taxid] = species_name
                return species_name
        
        # Move to parent
        current_taxid = taxid_to_parent.get(current_taxid)
        if not current_taxid or current_taxid == "1":  # Root of tree
            break
    
    return ""

def fetch_species_from_ncbi(accession: str, database_path: Path = None, cache_path: Path = None, nodes_file: Path = None, names_file: Path = None, max_retries: int = 3) -> str:
    """
    Fetch the actual species name from NCBI taxonomy using the accession.
    First checks local database/cache, then offline taxonomy files, finally API as last resort.
    
    Args:
        accession: NCBI accession number
        database_path: Path to comprehensive taxonomy database (SQLite or JSON)
        cache_path: Path to individual cache file
        nodes_file: Path to nodes.dmp file (for offline lookup)
        names_file: Path to names.dmp file (for offline lookup)
        max_retries: Maximum number of retry attempts for API calls
    
    Returns:
        Species name (e.g., "Measles morbillivirus") or empty string if not found
    """
    if not accession:
        return ""
    
    accession_base = accession.split('.')[0] if '.' in accession else accession
    
    # Step 1: Check comprehensive database (if it includes species)
    if database_path and database_path.exists():
        try:
            if database_path.suffix == '.db':
                # SQLite database
                conn = sqlite3.connect(str(database_path))
                cursor = conn.cursor()
                # Check if species column exists
                cursor.execute("PRAGMA table_info(accession_to_organism)")
                columns = [col[1] for col in cursor.fetchall()]
                has_species = 'species' in columns
                
                if has_species:
                    cursor.execute('SELECT species FROM accession_to_organism WHERE accession = ?', (accession,))
                    result = cursor.fetchone()
                    if result and result[0]:
                        conn.close()
                        return result[0]
                    
                    # Try unversioned
                    cursor.execute('SELECT species FROM accession_to_organism WHERE accession = ?', (accession_base,))
                    result = cursor.fetchone()
                    if result and result[0]:
                        conn.close()
                        return result[0]
                else:
                    # Database doesn't have species yet - try to get taxid and look up species
                    cursor.execute('SELECT organism FROM accession_to_organism WHERE accession = ?', (accession,))
                    result = cursor.fetchone()
                    if not result:
                        cursor.execute('SELECT organism FROM accession_to_organism WHERE accession = ?', (accession_base,))
                        result = cursor.fetchone()
                    
                    # If we have organism but not species, we need taxid - skip for now, will use offline files
                conn.close()
            else:
                # JSON database
                with open(database_path, 'r') as f:
                    database = json.load(f)
                    if accession in database:
                        entry = database[accession]
                        if isinstance(entry, dict) and 'species' in entry:
                            return entry['species']
                    if accession_base in database:
                        entry = database[accession_base]
                        if isinstance(entry, dict) and 'species' in entry:
                            return entry['species']
        except Exception as e:
            logger.debug(f"Could not read species from database: {e}")
    
    # Step 2: Check individual cache
    if cache_path:
        cache = load_taxonomy_cache(cache_path)
        if accession in cache:
            entry = cache[accession]
            if isinstance(entry, dict) and 'species' in entry:
                return entry['species']
        if accession_base in cache:
            entry = cache[accession_base]
            if isinstance(entry, dict) and 'species' in entry:
                return entry['species']
    
    # Step 3: Try offline lookup using taxonomy files (if available)
    # First, we need to get taxid from accession
    taxid = None
    if database_path and database_path.exists():
        try:
            # Try to get taxid from accession2taxid mapping (if we have it stored)
            # For now, we'll need to parse it from the downloaded file or use API
            pass
        except:
            pass
    
    # Step 4: Fall back to API (only if offline methods failed)
    for attempt in range(max_retries):
        try:
            time.sleep(0.1)  # Rate limiting
            
            # Get UID and TaxId from nuccore
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={accession_base}&retmode=json"
            with urllib.request.urlopen(search_url) as response:
                data = json.loads(response.read())
                if 'esearchresult' not in data or not data['esearchresult'].get('idlist'):
                    return ""
                uid = data['esearchresult']['idlist'][0]
            
            # Get summary with TaxId
            summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={uid}&retmode=json"
            with urllib.request.urlopen(summary_url) as response:
                data = json.loads(response.read())
                if 'result' not in data or uid not in data['result']:
                    return ""
                
                taxid = data['result'][uid].get('TaxId', '')
                if not taxid:
                    return ""
            
            # Try offline lookup first if files are available
            if nodes_file and names_file and nodes_file.exists() and names_file.exists():
                species = get_species_from_taxid(str(taxid), nodes_file, names_file)
                if species:
                    # Cache it
                    if cache_path:
                        cache = load_taxonomy_cache(cache_path)
                        if accession not in cache:
                            cache[accession] = {}
                        if not isinstance(cache[accession], dict):
                            cache[accession] = {'organism': cache[accession]}
                        cache[accession]['species'] = species
                        save_taxonomy_cache(cache_path, cache)
                    return species
            
            # Fallback to API taxonomy lookup
            tax_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml"
            with urllib.request.urlopen(tax_url) as response:
                import xml.etree.ElementTree as ET
                tree = ET.parse(response)
                root = tree.getroot()
                
                for taxon in root.findall('.//Taxon'):
                    for name in taxon.findall('.//ScientificName'):
                        rank_elem = taxon.find('.//Rank')
                        if rank_elem is not None and rank_elem.text == 'species':
                            species = name.text
                            # Cache it
                            if cache_path:
                                cache = load_taxonomy_cache(cache_path)
                                if accession not in cache:
                                    cache[accession] = {}
                                if not isinstance(cache[accession], dict):
                                    cache[accession] = {'organism': cache[accession]}
                                cache[accession]['species'] = species
                                save_taxonomy_cache(cache_path, cache)
                            return species
            
            return ""
        except Exception as e:
            if attempt < max_retries - 1:
                logger.debug(f"Attempt {attempt + 1} failed to fetch species for {accession}: {e}")
                time.sleep(1)
            else:
                logger.debug(f"Failed to fetch species for {accession} after {max_retries} attempts: {e}")
    
    return ""

def extract_viral_species(organism: str, accession: str = None, database_path: Path = None, cache_path: Path = None, nodes_file: Path = None, names_file: Path = None) -> str:
    """
    Extract viral species from organism name or fetch from NCBI taxonomy.
    For viruses, the organism name may include genotype/strain info (e.g., "Measles virus genotype B3"),
    but the actual species is "Measles morbillivirus".
    
    This function first tries to fetch the actual species from NCBI taxonomy if accession is provided.
    Otherwise, it attempts to parse the organism name to remove strain/genotype information.
    
    Examples:
        "Measles virus genotype B3" -> "Measles morbillivirus" (from NCBI taxonomy)
        "Canine morbillivirus" -> "Canine morbillivirus"
        "Lassa mammarenavirus" -> "Lassa mammarenavirus"
        "Influenza A virus (A/H1N1)" -> "Influenza A virus" (or actual species from NCBI)
    """
    # Try to fetch actual species from NCBI if accession is provided
    if accession:
        species = fetch_species_from_ncbi(accession, database_path=database_path, cache_path=cache_path, nodes_file=nodes_file, names_file=names_file)
        if species:
            return species
    
    # Fallback: try to parse organism name
    if not organism:
        return ""
    
    organism = organism.strip()
    
    # Remove common strain/isolate indicators in parentheses
    import re
    
    # Pattern to match strain info in parentheses at the end
    # e.g., "Influenza A virus (A/H1N1)" -> "Influenza A virus"
    strain_pattern = r'\s*\([^)]*\)\s*$'
    species = re.sub(strain_pattern, '', organism)
    
    # Remove common prefixes/suffixes that indicate genotypes/strains
    # e.g., "Measles virus genotype B3" -> try to extract just "Measles virus" or "Measles morbillivirus"
    # Remove "genotype X", "strain X", "isolate X", "variant X"
    strain_suffixes = [
        r'\s+genotype\s+[A-Z0-9]+.*$',
        r'\s+strain\s+.*$',
        r'\s+isolate\s+.*$',
        r'\s+variant\s+.*$',
    ]
    for suffix_pattern in strain_suffixes:
        species = re.sub(suffix_pattern, '', species, flags=re.IGNORECASE)
    
    # For "Measles virus genotype B3", try to convert to "Measles morbillivirus"
    # This is a heuristic - ideally we should use NCBI taxonomy
    if "measles virus" in species.lower() and "morbillivirus" not in species.lower():
        species = "Measles morbillivirus"
    
    return species.strip()

def deduplicate_by_species(stats_list: list) -> list:
    """Deduplicate by species name (from header), preferring RefSeq over GenBank."""
    if not stats_list:
        return []
    
    species_dict = {}
    for stat in stats_list:
        desc = stat.get("description", "")
        species = extract_species_from_header(desc)
        if not species:
            continue
        
        if species not in species_dict:
            species_dict[species] = stat
        else:
            # Prefer RefSeq over GenBank
            current_is_refseq = is_refseq(stat["description"])
            existing_is_refseq = is_refseq(species_dict[species]["description"])
            
            if current_is_refseq and not existing_is_refseq:
                species_dict[species] = stat
            elif current_is_refseq == existing_is_refseq:
                # If both same type, prefer one with more mapped reads
                if stat["mapped_reads"] > species_dict[species]["mapped_reads"]:
                    species_dict[species] = stat
    
    return list(species_dict.values())

def deduplicate_by_organism(stats_list: list, segment_database_path: Path = None) -> list:
    """
    Deduplicate by viral species (instead of organism), preferring RefSeq over GenBank.
    For segmented viruses, keeps the best reference for EACH segment.
    
    Args:
        stats_list: List of statistics dictionaries
        segment_database_path: Path to segment database (optional, for accurate segment detection)
    """
    if not stats_list:
        return []
    
    logger.info(f"Detecting segments for {len(stats_list)} references and deduplicating by viral species (checking database first, then NCBI if needed)...")
    
    # For segmented viruses, we need to group by viral_species + segment
    # For non-segmented, we group by viral_species only
    species_segment_dict = {}
    segment_db_hits = 0
    segment_api_calls = 0
    segment_fallback = 0
    
    # Debug: Track viral_species values to see if there are duplicates
    species_counts = {}
    
    for i, stat in enumerate(stats_list):
        # Use viral_species for deduplication (fallback to organism if not available)
        viral_species = stat.get("viral_species", "").strip()
        if not viral_species:
            # Fallback: extract from organism if not already set
            organism = stat.get("organism", "").strip()
            accession = stat.get("accession", "")
            if organism:
                # Note: We don't have taxonomy file paths here, so it will use API or parsing
                viral_species = extract_viral_species(organism, accession=accession)
                stat["viral_species"] = viral_species
            else:
                # If no organism/species, still keep it (might be important)
                # Use a unique key so it doesn't get deduplicated away
                key = f"NO_SPECIES||{stat.get('accession', 'unknown')}"
                if key not in species_segment_dict:
                    species_segment_dict[key] = stat
                continue
        
        # Normalize viral_species for consistent matching (case-insensitive, trimmed)
        viral_species = viral_species.strip()
        if not viral_species:
            # If still empty after normalization, use accession as fallback
            key = f"NO_SPECIES||{stat.get('accession', 'unknown')}"
            if key not in species_segment_dict:
                species_segment_dict[key] = stat
            continue
        
        # Extract segment information using accession (more reliable) or description (fallback)
        description = stat.get("description", "")
        accession = stat.get("accession", "")
        
        # Check database first (fast)
        segment = get_segment_from_database(accession, segment_database_path)
        if segment:
            segment_db_hits += 1
        else:
            # Try fetching from NCBI (slow, but accurate)
            segment = fetch_segment_from_ncbi(accession, segment_database_path)
            if segment:
                segment_api_calls += 1
            else:
                # Fall back to description parsing (fast, but less reliable)
                segment = extract_segment_from_description(description)
                if segment:
                    segment_fallback += 1
                    # Save parsed segment to database for future use
                    if segment_database_path:
                        save_segment_to_database(accession, segment, segment_database_path)
        
        # Log progress every 20 references
        if (i + 1) % 20 == 0:
            logger.info(f"  Processed {i + 1}/{len(stats_list)} references for segment detection... (DB: {segment_db_hits}, API: {segment_api_calls}, fallback: {segment_fallback})")
        
        # Create key: viral_species + segment (or just viral_species if not segmented)
        # Normalize key to ensure consistent matching (case-insensitive)
        viral_species_normalized = viral_species.strip().lower()
        if segment:
            segment_normalized = segment.strip().upper()  # Segments are usually uppercase (L, S, etc.)
            key = f"{viral_species_normalized}||SEGMENT||{segment_normalized}"
        else:
            key = f"{viral_species_normalized}||NO_SEGMENT"
        
        # Debug: Track species counts
        if viral_species_normalized not in species_counts:
            species_counts[viral_species_normalized] = 0
        species_counts[viral_species_normalized] += 1
        
        if key not in species_segment_dict:
            species_segment_dict[key] = stat
        else:
            # Prefer RefSeq over GenBank
            current_is_refseq = is_refseq(stat["description"])
            existing_is_refseq = is_refseq(species_segment_dict[key]["description"])
            
            current_reads = stat.get("mapped_reads", 0)
            existing_reads = species_segment_dict[key].get("mapped_reads", 0)
            
            # Debug: Log when Lassa virus references are being compared
            organism = stat.get("organism", "")
            if "lassa" in viral_species.lower() or "lassa" in organism.lower() or "lassa" in stat.get("description", "").lower():
                logger.debug(f"DEBUG: Comparing Lassa references for key {key}:")
                logger.debug(f"  Current: {stat.get('accession')} ({current_reads} reads, RefSeq: {current_is_refseq}, species: {viral_species})")
                logger.debug(f"  Existing: {species_segment_dict[key].get('accession')} ({existing_reads} reads, RefSeq: {existing_is_refseq}, species: {species_segment_dict[key].get('viral_species', 'N/A')})")
            
            # Never replace a reference with reads with one that has 0 reads
            if existing_reads > 0 and current_reads == 0:
                # Keep existing (has reads)
                if "lassa" in viral_species.lower() or "lassa" in organism.lower():
                    logger.warning(f"DEBUG: Keeping existing Lassa reference {species_segment_dict[key].get('accession')} ({existing_reads} reads) over {stat.get('accession')} (0 reads)")
                continue
            elif existing_reads == 0 and current_reads > 0:
                # Replace with current (has reads)
                if "lassa" in viral_species.lower() or "lassa" in organism.lower():
                    logger.warning(f"DEBUG: Replacing Lassa reference {species_segment_dict[key].get('accession')} (0 reads) with {stat.get('accession')} ({current_reads} reads)")
                species_segment_dict[key] = stat
            elif current_is_refseq and not existing_is_refseq:
                species_segment_dict[key] = stat
            elif current_is_refseq == existing_is_refseq:
                # If both same type, prefer one with more mapped reads
                if current_reads > existing_reads:
                    species_segment_dict[key] = stat
                # If equal reads (including both 0), keep existing
    
    logger.info(f"Segment detection and species-based deduplication complete (DB: {segment_db_hits}, API: {segment_api_calls}, fallback: {segment_fallback})")
    
    # Debug: Log species that had multiple entries
    duplicates_found = {species: count for species, count in species_counts.items() if count > 1}
    if duplicates_found:
        logger.info(f"DEBUG: Found {len(duplicates_found)} viral species with multiple references before deduplication:")
        for species, count in sorted(duplicates_found.items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  - {species}: {count} reference(s)")
    
    logger.info(f"Deduplication result: {len(stats_list)} references -> {len(species_segment_dict)} unique species/segments")
    return list(species_segment_dict.values())

def safe_stem(value: str, max_len: int = 80) -> str:
    """
    Create a safe filename stem from a value.
    Removes/replaces problematic characters and limits length.
    """
    # Remove or replace problematic characters
    cleaned = re.sub(r'[^\w\-_\.]', '_', value)
    # Limit length
    if len(cleaned) > max_len:
        suffix = hashlib.md5(cleaned.encode()).hexdigest()[:10]
        cleaned = cleaned[:max_len-10]
    return f"{cleaned}_{suffix}"

def ensure_minimap2_index(database_path: Path) -> Path:
    """
    Ensure a minimap2 index (.mmi) exists for a FASTA database.
    - If `database_path` is already a .mmi, returns it.
    - If FASTA(.gz), builds `<name>.mmi` alongside the FASTA.
    Rebuilds if missing or older than the FASTA.
    """
    database_path = Path(database_path)
    if str(database_path).endswith(".mmi"):
        return database_path

    if str(database_path).endswith(".gz"):
        logger.error("minimap2 cannot build an index directly from a .gz FASTA. Please provide an uncompressed .fasta.")
        sys.exit(1)

    index_path = database_path.with_suffix(database_path.suffix + ".mmi")  # e.g. .fasta.mmi
    # Consider 0-byte index as invalid (e.g. interrupted build)
    index_exists = index_path.exists()
    index_size = index_path.stat().st_size if index_exists else 0
    needs_build = (not index_exists) or (index_size == 0) or (index_path.stat().st_mtime < database_path.stat().st_mtime)

    if needs_build:
        reason = "missing" if not index_exists else ("empty" if index_size == 0 else "outdated")
        logger.info(f"Building minimap2 index ({reason}): {index_path}")
        t0 = time.time()
        cmd = ["minimap2", "-d", str(index_path), str(database_path)]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to build minimap2 index.\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
            sys.exit(1)
        logger.info(f"Index built in {time.time() - t0:.1f}s")
    else:
        logger.info(f"Using existing minimap2 index: {index_path}")

    return index_path

def extract_accessions_from_fasta(database_fasta: Path) -> set:
    """
    Extract all accessions from a FASTA file.
    Returns a set of accessions (both versioned and unversioned).
    """
    accessions = set()
    
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")
    
    with fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = extract_accession_from_header(header)
                if acc:
                    accessions.add(acc)
                    # Also add versioned and unversioned versions
                    if '.' in acc:
                        accessions.add(acc.split('.')[0])  # Unversioned
                    else:
                        # Try to find versioned version (though we can't know the version number)
                        # Just keep the unversioned one
                        pass
    
    return accessions

def build_header_mapping(database_fasta: Path) -> dict:
    """
    Build a mapping from accession to full header description.
    This helps match SAM headers (which may be truncated) to full FASTA headers.
    """
    header_map = {}
    
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")
    
    with fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = extract_accession_from_header(header)
                if acc:
                    header_map[acc] = header
    
    return header_map

def build_organism_mapping(database_fasta: Path) -> dict:
    """
    Build a mapping from accession to organism name by parsing FASTA headers.
    Handles both pipe-delimited format (RVDB/RefSeq) and simple NCBI format.
    Returns dict: accession -> organism_name
    """
    organism_map = {}
    
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")
    
    with fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = extract_accession_from_header(header)
                if acc:
                    # Try pipe-delimited format first (RVDB/RefSeq): acc|GENBANK|accession|description|organism|VRL|date
                    # Organism is typically at index 4 (5th part)
                    parts = header.split("|")
                    if len(parts) >= 5:
                        organism = parts[4].strip()
                        if organism:  # Only add if not empty
                            organism_map[acc] = organism
                    elif len(parts) >= 4:
                        # Fallback: try parts[-2] if format is different
                        organism = parts[-2].strip()
                        if organism:
                            organism_map[acc] = organism
                    else:
                        # Simple NCBI format: >accession description
                        # For custom accessions downloaded from NCBI, the header is usually:
                        # >accession description
                        # We can't extract organism from this format reliably, so leave it empty
                        # The organism will be fetched via NCBI API in enrich_stats_with_organism
                        pass
    
    return organism_map

def build_ref_length_mapping(database_fasta: Path) -> dict:
    """
    Build a mapping from accession to reference length by parsing FASTA file.
    Returns dict: accession -> length
    """
    ref_lengths = {}
    current_acc = None
    current_length = 0
    
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")
    
    with fh:
        for line in fh:
            if line.startswith(">"):
                # Save previous sequence length
                if current_acc and current_length > 0:
                    ref_lengths[current_acc] = current_length
                
                # Start new sequence
                header = line[1:].strip()
                current_acc = extract_accession_from_header(header)
                current_length = 0
            else:
                # Count sequence length (remove whitespace)
                current_length += len(line.strip())
        
        # Save last sequence
        if current_acc and current_length > 0:
            ref_lengths[current_acc] = current_length
    
    return ref_lengths

def count_reads(fastq_file: Path) -> int:
    """Count number of reads in a FASTQ file (compressed or uncompressed)."""
    count = 0
    if str(fastq_file).endswith(".gz"):
        import gzip
        fh = gzip.open(fastq_file, "rt")
    else:
        fh = open(fastq_file, "r")
    
    with fh:
        for line in fh:
            if line.startswith("@"):
                # Check if it's a read header (not quality score)
                if len(line.split()) > 0:
                    count += 1
    
    # FASTQ has 4 lines per read, but we count headers
    # Actually, we should count every @ line that's a header
    # Let me fix this
    count = 0
    if str(fastq_file).endswith(".gz"):
        import gzip
        fh = gzip.open(fastq_file, "rt")
    else:
        fh = open(fastq_file, "r")
    
    with fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                count += 1
    
    return count

def parse_sam_per_reference_stats(sam_file: Path, min_identity: float = 0.0):
    """
    Parse SAM file and calculate per-reference statistics.
    Only counts alignments that meet the minimum identity threshold.
    Returns:
    - stats_by_ref: dict mapping SAM header -> stats dict
    - ref_lengths: dict mapping SAM header -> reference length
    - ref_headers: dict mapping SAM header -> full header string
    """
    stats_by_ref = {}
    ref_lengths = {}
    ref_headers = {}
    
    def cigar_aligned_bases(cigar: str) -> tuple:
        """Parse CIGAR string and return (reference_consumed, aligned_M_bases)."""
        num = ""
        total = 0
        aligned_m = 0
        for ch in cigar:
            if ch.isdigit():
                num += ch
                continue
            if not num:
                continue
            n = int(num)
            num = ""
            if ch in ("M", "=", "X"):
                aligned_m += n
                total += n
            elif ch in ("D", "N"):
                total += n
        return total, aligned_m
    
    with open(sam_file, "r") as f:
        for line in f:
            if line.startswith("@SQ"):
                # Parse reference length from @SQ header
                # @SQ	SN:header	LN:length
                parts = line.strip().split("\t")
                header = None
                length = None
                for part in parts:
                    if part.startswith("SN:"):
                        header = part[3:]
                    elif part.startswith("LN:"):
                        length = int(part[3:])
                if header and length:
                    ref_lengths[header] = length
                    ref_headers[header] = header
            elif line.startswith("@"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue
                
                read_name = fields[0]
                flag = int(fields[1])
                if flag & 0x4:  # unmapped
                    continue
                
                sam_header = fields[2]
                pos = int(fields[3]) - 1  # Convert to 0-based
                cigar = fields[5]
                
                if sam_header not in stats_by_ref:
                    stats_by_ref[sam_header] = {
                        "mapped_reads": 0,
                        "aligned_bases": 0,
                        "matches": 0,
                        "covered_positions": set(),  # Track covered positions for breadth calculation
                    }
                
                # Parse CIGAR
                ref_consumed, aligned_m = cigar_aligned_bases(cigar)
                
                # Get NM (number of mismatches) from optional fields
                nm = 0
                for opt in fields[11:]:
                    if opt.startswith("NM:i:"):
                        try:
                            nm = int(opt.split(":")[2])
                        except Exception:
                            nm = 0
                        break
                
                matches = aligned_m - nm
                
                # Calculate identity for THIS INDIVIDUAL READ ALIGNMENT
                # This is per-read filtering, not average-based filtering
                alignment_identity = (matches / aligned_m * 100) if aligned_m > 0 else 0.0
                
                # Only count THIS READ ALIGNMENT if it meets the minimum identity threshold
                # For RefSeq: only reads with >= 97% identity are counted
                # For RVDB: only reads with >= 80% identity are counted
                if alignment_identity >= min_identity:
                    stats_by_ref[sam_header]["mapped_reads"] += 1
                    # Use aligned_m (query length) for aligned_bases to match how matches is calculated
                    # This ensures identity = matches/aligned_bases <= 100% (matches <= aligned_m)
                    stats_by_ref[sam_header]["aligned_bases"] += aligned_m
                    stats_by_ref[sam_header]["matches"] += matches
                    # Track covered positions (for breadth calculation) - only for reads meeting threshold
                    for i in range(pos, pos + ref_consumed):
                        stats_by_ref[sam_header]["covered_positions"].add(i)
                # else: skip this individual read alignment (doesn't meet identity threshold)
                # This read alignment is NOT counted in mapped_reads, matches, aligned_bases, covered_positions, etc.
                # else: skip this individual read alignment (doesn't meet identity threshold)
                # This read alignment is NOT counted in mapped_reads, matches, etc.
    
    # Convert covered_positions set to count for each reference
    for sam_header in stats_by_ref:
        covered_set = stats_by_ref[sam_header]["covered_positions"]
        stats_by_ref[sam_header]["covered_positions"] = len(covered_set)
    
    return stats_by_ref, ref_lengths, ref_headers

def extract_accession_from_header(header: str) -> str:
    """
    Extract accession number from FASTA header.
    Handles both GENBANK and REFSEQ formats.
    """
    if not header:
        return ""
    
    parts = header.split("|")
    
    # GENBANK format: acc|GENBANK|accession|...
    if len(parts) >= 3 and parts[1].strip() == "GENBANK":
        return parts[2].strip()
    
    # REFSEQ format: acc|REFSEQ|accession|...
    if len(parts) >= 3 and parts[1].strip() == "REFSEQ":
        return parts[2].strip()
    
    # Try to find accession-like pattern (e.g., NC_123456.1, AY123456.1, etc.)
    import re
    # Pattern for common accession formats: letters, numbers, underscore, dot
    match = re.search(r'([A-Z]{1,2}_?\d+\.?\d*)', header)
    if match:
        return match.group(1)
    
    # Fallback: return first part after splitting by space or pipe
    if "|" in header:
        return header.split("|")[0].strip()
    elif " " in header:
        return header.split()[0].strip()
    
    return ""

def extract_selected_references(database_fasta: Path, selected_headers: list, out_fasta: Path) -> int:
    """
    Extract multiple reference sequences from database FASTA matching the given headers.
    Returns the number of references successfully extracted.
    """
    database_fasta = Path(database_fasta)
    out_fasta = Path(out_fasta)
    # Ensure parent directory exists
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    # Create sets for matching - use both full headers and accessions for better matching
    selected_headers_set = set(selected_headers)
    selected_accessions_set = set()
    for h in selected_headers:
        acc = extract_accession_from_header(h)
        if acc:
            selected_accessions_set.add(acc)
            # Also add base accession (without version)
            acc_base = acc.split('.')[0] if '.' in acc else acc
            selected_accessions_set.add(acc_base)
    
    found_count = 0
    write = False
    
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")
    
    with fh:
        with open(out_fasta, "w") as out:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip()
                    # Try matching on full header first
                    write = (header in selected_headers_set)
                    if not write:
                        # Fallback: try matching by accession (handles format variations)
                        header_acc = extract_accession_from_header(header)
                        if header_acc:
                            header_acc_base = header_acc.split('.')[0] if '.' in header_acc else header_acc
                            write = (header_acc in selected_accessions_set or header_acc_base in selected_accessions_set)
                    if write:
                        out.write(line)
                        found_count += 1
                else:
                    if write:
                        out.write(line)
    
    logger.info(f"Extracted {found_count}/{len(selected_headers)} selected references to {out_fasta}")
    
    # Debug: Log if any Lassa virus references were not found
    if found_count < len(selected_headers):
        missing = []
        for h in selected_headers:
            if h not in selected_headers_set or not any(extract_accession_from_header(h) == extract_accession_from_header(db_h) for db_h in [h]):  # Simplified check
                if "lassa" in h.lower():
                    missing.append(h)
        if missing:
            logger.warning(f"DEBUG: {len(missing)} Lassa virus reference(s) may not have been extracted:")
            for m in missing[:5]:  # Show first 5
                logger.warning(f"  - {m}")
    return found_count

def remap_to_selected_references(sample_fastq: Path, selected_refs_fasta: Path, output_sam: Path, min_identity: float = 80.0, threads: int = 1) -> dict:
    """
    Re-map all reads to the selected references to get accurate mapped_reads counts.
    Creates minimap2 index for the selected references database for speed.
    For low-input samples (< 100k reads), uses secondary alignments to better detect multi-mapping reads.
    Returns a dictionary mapping reference descriptions to updated stats (mapped_reads, avg_identity, etc.).
    """
    
    # Count references in the selected FASTA
    ref_count = 0
    with open(selected_refs_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                ref_count += 1
    
    # Count reads in sample to determine if we should use secondary alignments
    logger.info(f"Counting reads in sample to determine mapping strategy...")
    total_reads = count_reads(sample_fastq)
    use_secondary = total_reads < 100000  # Use secondary alignments for low-input samples (< 100k reads)
    
    if use_secondary:
        logger.info(f"Low-input sample detected ({total_reads:,} reads < 100k). Using secondary alignments in re-mapping phase (similar to initial mapping).")
    else:
        logger.info(f"Large sample detected ({total_reads:,} reads >= 100k). Using primary alignments only in re-mapping phase.")
    
    # Create minimap2 index for selected references (for speed)
    logger.info(f"Creating minimap2 index for {ref_count} selected references (this speeds up re-mapping)...")
    ref_index = ensure_minimap2_index(selected_refs_fasta)
    logger.info(f"Index created successfully. Using indexed database for re-mapping.")
    
    # Don't use -p parameter - let minimap2 report all alignments, then filter by identity during SAM parsing
    # This gives us better control and avoids missing valid hits due to strict filtering during mapping
    
    # Detect if RefSeq database (for stricter parameters to reduce false positives)
    selected_refs_path_str = Path(selected_refs_fasta).as_posix().lower()
    is_refseq = "refseq" in selected_refs_path_str
    
    # Build minimap2 command
    # For low-input samples (< 100k reads), allow secondary alignments (like initial mapping)
    # For large samples (>= 100k reads), use primary alignments only to reduce false positives
    minimap_cmd = [
        "minimap2",
        "-ax", "map-ont",
        # Don't use -p here - filter by identity during SAM parsing instead
        "-f", "0.0002",
        "-I", "8G",
        "--max-chain-skip", "25",
        "-t", str(threads),
    ]
    
    # Add alignment reporting parameters based on read count
    if use_secondary:
        # Low-input sample: For segmented viruses and multi-reference scenarios,
        # map to each reference separately to ensure each gets its reads counted correctly.
        # When mapping to multiple refs together, minimap2 prioritizes best match, causing
        # lower-identity references (e.g., Lassa 85%) to lose reads to higher-identity ones (e.g., Measles 96%).
        # By mapping separately, each reference gets primary alignments, preserving read counts.
        logger.info(f"Low-input sample with secondary alignments: mapping to each reference separately to preserve read counts...")
        
        # Map to each reference separately and combine results
        updated_stats = {}
        ref_descriptions = {}  # sam_header -> full description
        
        # Extract each reference and map separately
        with open(selected_refs_fasta, 'r') as f:
            current_ref_lines = []
            current_header = None
            
            for line in f:
                if line.startswith('>'):
                    # Process previous reference if exists
                    if current_header and current_ref_lines:
                        # Create temporary FASTA for this single reference
                        import tempfile
                        temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
                        temp_fasta.write(current_header)
                        temp_fasta.writelines(current_ref_lines)
                        temp_fasta.close()
                        temp_fasta_path = Path(temp_fasta.name)
                        
                        # Extract accession and description
                        acc = extract_accession_from_header(current_header[1:].strip())
                        desc = current_header[1:].strip()
                        ref_descriptions[acc] = desc
                        
                        # Debug: Log reference being mapped
                        if "lassa" in desc.lower() or "AY628204.1" in acc:
                            logger.info(f"DEBUG: Mapping to Lassa reference separately: {acc}")
                            # Check reference length
                            seq_len = sum(len(line.strip()) for line in current_ref_lines if not line.startswith('>'))
                            logger.info(f"DEBUG: Lassa reference sequence length: {seq_len:,} bp")
                        
                        # Create index for this single reference
                        ref_index_single = ensure_minimap2_index(temp_fasta_path)
                        
                        # Map to this single reference (primary alignments only - each ref gets its own primaries)
                        temp_sam = tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False)
                        temp_sam.close()
                        temp_sam_path = Path(temp_sam.name)
                        
                        minimap_cmd_single = [
                            "minimap2",
                            "-ax", "map-ont",
                            "-N", "1",  # Primary alignments only (but each ref gets its own primaries)
                            "--secondary=no",
                            "-f", "0.0002",
                            "-I", "8G",
                            "--max-chain-skip", "25",
                            "-t", str(threads),
                            str(ref_index_single),
                            str(sample_fastq),
                        ]
                        
                        try:
                            with open(temp_sam_path, "w") as out:
                                result = subprocess.run(minimap_cmd_single, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
                            
                            # Debug: Count alignments in SAM before parsing
                            if "lassa" in desc.lower() or "AY628204.1" in acc:
                                sam_line_count = sum(1 for line in open(temp_sam_path) if not line.startswith('@') and line.strip())
                                logger.info(f"DEBUG: Lassa SAM file has {sam_line_count:,} alignment lines (before identity filtering)")
                            
                            # Parse SAM for this single reference
                            stats_single, ref_lengths_single, ref_headers_single = parse_sam_per_reference_stats(temp_sam_path, min_identity=min_identity)
                            
                            # Debug: Log stats after parsing
                            if "lassa" in desc.lower() or "AY628204.1" in acc:
                                for sh, s in stats_single.items():
                                    logger.info(f"DEBUG: Lassa stats after SAM parsing: {s['mapped_reads']:,} reads, "
                                               f"identity={(s['matches']/s['aligned_bases']*100) if s['aligned_bases'] > 0 else 0:.2f}%")
                            
                            # Add to combined stats
                            for sam_header, s in stats_single.items():
                                aligned_i = s["aligned_bases"]
                                matches_i = s["matches"]
                                identity = (matches_i / aligned_i * 100) if aligned_i > 0 else 0
                                
                                ref_len = ref_lengths_single.get(sam_header)
                                if ref_len and ref_len > 0:
                                    coverage_depth = aligned_i / ref_len
                                    covered_pos = s.get("covered_positions", 0)
                                    coverage_breadth = covered_pos / ref_len if covered_pos > 0 else 0.0
                                else:
                                    coverage_depth = 0.0
                                    coverage_breadth = 0.0
                                
                                # Use full description as key
                                updated_stats[desc] = {
                                    "mapped_reads": s["mapped_reads"],
                                    "avg_identity": identity,
                                    "coverage_depth": coverage_depth,
                                    "coverage_breadth": coverage_breadth,
                                }
                            
                            # Clean up
                            temp_fasta_path.unlink()
                            temp_sam_path.unlink()
                            ref_index_single.unlink() if ref_index_single.exists() else None
                            
                        except subprocess.CalledProcessError as e:
                            logger.warning(f"Failed to map to reference {acc}: {e.stderr}")
                            # Clean up on error
                            if temp_fasta_path.exists():
                                temp_fasta_path.unlink()
                            if temp_sam_path.exists():
                                temp_sam_path.unlink()
                    
                    # Start new reference
                    current_header = line
                    current_ref_lines = []
                else:
                    current_ref_lines.append(line)
            
            # Process last reference
            if current_header and current_ref_lines:
                import tempfile
                temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
                temp_fasta.write(current_header)
                temp_fasta.writelines(current_ref_lines)
                temp_fasta.close()
                temp_fasta_path = Path(temp_fasta.name)
                
                acc = extract_accession_from_header(current_header[1:].strip())
                desc = current_header[1:].strip()
                ref_descriptions[acc] = desc
                
                ref_index_single = ensure_minimap2_index(temp_fasta_path)
                temp_sam = tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False)
                temp_sam.close()
                temp_sam_path = Path(temp_sam.name)
                
                minimap_cmd_single = [
                    "minimap2",
                    "-ax", "map-ont",
                    "-N", "1",
                    "--secondary=no",
                    "-f", "0.0002",
                    "-I", "8G",
                    "--max-chain-skip", "25",
                    "-t", str(threads),
                    str(ref_index_single),
                    str(sample_fastq),
                ]
                
                try:
                    with open(temp_sam_path, "w") as out:
                        subprocess.run(minimap_cmd_single, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
                    
                    stats_single, ref_lengths_single, ref_headers_single = parse_sam_per_reference_stats(temp_sam_path, min_identity=min_identity)
                    
                    for sam_header, s in stats_single.items():
                        aligned_i = s["aligned_bases"]
                        matches_i = s["matches"]
                        identity = (matches_i / aligned_i * 100) if aligned_i > 0 else 0
                        
                        ref_len = ref_lengths_single.get(sam_header)
                        if ref_len and ref_len > 0:
                            coverage_depth = aligned_i / ref_len
                            covered_pos = s.get("covered_positions", 0)
                            coverage_breadth = covered_pos / ref_len if covered_pos > 0 else 0.0
                        else:
                            coverage_depth = 0.0
                            coverage_breadth = 0.0
                        
                        updated_stats[desc] = {
                            "mapped_reads": s["mapped_reads"],
                            "avg_identity": identity,
                            "coverage_depth": coverage_depth,
                            "coverage_breadth": coverage_breadth,
                        }
                    
                    temp_fasta_path.unlink()
                    temp_sam_path.unlink()
                    ref_index_single.unlink() if ref_index_single.exists() else None
                    
                except subprocess.CalledProcessError as e:
                    logger.warning(f"Failed to map to reference {acc}: {e.stderr}")
                    if temp_fasta_path.exists():
                        temp_fasta_path.unlink()
                    if temp_sam_path.exists():
                        temp_sam_path.unlink()
        
        # Build header mapping for description matching
        selected_header_map = build_header_mapping(selected_refs_fasta)
        
        # Create combined SAM file from individual mappings (needed for create_per_reference_outputs)
        # We'll collect all individual SAM files and combine them
        individual_sam_files = []
        with open(selected_refs_fasta, 'r') as f:
            current_ref_lines = []
            current_header = None
            temp_sam_files = []
            
            for line in f:
                if line.startswith('>'):
                    if current_header and current_ref_lines:
                        # Find corresponding temp SAM file (we'll need to track these)
                        # For now, we'll recreate the mapping to get the SAM files
                        # This is inefficient but ensures correctness
                        pass
                    current_header = line
                    current_ref_lines = []
                else:
                    current_ref_lines.append(line)
        
        # Actually, a simpler approach: re-map to all references together to create the combined SAM
        # But use the separate mapping stats we already computed
        # This is a bit redundant but ensures we have the combined SAM file
        logger.info("Creating combined SAM file for per-reference outputs...")
        minimap_cmd_combined = [
            "minimap2",
            "-ax", "map-ont",
            "-N", "1",  # Primary alignments only
            "--secondary=no",
            "-f", "0.0002",
            "-I", "8G",
            "--max-chain-skip", "25",
            "-t", str(threads),
            str(ref_index),
            str(sample_fastq),
        ]
        try:
            with open(output_sam, "w") as out:
                subprocess.run(minimap_cmd_combined, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            logger.warning(f"Failed to create combined SAM file: {e.stderr}")
            # Continue anyway - per-reference outputs might fail but stats are correct
        
        logger.info(f"Re-mapping complete (mapped to each reference separately for accurate stats). Updated stats for {len(updated_stats)} references.")
        
        # Debug: Log Lassa virus in updated_stats
        lassa_in_updated = [k for k in updated_stats.keys() if "lassa" in k.lower() or "AY628204.1" in k]
        if lassa_in_updated:
            logger.info(f"DEBUG: Found Lassa virus in updated_stats: {len(lassa_in_updated)} entry(ies)")
            for lassa_key in lassa_in_updated:
                lassa_data = updated_stats[lassa_key]
                logger.info(f"  - {lassa_key[:80]}: {lassa_data['mapped_reads']:,} reads, "
                           f"identity={lassa_data['avg_identity']:.2f}%, "
                           f"breadth={lassa_data['coverage_breadth']:.4f}")
        else:
            logger.warning(f"DEBUG: Lassa virus NOT found in updated_stats after remapping!")
            logger.warning(f"DEBUG: Available keys in updated_stats: {list(updated_stats.keys())[:5]}")
        
        return updated_stats
    else:
        # Large sample: use primary alignments only to reduce false positives and speed up
        minimap_cmd.extend([
            "-N", "1",  # Only report 1 alignment per read (primary alignment only)
            "--secondary=no",  # No secondary alignments
        ])
        logger.info(f"Re-mapping with primary alignments only (for large sample with {total_reads:,} reads)...")
    
    minimap_cmd.extend([
        str(ref_index),  # Use index instead of FASTA for speed
        str(sample_fastq),
    ])
    
    logger.info(f"Re-mapping all reads to selected references using indexed database (will filter by identity >= {min_identity}% during SAM parsing, {'with secondary alignments' if use_secondary else 'primary alignments only'})...")
    try:
        with open(output_sam, "w") as out:
            subprocess.run(minimap_cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"minimap2 re-mapping failed:\n{e.stderr}")
        return {}
    
    # Build header mapping from selected references FASTA to get full descriptions
    selected_header_map = build_header_mapping(selected_refs_fasta)
    
    # Parse SAM to get accurate counts per reference (filter by identity during parsing)
    stats_by_ref, ref_lengths, ref_headers = parse_sam_per_reference_stats(output_sam, min_identity=min_identity)
    
    # Convert to dictionary with descriptions as keys (matching the format from curated_stats_dedup)
    updated_stats = {}
    for sam_header, s in stats_by_ref.items():
        aligned_i = s["aligned_bases"]
        matches_i = s["matches"]
        identity = (matches_i / aligned_i * 100) if aligned_i > 0 else 0
        
        # Get accession and full description from selected references FASTA
        acc = extract_accession_from_header(sam_header)
        # Try to get full description from header map (using accession)
        desc = selected_header_map.get(acc, ref_headers.get(sam_header, sam_header))
        
        # If still not found, try to match by SAM header directly
        if desc == sam_header:
            # Try to find in selected_header_map by matching any header that contains this accession
            for full_header, full_desc in selected_header_map.items():
                if acc in full_header or acc in full_desc:
                    desc = full_desc
                    break
        
        ref_len = ref_lengths.get(sam_header)
        if ref_len and ref_len > 0:
            coverage_depth = aligned_i / ref_len
            covered_pos = s.get("covered_positions", 0)
            coverage_breadth = covered_pos / ref_len if covered_pos > 0 else 0.0
        else:
            coverage_depth = 0.0
            coverage_breadth = 0.0
        
        updated_stats[desc] = {
            "mapped_reads": s["mapped_reads"],
            "avg_identity": identity,
            "coverage_depth": coverage_depth,
            "coverage_breadth": coverage_breadth,
        }
    
    logger.info(f"Re-mapping complete. Updated stats for {len(updated_stats)} references.")
    
    # Debug: Log Lassa virus in updated_stats
    lassa_in_updated = [k for k in updated_stats.keys() if "lassa" in k.lower() or "AY628204.1" in k]
    if lassa_in_updated:
        logger.info(f"DEBUG: Found Lassa virus in updated_stats: {len(lassa_in_updated)} entry(ies)")
        for lassa_key in lassa_in_updated:
            lassa_data = updated_stats[lassa_key]
            logger.info(f"  - {lassa_key[:80]}: {lassa_data['mapped_reads']:,} reads, "
                       f"identity={lassa_data['avg_identity']:.2f}%, "
                       f"breadth={lassa_data['coverage_breadth']:.4f}")
    else:
        logger.warning(f"DEBUG: Lassa virus NOT found in updated_stats after remapping!")
        logger.warning(f"DEBUG: Available keys in updated_stats: {list(updated_stats.keys())[:5]}")
    
    logger.debug(f"Updated stats keys: {list(updated_stats.keys())[:3]}...")  # Debug: show first few keys
    return updated_stats

def create_per_reference_outputs(sample_name: str, curated_descriptions: list, selected_refs_fasta: Path, remap_sam: Path, sample_fastq: Path, sample_dir: Path):
    """
    Create per-reference folders with SAM, BAM, reference FASTA, and mapped reads FASTQ.
    For each reference in curated_descriptions.json, creates a folder named after its accession.
    """
    if not curated_descriptions:
        logger.info("No curated references to create per-reference outputs")
        return
    
    if not selected_refs_fasta.exists():
        logger.warning(f"Selected references FASTA not found: {selected_refs_fasta}. Skipping per-reference outputs.")
        return
    
    if not remap_sam.exists():
        logger.warning(f"Remapped SAM file not found: {remap_sam}. Skipping per-reference outputs.")
        return
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Creating per-reference output folders for {len(curated_descriptions)} reference(s)...")
    
    # Build a mapping of description -> accession for quick lookup
    desc_to_acc = {stat.get("description", ""): stat.get("accession", "") for stat in curated_descriptions}
    
    # Read remapped SAM file and group by reference
    logger.info("Reading remapped SAM file...")
    sam_by_reference = {}  # accession -> list of SAM lines
    read_ids_by_reference = {}  # accession -> set of read IDs
    
    with open(remap_sam, 'r') as f:
        for line in f:
            if line.startswith('@'):
                # Header line - store for later use
                continue
            if not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            
            # Get reference name from SAM line (3rd field)
            ref_name = parts[2]
            if ref_name == '*':
                continue
            
            # Find accession for this reference description
            accession = None
            for desc, acc in desc_to_acc.items():
                if ref_name in desc or desc in ref_name:
                    accession = acc
                    break
            
            if not accession:
                # Try to extract accession from ref_name directly
                accession = extract_accession_from_header(ref_name)
            
            if accession:
                if accession not in sam_by_reference:
                    sam_by_reference[accession] = []
                    read_ids_by_reference[accession] = set()
                sam_by_reference[accession].append(line)
                # Get read ID (1st field) - handle /1, /2 suffixes for paired-end reads
                read_id = parts[0].split('/')[0]  # Remove /1, /2 suffix if present
                read_ids_by_reference[accession].add(read_id)
    
    logger.info(f"Found alignments for {len(sam_by_reference)} reference(s)")
    
    # Process each curated reference
    for stat in curated_descriptions:
        accession = stat.get("accession", "")
        description = stat.get("description", "")
        
        if not accession:
            logger.warning(f"Skipping reference without accession: {description[:50]}...")
            continue
        
        # Create folder for this accession
        acc_dir = sample_dir / accession
        acc_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Processing {accession}...")
        
        # 1. Extract reference FASTA
        ref_fasta = acc_dir / f"{accession}.fasta"
        if not extract_fasta_record(selected_refs_fasta, description, ref_fasta):
            logger.warning(f"  Could not extract reference FASTA for {accession}")
        
        # 2. Filter SAM file for this reference
        if accession in sam_by_reference:
            sam_lines = sam_by_reference[accession]
            ref_sam = acc_dir / f"{accession}.sam"
            
            # Write SAM header first (from original remapped.sam)
            with open(remap_sam, 'r') as f_in:
                with open(ref_sam, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('@'):
                            f_out.write(line)
                        else:
                            break
            
            # Write alignments for this reference
            with open(ref_sam, 'a') as f:
                for line in sam_lines:
                    f.write(line)
            
            logger.info(f"  Created {ref_sam.name} ({len(sam_lines)} alignments)")
            
            # 3. Convert SAM to coordinate-sorted BAM (+ index) for easy visualization
            # Note: `samtools index` requires coordinate-sorted BAM.
            # The sorted BAM will replace any existing BAM file with the same name.
            ref_bam = acc_dir / f"{accession}.bam"
            ref_bai = acc_dir / f"{accession}.bam.bai"
            try:
                # Check if samtools is available
                subprocess.run(['samtools', '--version'], check=True, capture_output=True)

                # Remove existing BAM and BAI files if they exist (will be replaced with sorted version)
                if ref_bam.exists():
                    ref_bam.unlink()
                if ref_bai.exists():
                    ref_bai.unlink()

                # Create coordinate-sorted BAM (samtools sort writes directly to output file, replacing if exists):
                # samtools view -bS ref.sam | samtools sort -o ref.bam -
                # The -o flag writes directly to ref.bam (replaces existing file)
                view_cmd = ['samtools', 'view', '-bS', str(ref_sam)]
                sort_cmd = ['samtools', 'sort', '-o', str(ref_bam), '-']

                view_p = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                try:
                    sort_p = subprocess.Popen(sort_cmd, stdin=view_p.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                finally:
                    # Allow view_p to receive SIGPIPE if sort exits early
                    if view_p.stdout:
                        view_p.stdout.close()

                _, view_err = view_p.communicate()
                _, sort_err = sort_p.communicate()

                if view_p.returncode != 0:
                    raise subprocess.CalledProcessError(view_p.returncode, view_cmd, output=None, stderr=view_err)
                if sort_p.returncode != 0:
                    raise subprocess.CalledProcessError(sort_p.returncode, sort_cmd, output=None, stderr=sort_err)

                # Verify BAM file was created and is sorted
                if not ref_bam.exists():
                    raise FileNotFoundError(f"Sorted BAM file was not created: {ref_bam}")

                # Index BAM (creates .bai next to the BAM, replaces existing index)
                subprocess.run(['samtools', 'index', str(ref_bam)], check=True, capture_output=True)

                if ref_bai.exists():
                    logger.info(f"  Created {ref_bam.name} and {ref_bai.name}")
                else:
                    # Some samtools versions might name differently or write elsewhere; still log success.
                    logger.info(f"  Created {ref_bam.name} and indexed it")
                # Remove SAM file since we have BAM
                ref_sam.unlink()
                logger.info(f"  Removed {ref_sam.name} (BAM file available)")
            except (subprocess.CalledProcessError, FileNotFoundError):
                logger.warning(f"  Could not create BAM file (samtools not available or failed)")
                # Still remove SAM file to save space (we have the data in FASTQ and can recreate if needed)
                if ref_sam.exists():
                    ref_sam.unlink()
                    logger.info(f"  Removed {ref_sam.name} (to save space, BAM creation failed but data available in FASTQ)")
            
            # 4. Extract reads from original FASTQ that mapped ONLY to this reference
            if accession in read_ids_by_reference:
                # Get read IDs that mapped to THIS specific reference
                read_ids_for_this_ref = read_ids_by_reference[accession]
                
                # Count how many reads mapped to this reference
                logger.info(f"  Extracting {len(read_ids_for_this_ref)} read(s) that mapped to {accession}...")
                
                ref_fastq = acc_dir / f"{accession}_mapped_reads.fastq"
                
                # Read original FASTQ and extract ONLY reads that mapped to this reference
                sample_fastq_path = Path(sample_fastq)
                is_gzipped = sample_fastq_path.suffix == '.gz'
                
                read_count = 0
                with open(ref_fastq, 'w') as f_out:
                    if is_gzipped:
                        import gzip
                        f_in = gzip.open(sample_fastq_path, 'rt')
                    else:
                        f_in = open(sample_fastq_path, 'r')
                    
                    with f_in:
                        current_read_id = None
                        current_read_lines = []
                        for line in f_in:
                            if line.startswith('@'):
                                # Save previous read if it mapped to THIS reference
                                if current_read_id and current_read_id in read_ids_for_this_ref:
                                    f_out.writelines(current_read_lines)
                                    read_count += 1
                                # Start new read
                                # Handle both @read_id and @read_id/1 format
                                read_id_line = line.strip()[1:]  # Remove @
                                current_read_id = read_id_line.split()[0].split('/')[0]  # Get read ID (handle /1, /2 suffixes)
                                current_read_lines = [line]
                            else:
                                current_read_lines.append(line)
                        
                        # Don't forget the last read
                        if current_read_id and current_read_id in read_ids_for_this_ref:
                            f_out.writelines(current_read_lines)
                            read_count += 1
                
                logger.info(f"  Created {ref_fastq.name} ({read_count} reads mapped to {accession} only)")
        else:
            logger.warning(f"  No alignments found for {accession} in remapped SAM")
    
    logger.info(f"Per-reference outputs created in {sample_dir}")
    
    # Clean up: Remove main remapped SAM file (per-reference SAM/BAM files are in their folders)
    if remap_sam.exists():
        remap_sam.unlink()
        logger.info(f"Removed main remapped SAM file: {remap_sam.name}")
    
    # Clean up: Remove initial mapping SAM file (sample_name.sam) - no longer needed
    initial_sam = sample_dir / f"{sample_name}.sam"
    if initial_sam.exists():
        initial_sam.unlink()
        logger.info(f"Removed initial mapping SAM file: {initial_sam.name}")

def find_best_reference_with_index(sample_fastq: Path, database_fasta: Path, output_dir: Path, sample_name: str, minimap_p: float = None, min_identity: float = 80.0, min_mapped_reads: int = 100, coverage_depth_threshold: float = 1.0, coverage_breadth_threshold: float = 0.1, threads: int = 1):
    """
    Find best reference using minimap2 with indexed database (single-pass mapping).
    Returns best_ref, best_stats, all_stats, filtered_stats, curated_stats.
    
    If minimap_p is None, it will be calculated from min_identity (min_identity / 100).
    This allows minimap2 to filter low-quality alignments during mapping (like metamaps --pi).
    
    Note: output_dir is already the sample directory (or database-specific subdirectory if multiple databases).
    """
    # output_dir is already the final directory (sample_dir or sample_dir/database)
    sample_dir = output_dir
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Build header mapping for matching SAM headers to full descriptions
    logger.info("Building header mapping from database...")
    header_map = build_header_mapping(Path(database_fasta))
    
    # Build organism mapping from database FASTA (accession -> organism name)
    logger.info("Building organism mapping from database...")
    organism_map = build_organism_mapping(Path(database_fasta))
    
    # Build reference length mapping from database FASTA (for cases where SAM headers don't have @SQ LN)
    logger.info("Building reference length mapping from database...")
    db_ref_lengths = build_ref_length_mapping(Path(database_fasta))
    
    # Ensure minimap2 index exists
    logger.info("Ensuring minimap2 index exists for database...")
    db_index = ensure_minimap2_index(Path(database_fasta))
    
    aln_sam = sample_dir / f"{sample_name}.sam"

    # Don't use -p parameter - let minimap2 report all alignments, then filter by identity during SAM parsing
    # The -p parameter filters during mapping which is too strict and can miss valid hits
    # Instead, we filter by identity during SAM parsing which gives us full control
    if minimap_p is not None:
        logger.info(f"Using minimap2 -p {minimap_p:.2f} (user-specified)")
    else:
        logger.info(f"Using minimap2 without -p filter (will filter by identity >= {min_identity}% during SAM parsing)")

    # Build minimap2 command
    # Use primary alignments only (-N 1) to reduce false positives
    minimap_cmd = [
        "minimap2",
        "-ax", "map-ont",
        "-N", "1",  # Only report 1 alignment per read (primary alignment only)
        "--secondary=no",  # No secondary alignments
        "-f", "0.0002",  # Filter top 0.02% frequent minimizers (faster, minimal sensitivity loss)
        "-I", "8G",  # Batch size for better memory efficiency
        "--max-chain-skip", "25",  # Limit chain extension for speed
        "-t", str(threads),  # Use specified number of threads
    ]
    
    # Only add -p if user explicitly specified it
    if minimap_p is not None:
        minimap_cmd.append("-p")
        minimap_cmd.append(str(minimap_p))
    
    minimap_cmd.extend([
        str(db_index),  # Use index for speed
        str(sample_fastq),
    ])
    
    if minimap_p is not None:
        logger.info(f"Running minimap2 (single pass against indexed DB, {threads} threads, -p {minimap_p:.2f} for alignment filtering)...")
    else:
        logger.info(f"Running minimap2 (single pass against indexed DB, {threads} threads, will filter by identity >= {min_identity}% during SAM parsing)...")
    try:
        with open(aln_sam, "w") as out:
            subprocess.run(minimap_cmd, check=True, stdout=out, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"minimap2 failed:\n{e.stderr}")
        sys.exit(1)

    stats_by_ref, ref_lengths, ref_headers = parse_sam_per_reference_stats(aln_sam, min_identity=min_identity)
    if not stats_by_ref:
        logger.warning("No reads mapped to any reference.")
        # Clean up SAM file before returning
        if aln_sam.exists():
            aln_sam.unlink()
            logger.info(f"Removed SAM file (no reads mapped): {aln_sam.name}")
        # Return empty/default values for all 5 return values
        return None, None, [], [], []

    # Calculate identity for all references and filter by minimum identity threshold
    all_stats = []
    
    # Debug: Track AY628205.1 specifically
    debug_accession = "AY628205.1"
    debug_found = False
    
    for sam_header, s in stats_by_ref.items():
        aligned_i = s["aligned_bases"]
        matches_i = s["matches"]
        identity = (matches_i / aligned_i * 100) if aligned_i > 0 else 0
        
        acc = extract_accession_from_header(sam_header)
        # Look up full header from FASTA using accession (SAM header might be truncated)
        desc = header_map.get(acc, ref_headers.get(sam_header, sam_header))
        
        # Debug: Log AY628205.1 if found
        if acc == debug_accession or debug_accession in sam_header:
            debug_found = True
            ref_len = ref_lengths.get(sam_header) or db_ref_lengths.get(acc, 0)
            coverage_depth = aligned_i / ref_len if ref_len > 0 else 0.0
            covered_pos = s.get("covered_positions", 0)
            coverage_breadth = covered_pos / ref_len if ref_len > 0 else 0.0
            logger.info(f"DEBUG: Found {debug_accession} in SAM parsing: {s['mapped_reads']} reads, "
                       f"identity={identity:.2f}%, depth={coverage_depth:.2f}, breadth={coverage_breadth:.4f}")
        
        # Calculate coverage depth (average depth) and breadth (horizontal coverage)
        ref_len = ref_lengths.get(sam_header)
        # If ref_length missing from SAM headers, try to get it from database FASTA
        if not ref_len or ref_len <= 0:
            ref_len = db_ref_lengths.get(acc, 0)
        
        if ref_len and ref_len > 0:
            coverage_depth = aligned_i / ref_len  # Average depth (can be > 1.0, e.g., 10X = 10.0)
            covered_pos = s.get("covered_positions", 0)
            coverage_breadth = covered_pos / ref_len if covered_pos > 0 else 0.0  # Horizontal coverage (0.0 to 1.0)
        else:
            # If reference length is still missing, log warning and set coverage to 0 (will be filtered out)
            logger.warning(f"Reference length missing for {acc} ({sam_header}) - coverage set to 0")
            coverage_depth = 0.0
            coverage_breadth = 0.0
        
        stat_entry = {
            "accession": acc,
            "description": desc,
            "mapped_reads": s["mapped_reads"],
            "coverage_depth": coverage_depth,  # Average depth (vertical coverage, can be > 1.0)
            "coverage_breadth": coverage_breadth,  # Horizontal coverage (fraction of positions covered, 0.0 to 1.0)
            "avg_identity": identity,
            "aligned_bases": aligned_i,  # Store for aggregation
            "matches": matches_i,  # Store for aggregation
            "ref_length": ref_len,  # Store for aggregation
            "covered_positions": covered_pos,  # Store for aggregation
        }
        all_stats.append(stat_entry)
    
    if not debug_found:
        logger.debug(f"DEBUG: {debug_accession} not found in SAM parsing results")
        
    # Aggregate by organism/species BEFORE filtering to combine reads that map to multiple references of same species
    # IMPORTANT: For segmented viruses (e.g., Lassa, Influenza), keep segments separate by including segment in key
    # This helps when reads are spread across many references (e.g., many MPOX references)
    # but prevents different segments from being aggregated together
    logger.info("Aggregating reads by organism/species before filtering (combining reads across multiple references of same species, but keeping segments separate)...")
    
    # Get segment database path for segment extraction
    database_dir = Path(database_fasta).parent
    segment_database_path = get_segment_database_path(database_dir)
    
    aggregated_stats = {}
    for stat in all_stats:
        # Use accession to look up organism from database mapping (more reliable than parsing description)
        accession = stat.get("accession", "")
        description = stat.get("description", "")
        organism_key = None
        
        if accession and accession in organism_map:
            # Look up organism from database mapping
            organism = organism_map[accession]
            organism_key = organism.strip().lower().replace(" ", "_")[:50]  # Limit length
        else:
            # Fallback: try to extract from description if not found in mapping
            desc = description.lower()
            if "monkeypox" in desc or "mpox" in desc:
                organism_key = "monkeypox_virus"
            elif "orthopoxvirus" in desc:
                organism_key = "orthopoxvirus"
            else:
                # Try parsing description as last resort
                parts = description.split("|")
                if len(parts) >= 5:
                    organism_key = parts[4].strip().lower().replace(" ", "_")[:50]
                elif len(parts) >= 4:
                    organism_key = parts[-2].strip().lower().replace(" ", "_")[:50]
                else:
                    # Final fallback: use accession as key (no aggregation)
                    organism_key = accession if accession else "unknown"
        
        # Extract segment information to keep segmented viruses separate
        segment = extract_segment_from_description(description, accession=accession, database_path=segment_database_path)
        
        # Include segment in aggregation key if present (prevents different segments from being aggregated)
        if segment:
            organism_key = f"{organism_key}__segment_{segment}"
        
        if organism_key not in aggregated_stats:
            aggregated_stats[organism_key] = {
                "accession": stat["accession"],  # Keep best accession (prefer RefSeq)
                "description": stat["description"],
                "mapped_reads": 0,
                "aligned_bases": 0,
                "matches": 0,
                "ref_length": 0,
                "covered_positions": set(),
                "references": [stat],  # Track all references for this organism
            }
        
        # Aggregate: sum reads, bases, matches
        aggregated_stats[organism_key]["mapped_reads"] += stat["mapped_reads"]
        aggregated_stats[organism_key]["aligned_bases"] += stat.get("aligned_bases", 0)
        aggregated_stats[organism_key]["matches"] += stat.get("matches", 0)
        # Use max ref_length (represents the reference genome length for coverage calculation)
        aggregated_stats[organism_key]["ref_length"] = max(
            aggregated_stats[organism_key]["ref_length"], 
            stat.get("ref_length", 0)
        )
        # For covered_positions, we need to properly combine them when aggregating
        # Since covered_positions is stored as integer (count) after parsing, we can't union exact positions
        # Instead, we use max() as approximation, but this can underestimate breadth when multiple references overlap
        # Note: This is a limitation - ideally we'd track positions as sets and union them
        current_covered = aggregated_stats[organism_key]["covered_positions"]
        stat_covered = stat.get("covered_positions", 0)
        
        # If both are sets (shouldn't happen after parsing, but handle it), union them
        if isinstance(current_covered, set) and isinstance(stat_covered, set):
            aggregated_stats[organism_key]["covered_positions"] = current_covered.union(stat_covered)
        elif isinstance(current_covered, set):
            # Current is set, stat is int - convert set to int and use max
            aggregated_stats[organism_key]["covered_positions"] = max(len(current_covered), stat_covered)
        elif isinstance(stat_covered, set):
            # Stat is set, current is int - convert set to int and use max
            aggregated_stats[organism_key]["covered_positions"] = max(current_covered, len(stat_covered))
        else:
            # Both are integers - use max (approximation, underestimates when references overlap)
            # TODO: This could be improved by tracking positions as sets during aggregation
            aggregated_stats[organism_key]["covered_positions"] = max(
                current_covered if current_covered else 0,
                stat_covered if stat_covered else 0
            )
        # Prefer RefSeq accession
        if "refseq" in stat.get("description", "").lower() and "refseq" not in aggregated_stats[organism_key]["description"].lower():
            aggregated_stats[organism_key]["accession"] = stat["accession"]
            aggregated_stats[organism_key]["description"] = stat["description"]
    
    # Convert aggregated stats back to list format with recalculated metrics
    aggregated_all_stats = []
    for organism_key, agg in aggregated_stats.items():
        # Get ref_length from aggregation, but try database FASTA if missing
        ref_len = agg.get("ref_length", 0)
        if ref_len <= 0:
            # Try to get ref_length from database FASTA using accession
            acc = agg.get("accession", "")
            if acc and acc in db_ref_lengths:
                ref_len = db_ref_lengths[acc]
            if ref_len <= 0:
                logger.warning(f"Reference length missing for aggregated entry {acc} - using max covered position as estimate")
                # As last resort, estimate ref_length from max covered position
                covered_pos_count = len(agg["covered_positions"]) if isinstance(agg["covered_positions"], set) else agg.get("covered_positions", 0)
                if covered_pos_count > 0:
                    # Estimate: assume breadth should be at least what we see, so ref_len >= covered_positions
                    ref_len = max(covered_pos_count, 1000)  # Minimum 1000bp estimate
                else:
                    ref_len = 1  # Only use 1 as absolute last resort
        
        aligned_bases = agg["aligned_bases"]
        matches = agg["matches"]
        
        # Recalculate aggregated metrics
        agg_identity = (matches / aligned_bases * 100) if aligned_bases > 0 else 0.0
        agg_coverage_depth = aligned_bases / ref_len if ref_len > 0 else 0.0
        covered_pos_count = len(agg["covered_positions"]) if isinstance(agg["covered_positions"], set) else agg.get("covered_positions", 0)
        agg_coverage_breadth = covered_pos_count / ref_len if ref_len > 0 else 0.0
        
        aggregated_all_stats.append({
            "accession": agg["accession"],
            "description": agg["description"],
            "mapped_reads": agg["mapped_reads"],
            "coverage_depth": agg_coverage_depth,
            "coverage_breadth": agg_coverage_breadth,
            "avg_identity": agg_identity,
        })
    
    logger.info(f"Aggregated {len(all_stats)} references -> {len(aggregated_all_stats)} organisms/species")
    # Replace all_stats with aggregated version
    all_stats = aggregated_all_stats
    
    # Apply both filters: identity >= min_identity AND mapped_reads >= min_mapped_reads
    filtered_all_stats = [s for s in all_stats if s["avg_identity"] >= min_identity and s["mapped_reads"] >= min_mapped_reads]
    
    if not filtered_all_stats:
        logger.warning(f"No references meet the filtering criteria (identity >= {min_identity:.1f}% AND mapped_reads >= {min_mapped_reads})")
        if all_stats:
            best_identity = max(s['avg_identity'] for s in all_stats)
            best_reads = max(s['mapped_reads'] for s in all_stats)
            logger.info(f"Best available: identity={best_identity:.2f}%, mapped_reads={best_reads}")
        # Return the best from all stats even if below threshold (with warning)
        filtered_all_stats = all_stats
    
    # Create curated JSON with hardcoded thresholds: 50% identity and 100 mapped reads
    curated_stats = [s for s in all_stats if s["avg_identity"] >= 50.0 and s["mapped_reads"] >= 100]
    
    # Further curate: only keep entries with avg_identity >= min_identity AND coverage_depth >= coverage_depth_threshold AND coverage_breadth >= coverage_breadth_threshold
    # Filter out entries with missing/invalid coverage or below thresholds
    initial_count = len(curated_stats)
    curated_stats = [
        s for s in curated_stats 
        if s["avg_identity"] >= min_identity 
        and s.get("coverage_depth", 0.0) >= coverage_depth_threshold 
        and s.get("coverage_depth", 0.0) > 0.0
        and s.get("coverage_breadth", 0.0) >= coverage_breadth_threshold
    ]
    
    
    filtered_out = initial_count - len(curated_stats)
    if filtered_out > 0:
        logger.info(f"After strict curation (identity >= {min_identity}%, coverage_depth >= {coverage_depth_threshold}, coverage_breadth >= {coverage_breadth_threshold}): {len(curated_stats)} references remain (filtered out {filtered_out})")
    else:
        logger.info(f"After strict curation (identity >= {min_identity}%, coverage_depth >= {coverage_depth_threshold}, coverage_breadth >= {coverage_breadth_threshold}): {len(curated_stats)} references remain")
    
    if not filtered_all_stats:
        logger.warning(f"No references meet the filtering criteria (identity >= {min_identity:.1f}% AND mapped_reads >= {min_mapped_reads})")
        if all_stats:
            best_identity = max(s['avg_identity'] for s in all_stats)
            best_reads = max(s['mapped_reads'] for s in all_stats)
            logger.info(f"Best available: identity={best_identity:.2f}%, mapped_reads={best_reads}")
        # Return the best from all stats even if below threshold (with warning)
        filtered_all_stats = all_stats
    else:
        logger.info(f"Filtered to {len(filtered_all_stats)} references (identity >= {min_identity:.1f}% AND mapped_reads >= {min_mapped_reads}) out of {len(all_stats)} total")
    
    # Log curated stats info (handle empty case gracefully)
    if curated_stats:
        logger.info(f"Curated JSON: {len(curated_stats)} references (identity >= {min_identity}% AND coverage_depth >= {coverage_depth_threshold} AND coverage_breadth >= {coverage_breadth_threshold} from initial set of identity >= 50.0% AND mapped_reads >= 100) out of {len(all_stats)} total")
    else:
        logger.info(f"Curated JSON: 0 references (identity >= {min_identity}% AND coverage_depth >= {coverage_depth_threshold} AND coverage_breadth >= {coverage_breadth_threshold} from initial set of identity >= 50.0% AND mapped_reads >= 100) out of {len(all_stats)} total")
        logger.warning("No references meet curated criteria - results may be limited")
    
    # Choose best by mapped read count from filtered results (now working with aggregated stats)
    if filtered_all_stats:
        best_stat = max(filtered_all_stats, key=lambda s: s["mapped_reads"])
    elif all_stats:
        best_stat = max(all_stats, key=lambda s: s["mapped_reads"])
    else:
        # Fallback: return None if no stats available
        return None, None, [], [], []
    
    best_acc = best_stat["accession"]
    best_desc = best_stat["description"]
    
    # Calculate best stats from aggregated entry
    best_stats = {
        "mapped_reads": best_stat["mapped_reads"],
        "coverage": best_stat["coverage_depth"],  # Already calculated in aggregated stats
        "avg_identity": best_stat["avg_identity"],
        "coverage_depth": best_stat["coverage_depth"],
        "coverage_breadth": best_stat["coverage_breadth"],
    }

    best_ref = {"accession": best_acc, "description": best_desc, "sequence": None}
    # Return all stats, filtered stats, and curated stats (50% identity, 100 mapped reads)
    return best_ref, best_stats, all_stats, filtered_all_stats, curated_stats

def get_cpu_count():
    """Get the number of available CPU cores."""
    try:
        cpu_count = multiprocessing.cpu_count()
        # Use all CPUs for minimap2 (it handles threading efficiently)
        threads = cpu_count
        logger.info(f"Using {threads} CPU cores (all available)")
        return threads
    except Exception as e:
        logger.warning(f"Could not determine CPU count: {e}, using 1 thread")
        return 1

def extract_fasta_record(database_fasta: Path, target_header: str, out_fasta: Path) -> bool:
    """
    Stream-scan a FASTA and write the record matching `target_header` to `out_fasta`.
    Matches on the full header (everything after '>' up to newline, stripped).
    Returns True if found.
    """
    database_fasta = Path(database_fasta)
    if str(database_fasta).endswith(".gz"):
        import gzip
        fh = gzip.open(database_fasta, "rt")
    else:
        fh = open(database_fasta, "r")

    found = False
    write = False
    with fh:
        with open(out_fasta, "w") as out:
            for line in fh:
                if line.startswith(">"):
                    header = line[1:].strip()
                    # Match on full header
                    write = (header == target_header)
                    if write:
                        out.write(line)
                        found = True
                else:
                    if write:
                        out.write(line)
    return found

def process_sample(sample_name, sample_fastq, database_fasta, output_dir, min_identity=80.0, min_mapped_reads=100, coverage_depth_threshold=1.0, coverage_breadth_threshold=0.1, threads=1):
    """
    Process a single sample: map ALL reads to database and find best reference.
    Simple workflow: no rarefaction, no unmapped extraction, no minority detection.
    """
    logger.info(f"\n{'='*60}")
    logger.info(f"Processing sample: {sample_name}")
    logger.info(f"{'='*60}")
    logger.info(f"Mapping ALL reads to database (no rarefaction)...")
    
    # Use the provided output_dir (which may already be database-specific if multiple databases)
    sample_dir = output_dir
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Find best reference using ALL reads (single-pass mapping using DB index)
    # minimap_p=None means it will be calculated from min_identity (like metamaps --pi)
    # Pass sample_dir directly (it's already the final directory, no need to add sample_name again)
    best_ref, best_stats, all_stats, filtered_stats, curated_stats = find_best_reference_with_index(
        sample_fastq,  # Use ALL reads, no rarefaction
        database_fasta,
        sample_dir,  # Use sample_dir directly (already database-specific if needed)
        sample_name,
        minimap_p=None,  # Will be calculated from min_identity (e.g., 95% identity -> -p 0.95)
        min_identity=min_identity,
        min_mapped_reads=min_mapped_reads,
        coverage_depth_threshold=coverage_depth_threshold,
        coverage_breadth_threshold=coverage_breadth_threshold,
        threads=threads
    )
    
    if not best_ref:
        logger.warning("No reads mapped to any reference for this sample. Skipping reference selection.")
        # Still save empty results so the user knows the sample was processed
        save_results(
            sample_name, None, None, all_stats, filtered_stats, curated_stats,
            sample_dir, database_fasta, sample_fastq, min_identity, min_mapped_reads,
            coverage_depth_threshold, coverage_breadth_threshold, threads
        )
        return None
    
    # Select best reference from curated set (using user-specified thresholds)
    if curated_stats:
        curated_best = max(curated_stats, key=lambda x: x["mapped_reads"])
        logger.info(f"\nBest reference from curated set (identity >= {min_identity}%, coverage_depth >= {coverage_depth_threshold}, coverage_breadth >= {coverage_breadth_threshold}):")
        logger.info(f"  Accession: {curated_best['accession']}")
        logger.info(f"  Mapped reads: {curated_best['mapped_reads']:,}")
        logger.info(f"  Identity: {curated_best['avg_identity']:.2f}%")
        logger.info(f"  Coverage depth: {curated_best['coverage_depth']:.2f}")
        logger.info(f"  Coverage breadth: {curated_best['coverage_breadth']:.4f}")
        best_ref = {"accession": curated_best["accession"], "description": curated_best["description"], "sequence": None}
        best_stats = {
            "mapped_reads": curated_best["mapped_reads"],
            "coverage_depth": curated_best["coverage_depth"],
            "coverage_breadth": curated_best["coverage_breadth"],
            "avg_identity": curated_best["avg_identity"]
        }
    else:
        logger.warning("No references in curated set. Using best from filtered set.")
    
    # Save results with curated JSON (re-mapping will happen inside save_results after deduplication)
    # Use sample_dir (which is already database-specific if multiple databases)
    save_results(sample_name, best_ref, best_stats, all_stats, filtered_stats, curated_stats, sample_dir, database_fasta, sample_fastq, min_identity, min_mapped_reads, coverage_depth_threshold, coverage_breadth_threshold, threads)
    
    logger.info(f"\n{'='*60}")
    logger.info("Sample processing complete!")
    if best_ref and best_stats:
        logger.info(f"  - Best reference: {best_ref['accession']}")
        logger.info(f"  - Mapped reads: {best_stats['mapped_reads']:,}")
        logger.info(f"  - Average identity: {best_stats['avg_identity']:.2f}%")
    else:
        logger.info("  - No reads mapped to any reference")
    logger.info(f"{'='*60}")
    
    return best_ref

def save_results(sample_name, best_ref, best_stats, all_stats, filtered_stats, curated_stats, output_dir, database_fasta, sample_fastq, min_identity, min_mapped_reads, coverage_depth_threshold, coverage_breadth_threshold, threads=1):
    """Save mapping statistics and best reference to JSON file and save best reference sequence as FASTA."""
    # output_dir is already the sample directory (or database-specific subdirectory if multiple databases)
    sample_dir = output_dir
    
    # Determine which database we're using based on the database file path
    # This is more reliable than parsing descriptions
    database_fasta_path = Path(database_fasta)
    database_path_str = str(database_fasta_path).lower()
    if "refseq" in database_path_str:
        actual_database_source = "RefSeq"
    elif "rvdb" in database_path_str:
        actual_database_source = "RVDB"
    else:
        actual_database_source = None  # Will use description-based detection as fallback
    
    # Get database directory (where database FASTA is located)
    database_dir = database_fasta_path.parent
    
    # Get taxonomy paths (stored in database directory, not output directory)
    cache_path = get_taxonomy_cache_path(database_dir)
    database_path = get_taxonomy_database_path(database_dir)
    segment_database_path = get_segment_database_path(database_dir)
    
    # Get NCBI taxonomy file paths for offline species lookup
    ncbi_dir = get_ncbi_files_dir(database_dir)
    nodes_file = ncbi_dir / "nodes.dmp"
    names_file = ncbi_dir / "names.dmp"
    
    # Check if comprehensive database exists and is not empty, if not build it from NCBI pre-made files
    database_exists = database_path.exists()
    database_empty = False
    if database_exists:
        try:
            # Check SQLite database
            if database_path.suffix == '.db':
                conn = sqlite3.connect(str(database_path))
                cursor = conn.cursor()
                cursor.execute('SELECT COUNT(*) FROM accession_to_organism')
                count = cursor.fetchone()[0]
                conn.close()
                database_empty = (count == 0)
            else:
                # Check JSON database (legacy)
                with open(database_path, 'r') as f:
                    data = json.load(f)
                    database_empty = len(data) == 0
        except (sqlite3.Error, json.JSONDecodeError, Exception) as e:
            logger.debug(f"Error checking database: {e}")
            database_empty = True
    
    if not database_exists or database_empty:
        if database_empty:
            logger.info("="*60)
            logger.info("Comprehensive taxonomy database exists but is empty. Rebuilding...")
        else:
            logger.info("="*60)
            logger.info("Comprehensive taxonomy database not found.")
            logger.info("Downloading NCBI pre-made taxonomy files (one-time download, ~few GB)...")
            logger.info(f"Files will be stored in: {database_dir / 'taxonomy_cache'}")
            logger.info("This will work offline forever after this download completes.")
            logger.info("="*60)
        try:
            build_taxonomy_database_from_ncbi_download(database_dir, database_fasta=database_fasta_path)
            logger.info("Taxonomy database built successfully from NCBI files!")
            # Reload the database path after building
            database_path = get_taxonomy_database_path(database_dir)
        except Exception as e:
            logger.error(f"Failed to build taxonomy database: {e}")
            logger.warning("Continuing without comprehensive taxonomy database (will use API calls)")
    
    # sample_dir is already set at the beginning of save_results (line 2149)
    # Don't overwrite it here - output_dir is already the final directory
    
    # Enrich curated_stats with organism names
    curated_stats_with_organism = enrich_stats_with_organism(
        curated_stats.copy() if curated_stats else [],
        cache_path=cache_path,
        database_path=database_path,
        nodes_file=nodes_file,
        names_file=names_file
    )
    
    # Debug: Log Lassa virus stats before deduplication
    lassa_before = [s for s in curated_stats_with_organism if "lassa" in s.get("organism", "").lower() or "lassa" in s.get("description", "").lower()]
    if lassa_before:
        logger.info(f"DEBUG: Found {len(lassa_before)} Lassa virus reference(s) BEFORE deduplication:")
        for l in lassa_before:
            logger.info(f"  - {l.get('accession')}: {l.get('mapped_reads', 0)} reads, organism: {l.get('organism', 'N/A')}")
    
    # Deduplicate curated_stats by organism (preferring RefSeq)
    # For segmented viruses, keeps best reference for EACH segment
    curated_stats_dedup = deduplicate_by_organism(curated_stats_with_organism, segment_database_path=segment_database_path)
    
    # Debug: Log Lassa virus stats after deduplication
    lassa_after = [s for s in curated_stats_dedup if "lassa" in s.get("organism", "").lower() or "lassa" in s.get("description", "").lower()]
    if lassa_after:
        logger.info(f"DEBUG: Found {len(lassa_after)} Lassa virus reference(s) AFTER deduplication:")
        for l in lassa_after:
            logger.info(f"  - {l.get('accession')}: {l.get('mapped_reads', 0)} reads, organism: {l.get('organism', 'N/A')}")
    elif lassa_before:
        logger.warning(f"DEBUG: Lassa virus was present before deduplication but disappeared after!")
    
    # Log segment information for segmented viruses
    segments_found = {}
    for stat in curated_stats_dedup:
        organism = stat.get("organism", "").strip()
        accession = stat.get("accession", "")
        segment = extract_segment_from_description(stat.get("description", ""), accession=accession, database_path=segment_database_path)
        if segment:
            if organism not in segments_found:
                segments_found[organism] = []
            segments_found[organism].append(segment)
    
    logger.info(f"Deduplicating curated_reference_stats by VIRAL SPECIES (preferring REFSEQ over GenBank)...")
    logger.info(f"  - For segmented viruses (e.g., Lassa, Influenza), keeping best reference for EACH segment")
    if segments_found:
        for organism, segments in segments_found.items():
            logger.info(f"  - {organism}: keeping segments {', '.join(sorted(set(segments)))}")
    logger.info(f"  - Curated stats: {len(curated_stats_with_organism)} -> {len(curated_stats_dedup)} (deduplicated by VIRAL_SPECIES+SEGMENT, filtered: identity >= {min_identity}%, coverage_depth >= {coverage_depth_threshold}, coverage_breadth >= {coverage_breadth_threshold})")
    
    # Deduplicate all_stats and filtered_stats by species
    # Create deep copies to avoid sharing references with curated_stats
    import copy
    all_stats_dedup = [copy.deepcopy(s) for s in deduplicate_by_species(all_stats)]
    filtered_stats_dedup = [copy.deepcopy(s) for s in deduplicate_by_species(filtered_stats)]
    logger.info(f"Deduplicating all_stats and filtered_stats by species (keeping best per species)...")
    logger.info(f"  - All stats: {len(all_stats)} -> {len(all_stats_dedup)} (deduplicated by species)")
    logger.info(f"  - Filtered stats: {len(filtered_stats)} -> {len(filtered_stats_dedup)} (deduplicated by species)")
    
    # Create a deep copy of curated_stats_dedup BEFORE threshold filtering
    # This preserves ALL deduplicated references (including those that might not meet thresholds)
    # This will be used for reference_selection.json (should have original counts, not re-mapped)
    curated_stats_dedup_before_filter = copy.deepcopy(curated_stats_dedup)
    
    # Filter curated_stats_dedup by thresholds again (after deduplication)
    # This filtered version is used for re-mapping and curated_descriptions.json
    curated_stats_dedup = [
        s for s in curated_stats_dedup
        if s["avg_identity"] >= min_identity
        and s.get("coverage_depth", 0.0) >= coverage_depth_threshold
        and s.get("coverage_breadth", 0.0) >= coverage_breadth_threshold
    ]
    
    # Use the pre-filtered version for reference_selection.json (shows all deduplicated references with reads)
    curated_stats_dedup_original = curated_stats_dedup_before_filter
    
    # NOW re-map reads to FINAL deduplicated references for accurate counts
    if curated_stats_dedup:
        logger.info(f"\n{'='*60}")
        logger.info("Re-mapping all reads to FINAL deduplicated references for accurate read counts...")
        logger.info(f"Final reference set: {len(curated_stats_dedup)} references (after organism deduplication)")
        
        # Extract FINAL deduplicated references from database
        selected_headers = [stat["description"] for stat in curated_stats_dedup]
        # Ensure sample_dir exists and use it directly (no extra sample_name subdirectory)
        sample_dir.mkdir(parents=True, exist_ok=True)
        selected_refs_fasta = sample_dir / f"{sample_name}_selected_references.fasta"
        extracted_count = extract_selected_references(Path(database_fasta), selected_headers, selected_refs_fasta)
        
        if extracted_count > 0:
            # Re-map all reads to FINAL selected references (with indexing for speed)
            remap_sam = sample_dir / f"{sample_name}_remapped.sam"
            updated_stats = remap_to_selected_references(
                sample_fastq,
                selected_refs_fasta,
                remap_sam,
                min_identity=min_identity,
                threads=threads
            )
            
            # Update curated_stats_dedup with accurate counts
            # Match by description first, then by accession if description doesn't match
            # Also try partial description matching in case of format differences
            # For low-input samples, preserve original stats if remapping causes breadth to drop below threshold
            total_reads = count_reads(sample_fastq)
            is_low_input = total_reads < 100000
            
            for stat in curated_stats_dedup:
                desc = stat["description"]
                acc = stat.get("accession", "")
                updated = False
                is_lassa = "lassa" in desc.lower() or "lassa" in stat.get("organism", "").lower() or acc == "AY628204.1"
                
                # Store original stats before updating
                original_breadth = stat.get("coverage_breadth", 0.0)
                original_reads = stat.get("mapped_reads", 0)
                original_depth = stat.get("coverage_depth", 0.0)
                original_identity = stat.get("avg_identity", 0.0)
                
                # Try matching by exact description first
                if desc in updated_stats:
                    old_reads = stat["mapped_reads"]
                    new_reads = updated_stats[desc]["mapped_reads"]
                    new_breadth = updated_stats[desc]["coverage_breadth"]
                    new_depth = updated_stats[desc]["coverage_depth"]
                    new_identity = updated_stats[desc]["avg_identity"]
                    
                    # For low-input samples: if original breadth met threshold but remapped doesn't,
                    # preserve original stats to avoid losing valid detections
                    if is_low_input and original_breadth >= coverage_breadth_threshold and new_breadth < coverage_breadth_threshold:
                        logger.warning(f"DEBUG: {acc} remapping dropped breadth from {original_breadth:.4f} to {new_breadth:.4f} (below threshold {coverage_breadth_threshold}). Preserving original stats for low-input sample.")
                        # Keep original stats - don't update
                        continue
                    
                    stat["mapped_reads"] = new_reads
                    stat["avg_identity"] = new_identity
                    stat["coverage_depth"] = new_depth
                    stat["coverage_breadth"] = new_breadth
                    logger.info(f"  {stat['accession']}: {old_reads:,} -> {stat['mapped_reads']:,} mapped reads")
                    updated = True
                else:
                    # Try matching by accession (in case description format differs slightly)
                    for updated_desc, updated_data in updated_stats.items():
                        updated_acc = extract_accession_from_header(updated_desc)
                        if updated_acc == acc and acc:  # Only match if accession is not empty
                            old_reads = stat["mapped_reads"]
                            stat["mapped_reads"] = updated_data["mapped_reads"]
                            stat["avg_identity"] = updated_data["avg_identity"]
                            stat["coverage_depth"] = updated_data["coverage_depth"]
                            stat["coverage_breadth"] = updated_data["coverage_breadth"]
                            logger.info(f"  {stat['accession']}: {old_reads:,} -> {stat['mapped_reads']:,} mapped reads (matched by accession)")
                            updated = True
                            break
                    
                    # If still not matched, try partial description matching (for cases where description format differs)
                    if not updated:
                        # Try to find a description that contains the accession or key parts of the description
                        desc_parts = desc.split("|")
                        desc_accession = extract_accession_from_header(desc)
                        
                        for updated_desc, updated_data in updated_stats.items():
                            # Check if accessions match
                            if desc_accession and desc_accession == extract_accession_from_header(updated_desc):
                                old_reads = stat["mapped_reads"]
                                stat["mapped_reads"] = updated_data["mapped_reads"]
                                stat["avg_identity"] = updated_data["avg_identity"]
                                stat["coverage_depth"] = updated_data["coverage_depth"]
                                stat["coverage_breadth"] = updated_data["coverage_breadth"]
                                logger.info(f"  {stat['accession']}: {old_reads:,} -> {stat['mapped_reads']:,} mapped reads (matched by accession from description)")
                                updated = True
                                break
                            
                            # Try matching by key parts of description (e.g., organism name, segment)
                            if len(desc_parts) >= 4 and len(updated_desc.split("|")) >= 4:
                                # Compare key parts: accession and organism/segment info
                                if (desc_parts[2] if len(desc_parts) > 2 else "") in updated_desc:
                                    # Additional check: make sure accessions match
                                    if desc_accession == extract_accession_from_header(updated_desc):
                                        old_reads = stat["mapped_reads"]
                                        stat["mapped_reads"] = updated_data["mapped_reads"]
                                        stat["avg_identity"] = updated_data["avg_identity"]
                                        stat["coverage_depth"] = updated_data["coverage_depth"]
                                        stat["coverage_breadth"] = updated_data["coverage_breadth"]
                                        logger.info(f"  {stat['accession']}: {old_reads:,} -> {stat['mapped_reads']:,} mapped reads (matched by partial description)")
                                        updated = True
                                        break
                
                if not updated:
                    # If not found in updated_stats, keep original stats (don't set to 0)
                    # Debug: Log available updated_stats keys for Lassa virus
                    if is_lassa:
                        logger.warning(f"  DEBUG: {stat['accession']} (Lassa virus) not found in re-mapped stats")
                        logger.warning(f"  DEBUG: Looking for description: {desc[:100]}")
                        logger.warning(f"  DEBUG: Looking for accession: {acc}")
                        logger.warning(f"  DEBUG: Available updated_stats keys (first 10):")
                        for i, key in enumerate(list(updated_stats.keys())[:10]):
                            logger.warning(f"    {i+1}. {key[:100]}")
                        # Check if any key contains the accession
                        matching_keys = [k for k in updated_stats.keys() if acc in k or acc.split('.')[0] in k]
                        if matching_keys:
                            logger.warning(f"  DEBUG: Found {len(matching_keys)} key(s) containing accession {acc}:")
                            for mk in matching_keys[:3]:
                                logger.warning(f"    - {mk[:100]}")
                    # This can happen if the reference wasn't extracted correctly or description format differs
                    if is_lassa:
                        logger.warning(f"  {stat['accession']}: No updated stats found - keeping original stats. Description '{desc[:60]}...' not in updated_stats")
                    # Don't modify stat - keep original mapped_reads, identity, etc.
                    # The reference will appear with original stats in curated_descriptions.json
                elif is_lassa:
                    # Debug: Log Lassa virus stats after remapping
                    logger.info(f"  DEBUG: Lassa virus {stat['accession']} stats after remapping:")
                    logger.info(f"    - mapped_reads: {stat['mapped_reads']:,}")
                    logger.info(f"    - avg_identity: {stat['avg_identity']:.2f}%")
                    logger.info(f"    - coverage_depth: {stat['coverage_depth']:.2f}")
                    logger.info(f"    - coverage_breadth: {stat['coverage_breadth']:.4f}")
            
            logger.info(f"Re-mapping complete. Updated counts for {len([s for s in curated_stats_dedup if s['description'] in updated_stats])} references.")
            
            # Note: We do NOT update best_ref and best_stats here - keep the original initial mapping stats
            # Only curated_stats_dedup gets updated with re-mapped stats (for curated_descriptions.json)
        else:
            logger.warning("Could not extract final selected references. Skipping re-mapping.")
    else:
        logger.warning("No curated references meeting thresholds after deduplication. Skipping re-mapping.")
        # Clean up initial SAM file even if no curated references
        initial_sam = sample_dir / f"{sample_name}.sam"
        if initial_sam.exists():
            initial_sam.unlink()
            logger.info(f"Removed initial mapping SAM file: {initial_sam.name} (no curated references)")
    
    # Save curated descriptions JSON (simplified for easy review)
    # Filter by identity, coverage_depth, and coverage_breadth thresholds
    # Note: min_identity is database-specific: 97% for RefSeq (to reduce false positives),
    #       80% for RVDB (which works well with its curated references)
    curated_descriptions = []
    
    # Debug: Log Lassa virus stats before final filtering
    lassa_before_final = [s for s in curated_stats_dedup if "lassa" in s.get("organism", "").lower() or "lassa" in s.get("description", "").lower() or s.get("accession", "") == "AY628204.1"]
    if lassa_before_final:
        logger.info(f"DEBUG: Lassa virus BEFORE final threshold filtering:")
        for l in lassa_before_final:
            logger.info(f"  - {l.get('accession')}: identity={l.get('avg_identity', 0):.2f}% (>= {min_identity}%), "
                       f"depth={l.get('coverage_depth', 0):.2f} (>= {coverage_depth_threshold}), "
                       f"breadth={l.get('coverage_breadth', 0):.4f} (>= {coverage_breadth_threshold})")
    
    for stat in curated_stats_dedup:
        identity = stat.get("avg_identity", 0.0)
        depth = stat.get("coverage_depth", 0.0)
        breadth = stat.get("coverage_breadth", 0.0)
        
        # Debug: Log why Lassa virus is being filtered out
        is_lassa = "lassa" in stat.get("organism", "").lower() or "lassa" in stat.get("description", "").lower() or stat.get("accession", "") == "AY628204.1"
        
        if (identity >= min_identity and
            depth >= coverage_depth_threshold and
            breadth >= coverage_breadth_threshold):
            
            # Extract segment information
            accession = stat.get("accession", "")
            description = stat.get("description", "")
            segment = extract_segment_from_description(description, accession=accession, database_path=segment_database_path)
            
            # Format segment for display (e.g., "L-segment", "S-segment", "1-segment", etc.)
            segment_display = ""
            if segment:
                # Check if segment is a single letter (L, S, M) or a number/name
                if len(segment) == 1 and segment.isalpha():
                    segment_display = f"{segment}-segment"
                else:
                    segment_display = f"{segment}-segment"
            
            # Detect which database this reference came from
            # If we know which database we mapped against, use that (more reliable)
            if actual_database_source:
                database_source = actual_database_source
            else:
                # Fallback to description-based detection if database path doesn't indicate source
                database_source = "Unknown"
                desc_lower = description.lower()
                
                # Check for RefSeq patterns (more comprehensive)
                if (is_refseq(description) or 
                    "refseq" in desc_lower or 
                    description.startswith("NC_") or 
                    description.startswith("NM_") or
                    description.startswith("NZ_") or
                    "|refseq|" in desc_lower):
                    database_source = "RefSeq"
                # Check for RVDB/GenBank patterns
                elif ("|GENBANK|" in description or 
                      description.startswith("acc|GENBANK|") or
                      "genbank" in desc_lower or
                      "gb_" in desc_lower or
                      "gi|" in desc_lower or
                      description.startswith("acc|") and "|" in description and not "refseq" in desc_lower):
                    database_source = "RVDB"
            
            curated_descriptions.append({
                "accession": accession,
                "description": description,
                "organism": stat.get("organism", ""),
                "viral_species": stat.get("viral_species", ""),  # Include viral species
                "mapped_reads": stat.get("mapped_reads", 0),
                "avg_identity": stat.get("avg_identity", 0.0),
                "coverage_depth": stat.get("coverage_depth", 0.0),
                "coverage_breadth": stat.get("coverage_breadth", 0.0),
                "segment": segment_display if segment_display else None,  # None for non-segmented viruses
                "database_source": database_source,  # Indicates which database this reference came from
            })
        elif is_lassa:
            # Debug: Log why Lassa virus was filtered out
            reasons = []
            if identity < min_identity:
                reasons.append(f"identity {identity:.2f}% < {min_identity}%")
            if depth < coverage_depth_threshold:
                reasons.append(f"depth {depth:.2f} < {coverage_depth_threshold}")
            if breadth < coverage_breadth_threshold:
                reasons.append(f"breadth {breadth:.4f} < {coverage_breadth_threshold}")
            logger.warning(f"DEBUG: Lassa virus {stat.get('accession')} filtered out from final JSON: {', '.join(reasons)}")
    
    # Debug: Log Lassa virus stats after final filtering
    lassa_after_final = [s for s in curated_descriptions if "lassa" in s.get("organism", "").lower() or "lassa" in s.get("description", "").lower() or s.get("accession", "") == "AY628204.1"]
    if lassa_after_final:
        logger.info(f"DEBUG: Lassa virus AFTER final threshold filtering: {len(lassa_after_final)} entry(ies) in final JSON")
    elif lassa_before_final:
        logger.warning(f"DEBUG: Lassa virus was present before final filtering but is missing from final JSON!")
    
    if curated_descriptions:
        # Separate by database source if multiple databases were used
        curated_by_database = {}
        for desc in curated_descriptions:
            db_source = desc.get("database_source", "Unknown")
            if db_source not in curated_by_database:
                curated_by_database[db_source] = []
            curated_by_database[db_source].append(desc)
        
        # Check if we're already in a database-specific folder (sample/RVDB/ or sample/RefSeq/)
        # If so, don't create subdirectories - just save files directly
        sample_dir_str = str(sample_dir)
        sample_dir_parts = sample_dir.parts
        # Check if the last part of the path is a database name
        is_already_db_specific = (len(sample_dir_parts) > 0 and 
                                  sample_dir_parts[-1] in ['RVDB', 'RefSeq', 'Custom']) or \
                                 sample_dir_str.endswith('/RVDB') or sample_dir_str.endswith('/RefSeq') or \
                                 sample_dir_str.endswith('\\RVDB') or sample_dir_str.endswith('\\RefSeq') or \
                                 '/RVDB/' in sample_dir_str or '/RefSeq/' in sample_dir_str or \
                                 '\\RVDB\\' in sample_dir_str or '\\RefSeq\\' in sample_dir_str
        
        # Only create subdirectories if both RVDB and RefSeq are present AND we're not already in a database-specific folder
        known_databases = [db for db in curated_by_database.keys() if db in ["RVDB", "RefSeq"]]
        if len(known_databases) > 1 and not is_already_db_specific:
            # Multiple databases - save separate files for each (no combined file)
            for db_source in known_databases:
                descs = curated_by_database[db_source]
                db_dir = sample_dir / db_source
                db_dir.mkdir(parents=True, exist_ok=True)
                curated_json_file = db_dir / f"{sample_name}_final_selected_references.json"
                with open(curated_json_file, 'w') as f:
                    json.dump(descs, f, indent=2)
                logger.info(f"Saved {db_source} final selected references ({len(descs)} entries) to {curated_json_file}")
                
                # Create per-reference outputs for this database
                db_selected_refs = db_dir / f"{sample_name}_selected_references.fasta"
                db_remap_sam = db_dir / f"{sample_name}_remapped.sam"
                if db_selected_refs.exists() and db_remap_sam.exists():
                    create_per_reference_outputs(sample_name, descs, db_selected_refs, db_remap_sam, sample_fastq, db_dir)
                    
                    # Clean up: Remove selected_references.fasta and its index (references are in per-reference folders)
                    if db_selected_refs.exists():
                        db_selected_refs.unlink()
                        logger.info(f"Removed {db_selected_refs.name} (references are in per-reference folders)")
                    db_selected_refs_index = db_dir / f"{sample_name}_selected_references.fasta.mmi"
                    if db_selected_refs_index.exists():
                        db_selected_refs_index.unlink()
                        logger.info(f"Removed {db_selected_refs_index.name}")
        else:
            # Single database or already in database-specific folder - save directly in sample_dir
            curated_json_file = sample_dir / f"{sample_name}_final_selected_references.json"
        with open(curated_json_file, 'w') as f:
            json.dump(curated_descriptions, f, indent=2)
            logger.info(f"Saved final selected references ({len(curated_descriptions)} entries total) to {curated_json_file}")
            
            # Create per-reference outputs
            selected_refs_fasta = sample_dir / f"{sample_name}_selected_references.fasta"
            remap_sam = sample_dir / f"{sample_name}_remapped.sam"
            if selected_refs_fasta.exists() and remap_sam.exists():
                create_per_reference_outputs(sample_name, curated_descriptions, selected_refs_fasta, remap_sam, sample_fastq, sample_dir)
    
                # Clean up: Remove selected_references.fasta and its index (references are in per-reference folders)
                if selected_refs_fasta.exists():
                    selected_refs_fasta.unlink()
                    logger.info(f"Removed {selected_refs_fasta.name} (references are in per-reference folders)")
                selected_refs_index = sample_dir / f"{sample_name}_selected_references.fasta.mmi"
                if selected_refs_index.exists():
                    selected_refs_index.unlink()
                    logger.info(f"Removed {selected_refs_index.name}")
            else:
                # Remapping didn't complete - clean up initial SAM file anyway
                initial_sam = sample_dir / f"{sample_name}.sam"
                if initial_sam.exists():
                    initial_sam.unlink()
                    logger.info(f"Removed initial mapping SAM file: {initial_sam.name} (remapping files not found)")
    else:
        # No curated descriptions - clean up initial SAM file
        initial_sam = sample_dir / f"{sample_name}.sam"
        if initial_sam.exists():
            initial_sam.unlink()
            logger.info(f"Removed initial mapping SAM file: {initial_sam.name} (no curated references meeting thresholds)")
    
    # Note: best_reference.fasta is not created - use selected_references.fasta instead
    # The JSON file identifies which reference is the best (highest mapped_reads)
    # If you need just the best one, extract it from selected_references.fasta using the accession from the JSON
    
    # Detect database source for best reference
    best_database_source = "Unknown"
    if best_ref:
        best_desc = best_ref.get("description", "")
        if is_refseq(best_desc):
            best_database_source = "RefSeq"
        elif "|GENBANK|" in best_desc or best_desc.startswith("acc|GENBANK|"):
            best_database_source = "RVDB"
        else:
            desc_lower = best_desc.lower()
            if "refseq" in desc_lower:
                best_database_source = "RefSeq"
            elif any(term in desc_lower for term in ["genbank", "gb_", "gi|"]):
                best_database_source = "RVDB"
    
    # Save main results JSON
    results = {
        "sample_name": sample_name,
        "database_used": str(database_fasta_path.name),  # Show which database file was used
        "best_reference": {
            "accession": best_ref.get("accession", "") if best_ref else "",
            "description": best_ref.get("description", "") if best_ref else "",
            "database_source": best_database_source,  # Which database this reference came from
        },
        "best_stats": {
            "mapped_reads": best_stats.get("mapped_reads", 0) if best_stats else 0,
            "avg_identity": best_stats.get("avg_identity", 0.0) if best_stats else 0.0,
            "coverage_depth": best_stats.get("coverage_depth", 0.0) if best_stats else 0.0,
            "coverage_breadth": best_stats.get("coverage_breadth", 0.0) if best_stats else 0.0,
            "note": "These stats are from the initial mapping to the full database (reads may be split across multiple references). For accurate re-mapped counts, see curated_descriptions.json"
        },
        "all_stats": all_stats_dedup,
        "filtered_stats": filtered_stats_dedup,
        "curated_reference_stats": curated_stats_dedup_original,  # Use original stats (before re-mapping)
        "filtering_criteria": {
            "min_identity": min_identity,
            "min_mapped_reads": min_mapped_reads,
            "coverage_depth_threshold": coverage_depth_threshold,
            "coverage_breadth_threshold": coverage_breadth_threshold,
        }
    }
    
    # Separate results by database source if multiple databases were used
    # Group stats by database source
    all_stats_by_db = {}
    filtered_stats_by_db = {}
    curated_stats_by_db = {}
    
    def detect_database_source(description: str) -> str:
        """Detect database source from description with improved patterns."""
        # If we know which database we mapped against, use that (more reliable)
        if actual_database_source:
            return actual_database_source
        
        # Fallback to description-based detection if database path doesn't indicate source
        if not description:
            return "Unknown"
        desc_lower = description.lower()
        
        # Check for RefSeq patterns (more comprehensive)
        if (is_refseq(description) or 
            "refseq" in desc_lower or 
            description.startswith("NC_") or 
            description.startswith("NM_") or
            description.startswith("NZ_") or
            "|refseq|" in desc_lower):
            return "RefSeq"
        # Check for RVDB/GenBank patterns
        elif ("|GENBANK|" in description or 
              description.startswith("acc|GENBANK|") or
              "genbank" in desc_lower or
              "gb_" in desc_lower or
              "gi|" in desc_lower or
              (description.startswith("acc|") and "|" in description and "refseq" not in desc_lower)):
            return "RVDB"
        else:
            return "Unknown"
    
    for stat in all_stats_dedup:
        desc = stat.get("description", "")
        db_source = detect_database_source(desc)
        if db_source not in all_stats_by_db:
            all_stats_by_db[db_source] = []
        all_stats_by_db[db_source].append(stat)
    
    for stat in filtered_stats_dedup:
        desc = stat.get("description", "")
        db_source = detect_database_source(desc)
        if db_source not in filtered_stats_by_db:
            filtered_stats_by_db[db_source] = []
        filtered_stats_by_db[db_source].append(stat)
    
    for stat in curated_stats_dedup_original:
        desc = stat.get("description", "")
        db_source = detect_database_source(desc)
        if db_source not in curated_stats_by_db:
            curated_stats_by_db[db_source] = []
        curated_stats_by_db[db_source].append(stat)
    
    # Check if we're already in a database-specific folder (sample/RVDB/ or sample/RefSeq/)
    # If so, don't create subdirectories - just save files directly
    sample_dir_str = str(sample_dir)
    is_already_db_specific = sample_dir_str.endswith('/RVDB') or sample_dir_str.endswith('/RefSeq') or sample_dir_str.endswith('\\RVDB') or sample_dir_str.endswith('\\RefSeq')
    
    # Only create subdirectories if both RVDB and RefSeq are present AND we're not already in a database-specific folder
    all_db_sources = set(list(all_stats_by_db.keys()) + list(filtered_stats_by_db.keys()) + list(curated_stats_by_db.keys()))
    known_databases = [db for db in all_db_sources if db in ["RVDB", "RefSeq"]]
    
    # Save separate JSON files for each database (only if both RVDB and RefSeq are detected AND not already in db folder)
    if len(known_databases) > 1 and not is_already_db_specific:
        for db_source in known_databases:
            db_dir = sample_dir / db_source
            db_dir.mkdir(parents=True, exist_ok=True)
            
            # Find best reference for this database
            db_best_ref = None
            db_best_stats = None
            if curated_stats_by_db.get(db_source):
                db_best = max(curated_stats_by_db[db_source], key=lambda x: x.get("mapped_reads", 0))
                db_best_ref = {
                    "accession": db_best.get("accession", ""),
                    "description": db_best.get("description", ""),
                    "database_source": db_source
                }
                db_best_stats = {
                    "mapped_reads": db_best.get("mapped_reads", 0),
                    "avg_identity": db_best.get("avg_identity", 0.0),
                    "coverage_depth": db_best.get("coverage_depth", 0.0),
                    "coverage_breadth": db_best.get("coverage_breadth", 0.0),
                }
            
            db_results = {
                "sample_name": sample_name,
                "database_source": db_source,
                "database_used": str(database_fasta_path.name),
                "best_reference": db_best_ref if db_best_ref else {
                    "accession": "",
                    "description": "",
                    "database_source": db_source
                },
                "best_stats": db_best_stats if db_best_stats else {
                    "mapped_reads": 0,
                    "avg_identity": 0.0,
                    "coverage_depth": 0.0,
                    "coverage_breadth": 0.0,
                },
                "all_stats": all_stats_by_db.get(db_source, []),
                "filtered_stats": filtered_stats_by_db.get(db_source, []),
                "curated_reference_stats": curated_stats_by_db.get(db_source, []),
                "filtering_criteria": {
                    "min_identity": min_identity,
                    "min_mapped_reads": min_mapped_reads,
                    "coverage_depth_threshold": coverage_depth_threshold,
                    "coverage_breadth_threshold": coverage_breadth_threshold,
                }
            }
            
            db_json_file = db_dir / f"{sample_name}_unfiltered_all_references.json"
            with open(db_json_file, 'w') as f:
                json.dump(db_results, f, indent=2)
            logger.info(f"Saved {db_source} unfiltered results ({len(curated_stats_by_db.get(db_source, []))} curated references) to {db_json_file}")
    else:
        # Single database - save in main sample directory (no subdirectories)
        json_file = sample_dir / f"{sample_name}_unfiltered_all_references.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Saved unfiltered results to {json_file}")
    
    # Final safety check: Ensure initial SAM file is always cleaned up
    initial_sam = sample_dir / f"{sample_name}.sam"
    if initial_sam.exists():
        initial_sam.unlink()
        logger.info(f"Removed initial mapping SAM file: {initial_sam.name} (final cleanup)")
    
    # Print summary
    if best_ref and best_stats:
        logger.info(f"\nBest reference: {best_ref['accession']}")
        logger.info(f"  Description: {best_ref['description']}")
        logger.info(f"  Mapped reads: {best_stats['mapped_reads']:,}")
        logger.info(f"  Average identity: {best_stats['avg_identity']:.2f}%")
        logger.info(f"  Coverage depth: {best_stats.get('coverage_depth', 0.0):.2f}")
        logger.info(f"  Coverage breadth: {best_stats.get('coverage_breadth', 0.0):.4f}")

def generate_html_visualization(output_dir: Path):
    """
    Generate separate comprehensive HTML reports for each database with:
    - Filterable and sortable tables per sample
    - Interactive plots (linked to table order)
    - Download functionality (CSV, PNG, PDF)
    - Heatmap showing viral species across all samples
    - Pipeline configuration display
    - Sample toggle functionality
    
    Creates separate HTML files: results_summary_RefSeq.html, results_summary_RVDB.html, etc.
    """
    output_dir = Path(output_dir)
    
    # Find all final_selected_references.json files
    final_json_files = list(output_dir.rglob("*_final_selected_references.json"))
    # Find all unfiltered JSON files for metadata
    unfiltered_json_files = list(output_dir.rglob("*_unfiltered_all_references.json"))
    
    # Find all sample directories
    all_sample_dirs = [d for d in output_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
    all_sample_names = set()
    for sample_dir in all_sample_dirs:
        has_json = any(sample_dir.glob("*.json"))
        # Check for database subdirectories (RefSeq, RVDB, Custom, or any other subdirectory with JSON files)
        has_db_subdirs = any((sample_dir / subdir).is_dir() and any((sample_dir / subdir).glob("*.json")) for subdir in sample_dir.iterdir() if (sample_dir / subdir).is_dir())
        if has_json or has_db_subdirs:
            all_sample_names.add(sample_dir.name)
    
    for json_file in unfiltered_json_files:
        rel_path = json_file.relative_to(output_dir)
        parts = rel_path.parts
        if len(parts) >= 1:
            all_sample_names.add(parts[0])
    
    if not all_sample_names:
        logger.warning("No samples found for HTML visualization")
        return
    
    # Collect data from all JSON files
    samples_data = {}
    samples_metadata = {}
    samples_with_results = set()
    
    # Process final selected references
    for json_file in sorted(final_json_files):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            rel_path = json_file.relative_to(output_dir)
            parts = rel_path.parts
            
            if len(parts) == 2:
                # Single database case: file is directly in sample directory
                sample_name = parts[0]
                # Try to detect database from file path first
                database = "Unknown"
                json_path_str = str(json_file).lower()
                if "refseq" in json_path_str:
                    database = "RefSeq"
                elif "rvdb" in json_path_str:
                    database = "RVDB"
                else:
                    # Try to get from JSON data metadata if available
                    db_source_found = None
                    if isinstance(data, dict):
                        db_source_found = data.get("database_source", "")
                    elif isinstance(data, list) and len(data) > 0:
                        # Check items for database_source - look at first few items
                        for item in data[:5]:  # Check first 5 items
                            if isinstance(item, dict):
                                db_source = item.get("database_source", "")
                                if db_source:
                                    db_source_found = db_source
                                    break
                    
                    if db_source_found:
                        if "refseq" in db_source_found.lower():
                            database = "RefSeq"
                        elif "rvdb" in db_source_found.lower() or "genbank" in db_source_found.lower():
                            database = "RVDB"
                    # If still unknown, will be processed later from unfiltered JSON metadata
            elif len(parts) == 3:
                # Multiple databases case: file is in database subdirectory
                sample_name = parts[0]
                database = parts[1]
            else:
                sample_name = parts[0] if parts else "Unknown"
                database = "Unknown"
            
            database_lower = database.lower()
            if "refseq" in database_lower:
                database = "RefSeq"
            elif "rvdb" in database_lower:
                database = "RVDB"
            elif database == "Unknown" or "custom" in database_lower:
                # Try to detect from JSON data metadata first, then file path
                json_path_str = str(json_file).lower()
                db_source_from_data = None
                
                # Check JSON data for database_source
                if isinstance(data, dict):
                    db_source_from_data = data.get("database_source", "")
                elif isinstance(data, list) and len(data) > 0:
                    # Check first few items for database_source
                    for item in data[:5]:
                        if isinstance(item, dict):
                            db_source = item.get("database_source", "")
                            if db_source:
                                db_source_from_data = db_source
                                break
                
                if db_source_from_data:
                    if "refseq" in db_source_from_data.lower():
                        database = "RefSeq"
                    elif "rvdb" in db_source_from_data.lower() or "genbank" in db_source_from_data.lower():
                        database = "RVDB"
                    elif "refseq" in json_path_str:
                        database = "RefSeq"
                    elif "rvdb" in json_path_str:
                        database = "RVDB"
                    else:
                        database = "Custom"
                elif "refseq" in json_path_str:
                    database = "RefSeq"
                elif "rvdb" in json_path_str:
                    database = "RVDB"
                else:
                    # Try one more time to check JSON data for database_source
                    # This handles cases where the first check didn't find it
                    if isinstance(data, list) and len(data) > 0:
                        for item in data[:10]:  # Check more items
                            if isinstance(item, dict):
                                db_source = item.get("database_source", "")
                                if db_source:
                                    if "refseq" in db_source.lower():
                                        database = "RefSeq"
                                        break
                                    elif "rvdb" in db_source.lower() or "genbank" in db_source.lower():
                                        database = "RVDB"
                                        break
                    # If still unknown after all checks, default to Custom
                    if database == "Unknown":
                        database = "Custom"
            else:
                # Keep the database name as-is (could be Custom or other custom name)
                pass
            
            logger.debug(f"Processing JSON file: {json_file}, sample: {sample_name}, database: {database}")
            
            if sample_name not in samples_data:
                samples_data[sample_name] = {}
            
            if database not in samples_data[sample_name]:
                samples_data[sample_name][database] = []
            
            if isinstance(data, list):
                samples_data[sample_name][database].extend(data)
            else:
                samples_data[sample_name][database].append(data)
            
            samples_with_results.add(sample_name)
                
        except Exception as e:
            logger.warning(f"Could not read {json_file}: {e}")
            continue
    
    # Process unfiltered JSON files for metadata
    for json_file in sorted(unfiltered_json_files):
        try:
            with open(json_file, 'r') as f:
                metadata = json.load(f)
            
            rel_path = json_file.relative_to(output_dir)
            parts = rel_path.parts
            
            if len(parts) == 2:
                # Single database case: file is directly in sample directory
                sample_name = parts[0]
                # Try to detect database from metadata or file path
                db_source = metadata.get("database_source", "")
                db_used = metadata.get("database_used", "")
                json_path_str = str(json_file).lower()
                
                # Check database_source first
                if db_source and "refseq" in db_source.lower():
                    database = "RefSeq"
                elif db_source and ("rvdb" in db_source.lower() or "genbank" in db_source.lower()):
                    database = "RVDB"
                # Check database_used field (filename) as fallback
                elif db_used and "refseq" in db_used.lower():
                    database = "RefSeq"
                elif db_used and ("rvdb" in db_used.lower() or "genbank" in db_used.lower()):
                    database = "RVDB"
                # Check file path
                elif "refseq" in json_path_str:
                    database = "RefSeq"
                elif "rvdb" in json_path_str:
                    database = "RVDB"
                else:
                    database = "Unknown"
            elif len(parts) == 3:
                # Multiple databases case: file is in database subdirectory
                sample_name = parts[0]
                database = parts[1]
            else:
                sample_name = parts[0] if parts else "Unknown"
                database = "Unknown"
            
            database_lower = database.lower()
            if "refseq" in database_lower:
                database = "RefSeq"
            elif "rvdb" in database_lower:
                database = "RVDB"
            elif database == "Unknown" or "custom" in database_lower:
                # Try to detect from metadata or file path
                db_source = metadata.get("database_source", "")
                db_used = metadata.get("database_used", "")
                json_path_str = str(json_file).lower()
                
                # Check database_source first
                if db_source and "refseq" in db_source.lower():
                    database = "RefSeq"
                elif db_source and ("rvdb" in db_source.lower() or "genbank" in db_source.lower()):
                    database = "RVDB"
                # Check database_used field (filename) as fallback
                elif db_used and "refseq" in db_used.lower():
                    database = "RefSeq"
                elif db_used and ("rvdb" in db_used.lower() or "genbank" in db_used.lower()):
                    database = "RVDB"
                # Check file path
                elif "refseq" in json_path_str:
                    database = "RefSeq"
                elif "rvdb" in json_path_str:
                    database = "RVDB"
                else:
                    # Only default to Custom if we're certain it's not RefSeq/RVDB
                    # If database_source is "Unknown", don't create Custom - skip it
                    if db_source and db_source.lower() not in ["unknown", ""] and "refseq" not in db_source.lower() and "rvdb" not in db_source.lower() and "genbank" not in db_source.lower():
                        database = "Custom"
                    else:
                        # Don't create Custom HTML for Unknown databases
                        database = "Unknown"
            else:
                # Keep the database name as-is (could be Custom or other custom name)
                pass
            
            if sample_name not in samples_metadata:
                samples_metadata[sample_name] = {}
            
            if database not in samples_metadata[sample_name]:
                db_source = metadata.get("database_source", database)
                if db_source and "refseq" in db_source.lower():
                    db_source = "RefSeq"
                elif db_source and ("rvdb" in db_source.lower() or "genbank" in db_source.lower()):
                    db_source = "RVDB"
                else:
                    db_source = database
                
                samples_metadata[sample_name][database] = {
                    "database_source": db_source,
                    "filtering_criteria": metadata.get("filtering_criteria", {})
                }
                
        except Exception as e:
            logger.warning(f"Could not read metadata from {json_file}: {e}")
            continue
    
    # Prepare data for JavaScript
    import json as json_module
    samples_json = {}
    for sample_name in all_sample_names:
        samples_json[sample_name] = {}
        if sample_name in samples_data:
            for database in samples_data[sample_name]:
                samples_json[sample_name][database] = {
                    "references": samples_data[sample_name][database],
                    "metadata": samples_metadata.get(sample_name, {}).get(database, {
                        "database_source": database,
                        "filtering_criteria": {}
                    })
                }
        else:
            if sample_name in samples_metadata:
                for database in samples_metadata[sample_name]:
                    samples_json[sample_name][database] = {
                        "references": [],
                        "metadata": samples_metadata[sample_name][database]
                    }
            else:
                samples_json[sample_name]["Unknown"] = {
                    "references": [],
                    "metadata": {
                        "database_source": "Unknown",
                        "filtering_criteria": {}
                    }
                }
    
    # Find all unique databases (exclude "Unknown")
    all_databases = set()
    for sample_name in samples_data:
        for db_name in samples_data[sample_name].keys():
            if db_name != "Unknown":
                all_databases.add(db_name)
    for sample_name in samples_metadata:
        for db_name in samples_metadata[sample_name].keys():
            if db_name != "Unknown":
                all_databases.add(db_name)
    
    if not all_databases:
        logger.warning("No databases found for HTML visualization")
        return
    
    # Generate separate HTML file for each database
    for database_name in sorted(all_databases):
        # Filter data for this database only
        db_samples_data = {}
        db_samples_metadata = {}
        db_samples_with_results = set()
        
        for sample_name in all_sample_names:
            if sample_name in samples_data and database_name in samples_data[sample_name]:
                db_samples_data[sample_name] = {database_name: samples_data[sample_name][database_name]}
                db_samples_with_results.add(sample_name)
            if sample_name in samples_metadata and database_name in samples_metadata[sample_name]:
                db_samples_metadata[sample_name] = {database_name: samples_metadata[sample_name][database_name]}
        
        # Build heatmap data for this database only
        heatmap_data = {}  # species -> {sample_name: True/False or count}
        all_species = set()
        for sample_name in db_samples_data:
            if database_name in db_samples_data[sample_name]:
                for ref in db_samples_data[sample_name][database_name]:
                    species = ref.get('viral_species', '') or ref.get('organism', 'Unknown')
                    if species:
                        all_species.add(species)
                        if species not in heatmap_data:
                            heatmap_data[species] = {}
                        heatmap_data[species][sample_name] = True
        
        # Prepare data for JavaScript (only this database)
        # Only include samples that have data for this database
        db_samples_json = {}
        db_samples_list = []
        for sample_name in all_sample_names:
            has_data = False
            if sample_name in db_samples_data and database_name in db_samples_data[sample_name]:
                db_samples_json[sample_name] = {
                    database_name: {
                        "references": db_samples_data[sample_name][database_name],
                        "metadata": db_samples_metadata.get(sample_name, {}).get(database_name, {
                            "database_source": database_name,
                            "filtering_criteria": {}
                        })
                    }
                }
                has_data = True
                db_samples_list.append(sample_name)
            elif sample_name in db_samples_metadata and database_name in db_samples_metadata[sample_name]:
                db_samples_json[sample_name] = {
                    database_name: {
                        "references": [],
                        "metadata": db_samples_metadata[sample_name][database_name]
                    }
                }
                has_data = True
                db_samples_list.append(sample_name)
        
        # Generate HTML for this database
        samples_json_str = json_module.dumps(db_samples_json).replace('</', '<\\/')
        heatmap_data_str = json_module.dumps(heatmap_data).replace('</', '<\\/')
        all_samples_list = sorted(db_samples_list)
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Virasign Results Summary - {database_name}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 40px 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1800px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
            width: 100%;
            box-sizing: border-box;
        }}
        @media (min-width: 1840px) {{
            .container {{
                width: 1800px;
            }}
        }}
        .content {{
            width: 100%;
            box-sizing: border-box;
            max-width: 100%;
            overflow-x: hidden;
            padding: 30px;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        .table-section {{
            width: 100%;
            max-width: 100%;
            overflow-x: auto;
            box-sizing: border-box;
        }}
        .sample-selector {{
            background: #f8f9fa;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .sample-selector label {{
            font-size: 1.2em;
            font-weight: 600;
            color: #333;
            margin-right: 15px;
        }}
        .sample-selector select {{
            font-size: 1.1em;
            padding: 10px 20px;
            border: 2px solid #667eea;
            border-radius: 6px;
            background: white;
            color: #333;
            cursor: pointer;
            min-width: 300px;
        }}
        .summary {{
            background: #f8f9fa;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .summary h2 {{
            margin-bottom: 15px;
            color: #333;
        }}
        .summary-stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .stat-box {{
            background: white;
            padding: 15px;
            border-radius: 6px;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stat-box .number {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }}
        .stat-box .label {{
            color: #6c757d;
            margin-top: 5px;
        }}
        .sample-section {{
            display: none;
            margin-bottom: 40px;
        }}
        .sample-section.active {{
            display: block;
        }}
        .metadata-section {{
            background: #f8f9fa;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 8px;
            border-left: 4px solid #28a745;
        }}
        .metadata-section h3 {{
            margin-bottom: 15px;
            color: #333;
            font-size: 1.3em;
        }}
        .metadata-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .metadata-item {{
            background: white;
            padding: 15px;
            border-radius: 6px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metadata-item .label {{
            font-weight: 600;
            color: #6c757d;
            font-size: 0.9em;
            margin-bottom: 5px;
        }}
        .metadata-item .value {{
            font-size: 1.1em;
            color: #333;
            font-family: 'Courier New', monospace;
        }}
        .table-section {{
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .table-header {{
            background: #f8f9fa;
            padding: 15px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 2px solid #dee2e6;
        }}
        .table-header h3 {{
            margin: 0;
            color: #333;
        }}
        .filter-row {{
            background: #f8f9fa;
            padding: 10px;
            display: grid;
            grid-template-columns: 1.2fr 1.5fr 1.5fr 0.8fr 1fr 1fr 1fr 1fr auto;
            gap: 10px;
            border-bottom: 1px solid #dee2e6;
            align-items: center;
        }}
        .filter-row input[type="text"] {{
            padding: 5px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 0.9em;
        }}
        .filter-row .checkbox-container {{
            display: flex;
            align-items: center;
            gap: 5px;
            white-space: nowrap;
        }}
        .filter-row .checkbox-container input[type="checkbox"] {{
            width: 18px;
            height: 18px;
            cursor: pointer;
        }}
        .filter-row .checkbox-container label {{
            font-size: 0.9em;
            color: #495057;
            cursor: pointer;
            user-select: none;
        }}
        table {{
            width: 100%;
            max-width: 100%;
            border-collapse: collapse;
            background: white;
            table-layout: auto;
        }}
        th {{
            background: #f8f9fa;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            color: #495057;
            border-bottom: 2px solid #dee2e6;
            position: sticky;
            top: 0;
        }}
        th.sortable {{
            cursor: pointer;
            user-select: none;
            position: relative;
            padding-right: 25px;
        }}
        th.sortable:hover {{
            background: #e9ecef;
        }}
        th.sortable::after {{
            content: ' ↕';
            position: absolute;
            right: 8px;
            opacity: 0.5;
            font-size: 0.8em;
        }}
        th.sortable.asc::after {{
            content: ' ↑';
            opacity: 1;
        }}
        th.sortable.desc::after {{
            content: ' ↓';
            opacity: 1;
        }}
        td {{
            padding: 12px;
            border-bottom: 1px solid #e9ecef;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .accession {{
            font-family: 'Courier New', monospace;
            font-weight: bold;
            color: #1976d2;
        }}
        .stats {{
            text-align: right;
        }}
        .charts-section {{
            display: grid;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            gap: 15px;
            margin: 30px 0;
            width: 100%;
            box-sizing: border-box;
        }}
        .chart-container {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            min-width: 0;
            max-width: 100%;
            overflow: hidden;
            width: 100%;
            box-sizing: border-box;
        }}
        .chart-container canvas {{
            max-width: 100% !important;
            height: auto !important;
            width: 100% !important;
        }}
        .chart-container > div {{
            width: 100% !important;
            max-width: 100% !important;
            box-sizing: border-box;
        }}
        @media (max-width: 900px) {{
            .charts-section {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
        @media (max-width: 600px) {{
            .charts-section {{
                grid-template-columns: 1fr;
            }}
        }}
        .chart-container h3 {{
            margin-bottom: 15px;
            color: #333;
            text-align: center;
        }}
        .heatmap-section {{
            background: white;
            padding: 20px;
            margin: 30px 0;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .heatmap-section h2 {{
            margin-bottom: 20px;
            color: #333;
        }}
        .heatmap-container {{
            overflow-x: auto;
        }}
        .heatmap-table {{
            border-collapse: collapse;
            width: 100%;
        }}
        .heatmap-table th {{
            background: #667eea;
            color: white;
            padding: 10px;
            text-align: center;
            border: 1px solid #555;
        }}
        .heatmap-table td {{
            padding: 10px;
            text-align: center;
            border: 1px solid #ddd;
        }}
        .heatmap-cell {{
            padding: 8px;
            border-radius: 4px;
        }}
        .heatmap-present {{
            background: #28a745;
            color: white;
            font-weight: bold;
        }}
        .heatmap-absent {{
            background: #f8f9fa;
            color: #999;
        }}
        button {{
            background: #667eea;
            color: white;
            border: none;
            padding: 8px 15px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.9em;
        }}
        button:hover {{
            background: #764ba2;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🔬 Virasign Results Summary - {database_name}</h1>
            <p>Viral Read ASSIGNment from nanopore sequencing</p>
        </div>
        <div class="content">
            <div class="sample-selector">
                <label for="sampleSelect">Select Sample:</label>
                <select id="sampleSelect" onchange="showSample(this.value)">
                    <option value="">-- Select a sample --</option>
"""
        
        # Add sample options for this database
        for sample_name in all_samples_list:
            has_results = sample_name in db_samples_with_results
            display_name = f"{sample_name} {'(no viruses detected)' if not has_results else ''}"
            html_content += f'                    <option value="{sample_name}">{display_name}</option>\n'
        
        html_content += f"""                </select>
            </div>
            
            <div class="summary">
                <h2>Summary - {database_name} Database</h2>
                <div class="summary-stats">
                    <div class="stat-box">
                        <div class="number">{len(all_samples_list)}</div>
                        <div class="label">Sample(s)</div>
                    </div>
                    <div class="stat-box">
                        <div class="number">{database_name}</div>
                        <div class="label">Database</div>
                    </div>
                </div>
            </div>
"""
        
        # Add sample sections for this database
        for idx, sample_name in enumerate(all_samples_list):
            has_results = sample_name in db_samples_with_results
            active_class = " active" if idx == 0 else ""
            
            html_content += f"""
            <div class="sample-section{active_class}" id="sample-{sample_name}">
                <div class="metadata-section">
                    <h3>📋 Pipeline Configuration</h3>
                    <div class="metadata-grid" id="metadata-{sample_name}">
                        <!-- Metadata will be populated by JavaScript -->
                    </div>
                </div>
"""
            
            if not has_results:
                html_content += f"""
                <div style="background: #fff3cd; border: 2px solid #ffc107; border-radius: 8px; padding: 30px; margin: 20px 0; text-align: center;">
                    <h2 style="color: #856404;">⚠️ No Viruses Detected</h2>
                    <p style="color: #856404;">Sample: {sample_name}</p>
                </div>
            </div>
"""
            else:
                # Add table section
                html_content += f"""
                <div class="table-section">
                    <div class="table-header">
                        <h3>Viral References</h3>
                        <button onclick="downloadTableAsCSV('{sample_name}')">📥 Download Table (CSV)</button>
                    </div>
                    <div class="filter-row" id="filters-{sample_name}">
                        <input type="text" placeholder="Filter Accession" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Filter Organism" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Filter Species" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Filter Segment" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Minimum Reads" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Minimum Identity" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Minimum Depth" onkeyup="applyAllFilters('{sample_name}')">
                        <input type="text" placeholder="Minimum Breadth" onkeyup="applyAllFilters('{sample_name}')">
                        <div class="checkbox-container">
                            <input type="checkbox" id="excludeUnknown-{sample_name}" onchange="applyAllFilters('{sample_name}')">
                            <label for="excludeUnknown-{sample_name}">Exclude Unknown</label>
                        </div>
                    </div>
                    <table id="table-{sample_name}">
                        <thead>
                            <tr>
                                <th>Accession</th>
                                <th>Organism</th>
                                <th>Viral Species</th>
                                <th>Segment</th>
                                <th class="stats sortable" data-sort="mapped_reads">Mapped Reads</th>
                                <th class="stats sortable" data-sort="identity">Identity (%)</th>
                                <th class="stats sortable" data-sort="depth">Coverage Depth</th>
                                <th class="stats sortable" data-sort="breadth">Coverage Breadth</th>
                            </tr>
                        </thead>
                        <tbody id="tbody-{sample_name}">
                            <!-- Rows will be populated by JavaScript -->
                        </tbody>
                    </table>
                </div>
                
                <div class="charts-section" id="charts-{sample_name}">
                    <!-- Charts will be populated by JavaScript -->
                </div>
            </div>
"""
    
        # Add heatmap section at the bottom (after all sample sections)
        html_content += """
            
            <!-- Heatmap Section (All Samples) - At Bottom -->
            <div class="heatmap-section">
                <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
                    <h2>Viral Species Heatmap (All Samples)</h2>
                    <div style="display: flex; gap: 10px;">
                        <button onclick="downloadHeatmapAsPNG()">📥 Download PNG</button>
                    </div>
                </div>
                <div class="heatmap-container" id="heatmap-container">
                    <!-- Heatmap will be generated by JavaScript -->
                </div>
            </div>
"""
        
        # Add JavaScript
        html_content += f"""
        </div>
    </div>
    <script>
        const samplesData = {samples_json_str};
        const heatmapData = {heatmap_data_str};
        const allSamples = {json_module.dumps(all_samples_list)};
        let charts = {{}};
        let currentSortOrder = {{}};
        
        // Show sample function
        function showSample(sampleName) {{
            document.querySelectorAll('.sample-section').forEach(section => {{
                section.classList.remove('active');
            }});
            
            if (!sampleName) return;
            
            const sampleSection = document.getElementById(`sample-${{sampleName}}`);
            if (sampleSection) {{
                sampleSection.classList.add('active');
                populateMetadata(sampleName);
                populateTable(sampleName);
                if (!charts[sampleName]) {{
                    createCharts(sampleName);
                    charts[sampleName] = true;
                }}
            }}
        }}
        
        // Populate metadata
        function populateMetadata(sampleName) {{
            const container = document.getElementById(`metadata-${{sampleName}}`);
            if (!container || !samplesData[sampleName]) return;
            
            let html = '';
            const databases = Object.keys(samplesData[sampleName]);
            
            for (let db of databases) {{
                const meta = samplesData[sampleName][db].metadata;
                const criteria = meta.filtering_criteria || {{}};
                
                // Use actual values, or defaults if missing/zero
                const dbSource = meta.database_source || db;
                const isRefSeq = dbSource.toLowerCase().includes('refseq');
                const minId = (criteria.min_identity && criteria.min_identity > 0) ? criteria.min_identity : (isRefSeq ? 95 : 80);
                const minReads = (criteria.min_mapped_reads && criteria.min_mapped_reads > 0) ? criteria.min_mapped_reads : 100;
                const depthThresh = (criteria.coverage_depth_threshold && criteria.coverage_depth_threshold > 0) ? criteria.coverage_depth_threshold : 1.0;
                const breadthThresh = (criteria.coverage_breadth_threshold && criteria.coverage_breadth_threshold > 0) ? criteria.coverage_breadth_threshold : 0.1;
                
                html += `
                    <div class="metadata-item">
                        <div class="label">Database Source</div>
                        <div class="value">${{dbSource}}</div>
                    </div>
                    <div class="metadata-item">
                        <div class="label">Min Identity</div>
                        <div class="value">${{minId}}%</div>
                    </div>
                    <div class="metadata-item">
                        <div class="label">Min Mapped Reads</div>
                        <div class="value">${{minReads}}</div>
                    </div>
                    <div class="metadata-item">
                        <div class="label">Coverage Depth Threshold</div>
                        <div class="value">${{depthThresh}}x</div>
                    </div>
                    <div class="metadata-item">
                        <div class="label">Coverage Breadth Threshold</div>
                        <div class="value">${{(breadthThresh * 100).toFixed(1)}}%</div>
                    </div>
                `;
            }}
            
            container.innerHTML = html;
        }}
        
        // Populate table
        function populateTable(sampleName) {{
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            if (!tbody || !samplesData[sampleName]) return;
            
            let allRefs = [];
            Object.keys(samplesData[sampleName]).forEach(db => {{
                if (samplesData[sampleName][db].references) {{
                    allRefs = allRefs.concat(samplesData[sampleName][db].references);
                }}
            }});
            
            // Sort by coverage breadth (default)
            allRefs.sort((a, b) => (b.coverage_breadth || 0) - (a.coverage_breadth || 0));
            
            let html = '';
            allRefs.forEach(ref => {{
                const accession = ref.accession || 'N/A';
                const organism = ref.organism || 'Unknown';
                const species = ref.viral_species || organism || 'N/A';
                const description = ref.description || 'N/A';
                const segment = ref.segment || '';
                const reads = ref.mapped_reads || 0;
                const identity = ref.avg_identity || 0;
                const depth = ref.coverage_depth || 0;
                const breadth = ref.coverage_breadth || 0;
                
                html += `
                    <tr data-accession="${{accession}}" data-organism="${{organism}}" data-species="${{species}}" 
                        data-segment="${{segment}}"
                        data-reads="${{reads}}" data-identity="${{identity}}" data-depth="${{depth}}" data-breadth="${{breadth}}">
                        <td class="accession"><a href="https://www.ncbi.nlm.nih.gov/nuccore/${{accession}}" target="_blank">${{accession}}</a></td>
                        <td>${{organism}}</td>
                        <td>${{species}}</td>
                        <td>${{segment || '—'}}</td>
                        <td class="stats">${{reads.toLocaleString()}}</td>
                        <td class="stats">${{identity.toFixed(2)}}%</td>
                        <td class="stats">${{depth.toFixed(2)}}x</td>
                        <td class="stats">${{(breadth * 100).toFixed(1)}}%</td>
                    </tr>
                `;
            }});
            
            tbody.innerHTML = html;
            
            // Initialize sorting
            initSorting(sampleName);
        }}
        
        // Apply all filters (including exclude Unknown checkbox)
        function applyAllFilters(sampleName) {{
            const filterRow = document.getElementById(`filters-${{sampleName}}`);
            if (!filterRow) return;
            
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            if (!tbody) return;
            
            const rows = tbody.querySelectorAll('tr');
            const inputs = filterRow.querySelectorAll('input[type="text"]');
            const excludeUnknownCheckbox = document.getElementById(`excludeUnknown-${{sampleName}}`);
            const excludeUnknown = excludeUnknownCheckbox && excludeUnknownCheckbox.checked;
            
            // Get all filter values
            const filterValues = Array.from(inputs).map(input => input.value.trim());
            
            rows.forEach(row => {{
                const cells = row.querySelectorAll('td');
                let shouldShow = true;
                
                // Apply each column filter
                for (let colIndex = 0; colIndex < filterValues.length && colIndex < cells.length; colIndex++) {{
                    const filterValue = filterValues[colIndex];
                    const filterLower = filterValue.toLowerCase();
                    const filterTrimmed = filterValue;
                    
                    // Columns 4-7 are numeric (reads, identity, depth, breadth)
                    const isNumericColumn = colIndex >= 4 && colIndex <= 7;
                    
                    if (isNumericColumn) {{
                        if (filterTrimmed !== '') {{
                            const minValue = parseFloat(filterTrimmed);
                            if (!isNaN(minValue)) {{
                                let cellValue = null;
                                if (colIndex === 4) {{
                                    // Mapped Reads
                                    const dataReads = row.getAttribute('data-reads');
                                    if (dataReads !== null && dataReads !== '' && !isNaN(parseFloat(dataReads))) {{
                                        cellValue = parseFloat(dataReads);
                                    }} else {{
                                        const cellText = cells[colIndex].textContent.replace(/[,\s]/g, '').trim();
                                        cellValue = parseFloat(cellText);
                                        if (isNaN(cellValue)) {{
                                            const altReads = row.dataset.reads;
                                            if (altReads) {{
                                                cellValue = parseFloat(altReads);
                                            }}
                                        }}
                                    }}
                                }} else if (colIndex === 5) {{
                                    // Identity
                                    cellValue = parseFloat(row.getAttribute('data-identity') || 0);
                                }} else if (colIndex === 6) {{
                                    // Depth
                                    cellValue = parseFloat(row.getAttribute('data-depth') || 0);
                                }} else if (colIndex === 7) {{
                                    // Breadth
                                    cellValue = parseFloat(row.getAttribute('data-breadth') || 0) * 100;
                                }}
                                if (cellValue === null || isNaN(cellValue) || cellValue < minValue) {{
                                    shouldShow = false;
                                    break;
                                }}
                            }}
                        }}
                    }} else {{
                        // Text matching for other columns
                        if (filterLower !== '') {{
                            const cellText = cells[colIndex].textContent.toLowerCase();
                            if (!cellText.includes(filterLower)) {{
                                shouldShow = false;
                                break;
                            }}
                        }}
                    }}
                }}
                
                // Apply "Exclude Unknown" filter if checkbox is checked
                if (shouldShow && excludeUnknown && cells.length > 2) {{
                    const speciesCell = cells[2];
                    if (speciesCell) {{
                        const speciesText = speciesCell.textContent.trim().toLowerCase();
                        if (speciesText === 'unknown') {{
                            shouldShow = false;
                        }}
                    }}
                }}
                
                row.style.display = shouldShow ? '' : 'none';
            }});
            
            // Update charts immediately after filtering
            if (charts[sampleName]) {{
                updateCharts(sampleName);
            }}
        }}
        
        // Filter table (kept for backward compatibility, but now calls applyAllFilters)
        function filterTable(sampleName, colIndex, filterValue) {{
            applyAllFilters(sampleName);
        }}
        
        // Initialize sorting
        function initSorting(sampleName) {{
            const table = document.getElementById(`table-${{sampleName}}`);
            if (!table) return;
            
            const sortableHeaders = table.querySelectorAll('th.sortable');
            sortableHeaders.forEach(header => {{
                // Remove any existing event listeners by cloning
                const newHeader = header.cloneNode(true);
                header.parentNode.replaceChild(newHeader, header);
                
                newHeader.addEventListener('click', function(e) {{
                    e.preventDefault();
                    e.stopPropagation();
                    const sortType = this.getAttribute('data-sort');
                    const isAsc = this.classList.contains('asc');
                    const isDesc = this.classList.contains('desc');
                    
                    // Get current sort state for this column
                    const currentSortType = currentSortOrder[sampleName]?.type;
                    
                    // Remove all sort classes from all sortable headers in this table
                    const allHeaders = table.querySelectorAll('th.sortable');
                    allHeaders.forEach(th => th.classList.remove('asc', 'desc'));
                    
                    // Determine new sort order
                    let newOrder;
                    if (sortType === currentSortType) {{
                        // Same column: toggle between asc and desc
                        if (isAsc) {{
                            newOrder = 'desc';
                        }} else if (isDesc) {{
                            newOrder = 'asc';
                        }} else {{
                            newOrder = 'asc';
                        }}
                    }} else {{
                        // Different column: always start with ascending
                        newOrder = 'asc';
                    }}
                    
                    this.classList.add(newOrder);
                    
                    // Store current sort state
                    currentSortOrder[sampleName] = {{
                        type: sortType,
                        order: newOrder
                    }};
                    
                    // Sort table
                    sortTable(sampleName, sortType, newOrder);
                }});
            }});
        }}
        
        // Sort table
        function sortTable(sampleName, sortType, order) {{
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            if (!tbody) return;
            
            const rows = Array.from(tbody.querySelectorAll('tr:not([style*="display: none"])'));
            
            rows.sort((a, b) => {{
                let aVal, bVal;
                if (sortType === 'mapped_reads') {{
                    aVal = parseFloat(a.getAttribute('data-reads') || 0);
                    bVal = parseFloat(b.getAttribute('data-reads') || 0);
                }} else if (sortType === 'identity') {{
                    aVal = parseFloat(a.getAttribute('data-identity') || 0);
                    bVal = parseFloat(b.getAttribute('data-identity') || 0);
                }} else if (sortType === 'depth') {{
                    aVal = parseFloat(a.getAttribute('data-depth') || 0);
                    bVal = parseFloat(b.getAttribute('data-depth') || 0);
                }} else if (sortType === 'breadth') {{
                    aVal = parseFloat(a.getAttribute('data-breadth') || 0);
                    bVal = parseFloat(b.getAttribute('data-breadth') || 0);
                }} else {{
                    return 0;
                }}
                
                return order === 'asc' ? aVal - bVal : bVal - aVal;
            }});
            
            rows.forEach(row => tbody.appendChild(row));
            
            // Update charts
            updateCharts(sampleName);
        }}
        
        // Create charts
        function createCharts(sampleName) {{
            const container = document.getElementById(`charts-${{sampleName}}`);
            if (!container || !samplesData[sampleName]) return;
            
            let allRefs = [];
            Object.keys(samplesData[sampleName]).forEach(db => {{
                if (samplesData[sampleName][db].references) {{
                    allRefs = allRefs.concat(samplesData[sampleName][db].references);
                }}
            }});
            
            if (allRefs.length === 0) return;
            
            // Get current table order
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            const rows = tbody ? Array.from(tbody.querySelectorAll('tr:not([style*="display: none"])')) : [];
            
            let labels = [];
            let reads = [];
            let identities = [];
            let depths = [];
            let breadths = [];
            
            if (rows.length > 0) {{
                rows.forEach(row => {{
                    labels.push(row.getAttribute('data-species') || 'N/A');
                    reads.push(parseFloat(row.getAttribute('data-reads') || 0));
                    identities.push(parseFloat(row.getAttribute('data-identity') || 0));
                    depths.push(parseFloat(row.getAttribute('data-depth') || 0));
                    breadths.push(parseFloat(row.getAttribute('data-breadth') || 0) * 100);
                }});
            }} else {{
                allRefs.sort((a, b) => (b.coverage_breadth || 0) - (a.coverage_breadth || 0));
                allRefs.forEach(ref => {{
                    labels.push(ref.viral_species || ref.organism || ref.accession || 'N/A');
                    reads.push(ref.mapped_reads || 0);
                    identities.push(ref.avg_identity || 0);
                    depths.push(ref.coverage_depth || 0);
                    breadths.push((ref.coverage_breadth || 0) * 100);
                }});
            }}
            
            container.innerHTML = `
                <div class="chart-container">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                        <h3>Mapped Reads</h3>
                        <div>
                            <button onclick="downloadChart('${{sampleName}}', 'reads', 'png')">📥 Download PNG</button>
                        </div>
                    </div>
                    <div style="position: relative; width: 100%; height: 250px;">
                        <canvas id="chart-reads-${{sampleName}}"></canvas>
                    </div>
                </div>
                <div class="chart-container">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                        <h3>Identity (%)</h3>
                        <div>
                            <button onclick="downloadChart('${{sampleName}}', 'identity', 'png')">📥 Download PNG</button>
                        </div>
                    </div>
                    <div style="position: relative; width: 100%; height: 250px;">
                        <canvas id="chart-identity-${{sampleName}}"></canvas>
                    </div>
                </div>
                <div class="chart-container">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                        <h3>Coverage Depth</h3>
                        <div>
                            <button onclick="downloadChart('${{sampleName}}', 'depth', 'png')">📥 Download PNG</button>
                        </div>
                    </div>
                    <div style="position: relative; width: 100%; height: 250px;">
                        <canvas id="chart-depth-${{sampleName}}"></canvas>
                    </div>
                </div>
                <div class="chart-container">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
                        <h3>Coverage Breadth (%)</h3>
                        <div>
                            <button onclick="downloadChart('${{sampleName}}', 'breadth', 'png')">📥 Download PNG</button>
                        </div>
                    </div>
                    <div style="position: relative; width: 100%; height: 250px;">
                        <canvas id="chart-breadth-${{sampleName}}"></canvas>
                    </div>
                </div>
            `;
            
            // Create charts
            createChart(`chart-reads-${{sampleName}}`, labels, reads, 'Mapped Reads', '#667eea');
            createChart(`chart-identity-${{sampleName}}`, labels, identities, 'Identity (%)', '#28a745', 100);
            createChart(`chart-depth-${{sampleName}}`, labels, depths, 'Coverage Depth', '#ffc107');
            createChart(`chart-breadth-${{sampleName}}`, labels, breadths, 'Coverage Breadth (%)', '#dc3545', 100);
        }}
        
        // Create single chart
        function createChart(canvasId, labels, data, title, color, maxY = null) {{
            const canvas = document.getElementById(canvasId);
            if (!canvas) return;
            
            new Chart(canvas, {{
                type: 'bar',
                data: {{
                    labels: labels,
                    datasets: [{{
                        data: data,
                        backgroundColor: color,
                        borderColor: color,
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    aspectRatio: 1.5,
                    plugins: {{
                        legend: {{ display: false }}
                    }},
                    scales: {{
                        x: {{
                            ticks: {{
                                maxRotation: 45,
                                minRotation: 45,
                                font: {{ size: 10 }},
                                callback: function(value, index, ticks) {{
                                    const label = this.getLabelForValue(value);
                                    if (label && label.length > 20) {{
                                        return label.substring(0, 17) + '...';
                                    }}
                                    return label;
                                }}
                            }},
                            afterFit: function(scale) {{
                                scale.height = Math.max(scale.height, 80);
                            }}
                        }},
                        y: {{
                            beginAtZero: true,
                            max: maxY,
                            ticks: {{
                                callback: function(value) {{
                                    return maxY === 100 ? value + '%' : value.toLocaleString();
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}
        
        // Update charts when table is sorted
        function updateCharts(sampleName) {{
            if (!charts[sampleName]) return;
            
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            const rows = Array.from(tbody.querySelectorAll('tr:not([style*="display: none"])'));
            
            let labels = [];
            let reads = [];
            let identities = [];
            let depths = [];
            let breadths = [];
            
            rows.forEach(row => {{
                labels.push(row.getAttribute('data-species') || 'N/A');
                reads.push(parseFloat(row.getAttribute('data-reads') || 0));
                identities.push(parseFloat(row.getAttribute('data-identity') || 0));
                depths.push(parseFloat(row.getAttribute('data-depth') || 0));
                breadths.push(parseFloat(row.getAttribute('data-breadth') || 0) * 100);
            }});
            
            const chartReads = Chart.getChart(document.getElementById(`chart-reads-${{sampleName}}`));
            const chartIdentity = Chart.getChart(document.getElementById(`chart-identity-${{sampleName}}`));
            const chartDepth = Chart.getChart(document.getElementById(`chart-depth-${{sampleName}}`));
            const chartBreadth = Chart.getChart(document.getElementById(`chart-breadth-${{sampleName}}`));
            
            if (chartReads) {{
                chartReads.data.labels = labels;
                chartReads.data.datasets[0].data = reads;
                chartReads.update();
            }}
            if (chartIdentity) {{
                chartIdentity.data.labels = labels;
                chartIdentity.data.datasets[0].data = identities;
                chartIdentity.update();
            }}
            if (chartDepth) {{
                chartDepth.data.labels = labels;
                chartDepth.data.datasets[0].data = depths;
                chartDepth.update();
            }}
            if (chartBreadth) {{
                chartBreadth.data.labels = labels;
                chartBreadth.data.datasets[0].data = breadths;
                chartBreadth.update();
            }}
        }}
        
        // Download functions
        function downloadTableAsCSV(sampleName) {{
            const tbody = document.getElementById(`tbody-${{sampleName}}`);
            if (!tbody) return;
            
            let csv = 'Accession,Organism,Viral Species,Segment,Mapped Reads,Identity (%),Coverage Depth,Coverage Breadth (%)\\n';
            const rows = tbody.querySelectorAll('tr:not([style*="display: none"])');
            
            rows.forEach(row => {{
                const cells = row.querySelectorAll('td');
                const rowData = [];
                cells.forEach(cell => {{
                    let text = cell.textContent.trim();
                    if (text.includes(',') || text.includes('"')) {{
                        text = '"' + text.replace(/"/g, '""') + '"';
                    }}
                    rowData.push(text);
                }});
                csv += rowData.join(',') + '\\n';
            }});
            
            const blob = new Blob(['\\ufeff' + csv], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            link.href = URL.createObjectURL(blob);
            link.download = `${{sampleName}}_results.csv`;
            link.click();
        }}
        
        function downloadChart(sampleName, chartType, format) {{
            const canvas = document.getElementById(`chart-${{chartType}}-${{sampleName}}`);
            if (!canvas) return;
            
            const chart = Chart.getChart(canvas);
            if (!chart) return;
            
            if (format === 'png') {{
                // High quality PNG: temporarily increase canvas resolution
                const originalWidth = canvas.width;
                const originalHeight = canvas.height;
                const scale = 3; // 3x resolution for high quality
                
                // Temporarily resize canvas for high-res export
                canvas.width = originalWidth * scale;
                canvas.height = originalHeight * scale;
                chart.resize();
                
                // Export at high resolution
                const url = chart.toBase64Image('image/png', 1.0);
                
                // Restore original canvas size
                canvas.width = originalWidth;
                canvas.height = originalHeight;
                chart.resize();
                
                const link = document.createElement('a');
                link.href = url;
                link.download = `${{sampleName}}_${{chartType}}.png`;
                link.click();
            }}
        }}
        
        // Create heatmap
        function createHeatmap() {{
            const container = document.getElementById('heatmap-container');
            if (!container) return;
            
            const speciesList = Object.keys(heatmapData).sort();
            
            let html = '<table class="heatmap-table"><thead><tr><th>Viral Species</th>';
            allSamples.forEach(sample => {{
                html += `<th>${{sample}}</th>`;
            }});
            html += '</tr></thead><tbody>';
            
            speciesList.forEach(species => {{
                html += `<tr><th>${{species}}</th>`;
                allSamples.forEach(sample => {{
                    const present = heatmapData[species][sample] || false;
                    html += `<td class="heatmap-cell ${{present ? 'heatmap-present' : 'heatmap-absent'}}">${{present ? '✓' : '—'}}</td>`;
                }});
                html += '</tr>';
            }});
            
            html += '</tbody></table>';
            container.innerHTML = html;
        }}
        
        function downloadHeatmapAsPNG() {{
            const container = document.getElementById('heatmap-container');
            if (!container) return;
            
            // High quality PNG: use higher scale
            const scale = 3; // 3x resolution for high quality
            html2canvas(container, {{
                scale: scale,
                useCORS: true,
                logging: false,
                backgroundColor: '#ffffff'
            }}).then(canvas => {{
                const url = canvas.toDataURL('image/png', 1.0);
                const link = document.createElement('a');
                link.href = url;
                link.download = 'heatmap_all_samples.png';
                link.click();
            }});
        }}
        
        // Initialize on page load
        window.addEventListener('DOMContentLoaded', function() {{
            createHeatmap();
            const firstSample = allSamples[0];
            if (firstSample) {{
                const select = document.getElementById('sampleSelect');
                if (select) select.value = firstSample;
                showSample(firstSample);
            }}
        }});
    </script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js"></script>
</body>
</html>
"""
        
        # Write HTML file for this database
        html_file = output_dir / f"results_summary_{database_name}.html"
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Generated HTML visualization for {database_name}: {html_file}")

def find_samples(input_dir):
    """Find all FASTQ samples in input directory (compressed or not)."""
    input_dir = Path(input_dir)
    samples = []
    
    for fastq_file in sorted(input_dir.glob("*.fastq*")):
        if fastq_file.is_file():
            samples.append(fastq_file)
    
    return samples

def main(args=None):
    """Main entry point for the reference selection pipeline."""
    parser = argparse.ArgumentParser(
        description="Select best reference sequence from database for each sample",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Basic usage with default RVDB database:
    virasign -i input_dir -o output_dir
  
  Use specific RVDB version:
    virasign -i input_dir -o output_dir -d RVDB --rvdb-version 31.0
  
  Use both databases with custom accessions:
    virasign -i input_dir -d RVDB,RefSeq -o output_dir -a PX852146.1,NC_123456.1 -t 16
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        dest="input",
        help="Input directory containing FASTQ files"
    )
    
    parser.add_argument(
        "-d", "--database",
        type=str,
        required=False,
        default="RVDB",
        dest="database",
        help="Path to reference database FASTA file, or database name (RVDB, RefSeq, or comma-separated: RVDB,RefSeq). Database names will be automatically downloaded to Databases/ directory. Default: RVDB"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        required=True,
        dest="output",
        help="Output directory for results"
    )
    
    parser.add_argument(
        "--min_identity",
        type=float,
        default=None,  # Will be set based on database type
        help="Minimum average identity percentage (default: 90.0 for RefSeq, 80.0 for other databases)"
    )
    
    parser.add_argument(
        "--min_mapped_reads",
        type=int,
        default=100,
        help="Minimum number of mapped reads (default: 100)"
    )
    
    parser.add_argument(
        "--coverage_depth_threshold",
        type=float,
        default=1.0,
        help="Minimum coverage depth threshold (default: 1.0)"
    )
    
    parser.add_argument(
        "--coverage_breadth_threshold",
        type=float,
        default=0.1,
        help="Minimum coverage breadth threshold (default: 0.1)"
    )
    
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Number of threads to use for minimap2 (default: 1)"
    )
    
    parser.add_argument(
        "-a", "--accession",
        type=str,
        default=None,
        dest="accessions",
        help="NCBI accession number(s) to download and merge with the database. Either: (1) comma-separated list (e.g., -a PX852146.1,NC_123456.1), or (2) text file with one accession per line (e.g., -a accessions.txt). Optional."
    )
    
    parser.add_argument(
        "--enable-clustering",
        action="store_true",
        dest="enable_clustering",
        help="Enable clustering for RVDB database (default: clustering disabled). Use --cluster_identity to set identity threshold."
    )
    
    parser.add_argument(
        "--cluster_identity",
        type=float,
        default=0.98,
        dest="cluster_identity",
        help="Identity threshold for RVDB clustering (default: 0.98, i.e., 98%%). Only used if clustering is enabled with --enable-clustering."
    )
    
    parser.add_argument(
        "--rvdb-version",
        type=str,
        default=None,
        dest="rvdb_version",
        help="RVDB database version to download (e.g., '30.0', '31.0', '29.0'). Default: 31.0. Only applies when using RVDB database. Examples: '31.0' downloads C-RVDBv31.0.fasta.gz"
    )
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    # Set default for enable_clustering (defaults to False if --enable-clustering not specified)
    if not hasattr(args, 'enable_clustering'):
        args.enable_clustering = False
    
    # Set up output directory structure
    # If --output . (current directory), create Virasign_output folder
    # Otherwise use the specified output directory
    # Handle case where current directory was deleted
    try:
        current_dir = Path.cwd()
    except (OSError, FileNotFoundError):
        # Current directory was deleted, use os.getcwd() with error handling
        try:
            current_dir = Path(os.getcwd())
        except (OSError, FileNotFoundError):
            # Fallback to home directory or /tmp
            current_dir = Path(os.environ.get('HOME', '/tmp'))
            print(f"Warning: Current working directory was deleted, using {current_dir} as base", file=sys.stderr)
    
    if args.output == ".":
        output_dir = current_dir / "Virasign_output"
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Clean output directory (remove all old sample directories and files, including log file)
    # Do this BEFORE setup_logging so the log file is overwritten
    if output_dir.exists():
        import shutil
        items_removed = 0
        for item in output_dir.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
                items_removed += 1
            else:
                item.unlink()
                items_removed += 1
        if items_removed > 0:
            # Use print since logger isn't set up yet
            print(f"Cleaned output directory: removed {items_removed} old item(s) from {output_dir}")
    
    # Set up logging (log file goes in output_dir, same as sample outputs)
    # This will create a fresh log file
    setup_logging(output_dir)
    logger.info(f"Cleaned output directory: {output_dir} (removed old sample directories and files)")
    
    # Parse accessions (either comma-separated string OR text file, not both)
    accessions_list = None
    if hasattr(args, 'accessions') and args.accessions:
        acc_arg = args.accessions.strip()
        
        # Check if it's a file path
        if Path(acc_arg).exists() and Path(acc_arg).is_file():
            logger.info(f"Reading accessions from file: {acc_arg}")
            accessions_list = []
            with open(acc_arg, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):  # Skip empty lines and comments
                        accessions_list.append(line)
            logger.info(f"  Found {len(accessions_list)} accession(s) in file")
        else:
            # Treat as comma-separated list
            accessions_list = [a.strip() for a in acc_arg.split(',') if a.strip()]
            logger.info(f"Parsed {len(accessions_list)} accession(s) from comma-separated list")
        
        if accessions_list:
            logger.info(f"Total accessions to download: {len(accessions_list)}")
            for i, acc in enumerate(accessions_list[:5], 1):  # Show first 5
                logger.info(f"  {i}. {acc}")
            if len(accessions_list) > 5:
                logger.info(f"  ... and {len(accessions_list) - 5} more")
    
    # Resolve database path(s) (handles database names like 'RVDB', 'RefSeq', etc.)
    # Can return a single Path or a list of Paths if multiple databases specified
    # Also handles downloading and merging accessions if provided
    try:
        database_result = resolve_database_path(
            args.database, 
            accessions=accessions_list,
            enable_clustering=args.enable_clustering,
            cluster_identity=args.cluster_identity,
            rvdb_version=getattr(args, 'rvdb_version', None)
        )
        if isinstance(database_result, list):
            database_fasta_paths = database_result
            logger.info(f"Resolved databases: {args.database} -> {len(database_fasta_paths)} separate databases")
            if accessions_list:
                logger.info(f"  (merged with {len(accessions_list)} accession(s))")
            for db_path in database_fasta_paths:
                logger.info(f"  - {db_path}")
        else:
            database_fasta_paths = [database_result]
            logger.info(f"Resolved database: {args.database} -> {database_result}")
            if accessions_list:
                logger.info(f"  (merged with {len(accessions_list)} accession(s))")
    except Exception as e:
        logger.error(f"Failed to resolve database '{args.database}': {e}")
        sys.exit(1)
    
    # Auto-detect database type and set default identity threshold if not specified
    # Don't set min_identity globally - set it per database in the loop below
    # This ensures RVDB uses 80% and RefSeq uses 95% when both are used
    user_specified_identity = args.min_identity  # Store user's explicit value if provided
    
    logger.info("="*60)
    logger.info("Virasign: Viral Read ASSIGNment from nanopore sequencing")
    logger.info("="*60)
    logger.info(f"Input directory: {args.input}")
    if len(database_fasta_paths) > 1:
        logger.info(f"Databases: {len(database_fasta_paths)} separate databases (will map against each separately)")
        logger.info(f"  - Identity thresholds will be set per database (95% for RefSeq, 80% for RVDB)")
    else:
        logger.info(f"Database: {database_fasta_paths[0]}")
        # Show which threshold will be used for single database
        db_path_str = str(database_fasta_paths[0]).lower()
        if user_specified_identity is not None:
            logger.info(f"Min identity: {user_specified_identity}% (user-specified)")
        elif "refseq" in db_path_str:
            logger.info(f"Min identity: 95% (RefSeq default)")
        else:
            logger.info(f"Min identity: 80% (default)")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Min mapped reads: {args.min_mapped_reads}")
    logger.info(f"Coverage depth threshold: {args.coverage_depth_threshold}")
    logger.info(f"Coverage breadth threshold: {args.coverage_breadth_threshold}")
    logger.info(f"Threads: {args.threads}")
    if "rvdb" in args.database.lower():
        clustering_status = "enabled" if args.enable_clustering else "disabled"
        if args.enable_clustering:
            logger.info(f"RVDB clustering: {clustering_status} (identity threshold: {int(args.cluster_identity*100)}%%)")
        else:
            logger.info(f"RVDB clustering: {clustering_status}")
    logger.info("="*60)
    
    # Find samples
    input_dir = Path(args.input)
    samples = find_samples(input_dir)
    
    if not samples:
        logger.error(f"No FASTQ files found in {input_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(samples)} sample(s)")
    
    # Process each sample
    for sample_file in samples:
        sample_name = sample_file.stem
        if sample_name.endswith('.fastq'):
            sample_name = sample_file.stem[:-6]  # Remove .fastq
        elif sample_name.endswith('.fq'):
            sample_name = sample_file.stem[:-3]  # Remove .fq
        
        # If multiple databases, process each separately (no combining)
        for database_fasta_path in database_fasta_paths:
            # Detect database type from path
            db_path_str = str(database_fasta_path).lower()
            if "refseq" in db_path_str:
                db_name = "RefSeq"
                is_current_db_refseq = True
            elif "rvdb" in db_path_str:
                db_name = "RVDB"
                is_current_db_refseq = False
            elif "custom" in db_path_str or any(part.endswith('.fasta') and not any(x in part.lower() for x in ['refseq', 'rvdb', 'complete', 'genomes']) for part in db_path_str.split('/')):
                # Custom database (accession-based or merged accessions)
                # Try to infer name from filename or use "Custom"
                db_name = database_fasta_path.stem
                # Clean up common suffixes
                if db_name.endswith('_with_accessions'):
                    db_name = db_name.replace('_with_accessions', '')
                if db_name.endswith('_merged'):
                    db_name = db_name.replace('_merged', '')
                if db_name.endswith('_database'):
                    db_name = db_name.replace('_database', '')
                # If it's still a single accession, use "Custom"
                if len(db_name.split('_')) == 1 and '.' in db_name:
                    db_name = "Custom"
                is_current_db_refseq = False
            else:
                db_name = "Custom"
                is_current_db_refseq = False
            
            # Set database-specific min_identity
            # Use user-specified value if provided, otherwise use database-specific defaults
            if user_specified_identity is not None:
                db_min_identity = user_specified_identity
                logger.info(f"Using user-specified identity threshold: {db_min_identity}% for {db_name}")
            else:
                if is_current_db_refseq:
                    db_min_identity = 95.0  # Higher threshold for RefSeq
                    logger.info(f"Detected {db_name} database - using identity threshold (95%) to reduce false positives")
                else:
                    db_min_identity = 80.0  # Default for RVDB and other databases
                    logger.info(f"Using default identity threshold (80%) for {db_name}")
            
            # Create database-specific subdirectory
            db_output_dir = output_dir / sample_name / db_name if len(database_fasta_paths) > 1 else output_dir / sample_name
            db_output_dir.mkdir(parents=True, exist_ok=True)
            
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing sample: {sample_name} with database: {db_name}")
            logger.info(f"Database file: {database_fasta_path}")
            logger.info(f"Output directory: {db_output_dir}")
            logger.info(f"{'='*60}")
        
            process_sample(
            sample_name,
            sample_file,
                str(database_fasta_path),
                db_output_dir,  # Use database-specific output directory
                min_identity=db_min_identity,  # Use database-specific threshold
            min_mapped_reads=args.min_mapped_reads,
            coverage_depth_threshold=args.coverage_depth_threshold,
                coverage_breadth_threshold=args.coverage_breadth_threshold,
                threads=args.threads
        )
    
    logger.info("\n" + "="*60)
    logger.info("All samples processed successfully!")
    logger.info("="*60)
    
    # Generate HTML visualization of results
    generate_html_visualization(output_dir)

if __name__ == "__main__":
    main()
