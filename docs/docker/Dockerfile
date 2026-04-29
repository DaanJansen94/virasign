FROM condaforge/miniforge3:latest

LABEL maintainer="DaanJansen94" \
      description="Virasign: Viral taxonomic classification for nanopore sequencing data" \
      version="0.0.7" \
      url="https://github.com/DaanJansen94/virasign"

COPY . /opt/virasign

RUN mamba install -y -c conda-forge -c bioconda \
        "python>=3.9,<3.14" \
        minimap2=2.24 \
        seqtk=1.3 \
        curl \
        "samtools>=1.17" \
        mmseqs2 \
        nextclade \
    && mamba clean -afy

RUN pip install --no-deps /opt/virasign

# Do not hard-set ENTRYPOINT to the tool binary.
# Nextflow runs containers via `/bin/bash -c <script>`, and a fixed ENTRYPOINT would turn that into
# `virasign /bin/bash -c ...` and fail. Use CMD so `virasign` remains the default while still
# allowing Nextflow (and users) to override the command.
CMD ["virasign"]
