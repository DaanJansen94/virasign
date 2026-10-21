FROM condaforge/miniforge3:latest

LABEL maintainer="DaanJansen94" \
      description="Virasign: Viral taxonomic classification for nanopore sequencing data" \
      version="0.1.2" \
      url="https://github.com/DaanJansen94/virasign"

RUN mamba create -y -n virasign -c conda-forge -c bioconda \
        "python>=3.10,<3.13" \
        minimap2=2.24 \
        seqtk=1.3 \
        curl \
        "samtools>=1.17" \
        mmseqs2 \
        nextclade \
        pip \
    && mamba clean -afy

ENV PATH=/opt/conda/envs/virasign/bin:$PATH

COPY . /opt/virasign

RUN pip install --no-deps /opt/virasign

# Do not hard-set ENTRYPOINT to the tool binary.
# Nextflow runs containers via `/bin/bash -c <script>`, and a fixed ENTRYPOINT would turn that into
# `virasign /bin/bash -c ...` and fail. Use CMD so `virasign` remains the default while still
# allowing Nextflow (and users) to override the command.
CMD ["virasign"]
