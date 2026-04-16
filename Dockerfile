FROM condaforge/miniforge3:latest

LABEL maintainer="DaanJansen94" \
      description="Virasign: Viral taxonomic classification for nanopore data" \
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

ENTRYPOINT ["virasign"]
