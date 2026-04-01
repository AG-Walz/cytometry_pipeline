FROM rocker/r-ver:4.3.3

LABEL maintainer="AG Walz" \
      description="Cytometry analysis pipeline"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libpng-dev libjpeg-turbo8-dev libcurl4-openssl-dev \
    libssl-dev libhdf5-dev libfontconfig1-dev libfreetype6-dev \
    libfribidi-dev libharfbuzz-dev libtiff5-dev \
    pkg-config make gcc g++ gfortran \
    python3 python3-pip python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Python for leidenalg (used by R leiden via reticulate)
RUN python3 -m pip install --break-system-packages leidenalg igraph numpy

# BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'flowCore', 'flowWorkspace', 'CytoML', 'openCyto', \
    'flowAI', 'PeacoQC', 'FlowSOM', 'flowVS', \
    'preprocessCore', 'flowTrans', 'Biobase', 'BiocGenerics', \
    'BiocParallel', 'cytolib', 'limma', 'diffcyt', \
    'flowDensity' \
), update=FALSE, ask=FALSE)"

# CRAN packages
RUN R -e "install.packages(c( \
    'data.table', 'ggplot2', 'patchwork', 'viridis', 'ggrepel', \
    'pheatmap', 'uwot', 'Rtsne', 'leiden', 'igraph', 'Matrix', \
    'lme4', 'lmerTest', 'emmeans', 'MASS', 'dbscan', 'mclust', \
    'survival', 'survminer', 'sva', 'harmony', \
    'tidyverse', 'future', 'future.apply', 'progress', \
    'ggridges', 'ggpubr', 'scales', 'openxlsx', 'jsonlite', \
    'optparse', 'reshape2', 'diptest', 'ggcorrplot', \
    'RColorBrewer', 'UpSetR', 'VennDiagram', 'eulerr', \
    'nnls', 'devtools' \
), repos='https://cloud.r-project.org')"

# CytoNorm (Bioconductor)
RUN R -e "BiocManager::install('CytoNorm', update=FALSE, ask=FALSE)" || true

# Rphenograph from GitHub
RUN R -e "devtools::install_github('JinmiaoChenLab/Rphenograph', upgrade='never')" || true

# Copy source code into container
COPY src/ /opt/cytometry/src/

ENV CYTOMETRY_SRC=/opt/cytometry/src
ENV RETICULATE_PYTHON=/usr/bin/python3
