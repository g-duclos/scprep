# Use the official R base image with Ubuntu
FROM rocker/r-ver:4.0.0

# Set maintainer
LABEL maintainer="grant.duclos@gmail.com"
LABEL description="Docker container for scprep R package - single-cell RNA-Seq data preprocessing"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and required R packages
RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install('Biobase')" \
    && R -e "install.packages(c('devtools', 'Matrix', 'remotes'))"

# Install Seurat (handling potential installation issues)
RUN R -e "install.packages('Seurat')"

# Create working directory
WORKDIR /scprep

# Copy package files
COPY . /scprep/

# Install the scprep package
RUN R -e "devtools::install('.', dependencies = TRUE)"

# Set default command to R
CMD ["R"]

# Expose port for RStudio (optional)
EXPOSE 8787