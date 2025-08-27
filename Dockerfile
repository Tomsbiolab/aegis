# -----------------------------------------------------------------------------
# STAGE 1: BASE IMAGE & SYSTEM DEPENDENCIES
# -----------------------------------------------------------------------------
FROM ubuntu:20.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install all required system dependencies for building code,
# downloading files, and running tools.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gfortran \
    autoconf \
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    zlib1g-dev \
    wget \
    git \
    unzip \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# STAGE 2: INSTALLING MINICONDA & MAMBA
# -----------------------------------------------------------------------------
# Install Miniconda to manage Python/Bioinformatics environments and packages.
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# Add Conda to the system PATH.
ENV PATH=$CONDA_DIR/bin:$PATH

# Accept Anaconda Terms of Service for the required channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install Mamba for faster package management
RUN conda install -n base -c conda-forge mamba -y

# -----------------------------------------------------------------------------
# STAGE 3: ENVIRONMENT CREATION & SOFTWARE INSTALLATION
# -----------------------------------------------------------------------------
# Install most packages directly with Mamba.
# This is more reliable and avoids compilation/dependency issues.
# Mamba will resolve compatible versions for Python 3.10.
RUN mamba create -n bio_env -c conda-forge -c bioconda -y \
    python=3.10 \
    pip \
    # Bioinformatics tools
    jcvi \
    last \
    diamond \
    minimap2 \
    liftoff \
    orthofinder \
    # Python libraries
    pandas \
    plotly \
    python-kaleido \
    biopython \
    matplotlib \
    scipy \
    colorlover \
    tqdm

# Activate the environment for subsequent commands
SHELL ["conda", "run", "-n", "bio_env", "/bin/bash", "-c"]

# Use pip only for packages not available on Conda or for local installs
ENV PYTHONUNBUFFERED=1

# Install LiftOn (requires numpy, networkx, etc., already installed by Mamba)
RUN git clone https://github.com/Kuanhao-Chao/LiftOn /opt/LiftOn && \
    cd /opt/LiftOn && \
    pip install .

# Install miniprot
RUN git clone https://github.com/lh3/miniprot /opt/miniprot && \
    cd /opt/miniprot && \
    make

ENV PATH="/opt/miniprot:${PATH}"

# -----------------------------------------------------------------------------
# STAGE 4: 'AEGIS' APPLICATION SETUP
# -----------------------------------------------------------------------------
# Switch back to the default shell for copy operations
SHELL ["/bin/bash", "-c"]
RUN git clone https://github.com/Tomsbiolab/aegis.git /aegis
WORKDIR /aegis
ENV PYTHONPATH="/aegis"

# Switch to Conda shell to install the local package
SHELL ["conda", "run", "-n", "bio_env", "/bin/bash", "-c"]
RUN pip install -e .

# -----------------------------------------------------------------------------
# STAGE 5: CONTAINER EXECUTION
# -----------------------------------------------------------------------------
WORKDIR /aegis
ENTRYPOINT ["conda", "run", "-n", "bio_env", "--no-capture-output"]
CMD ["/bin/bash"]
