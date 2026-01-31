# ============================================================
# BASE: Ubuntu 22.04 with GROMACS compiled for CPU/multithread
# ============================================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------
# System dependencies
# ------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-venv \
    libfftw3-dev \
    libx11-dev \
    libxml2-dev \
    zlib1g-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
# Install GROMACS (CPU build)
# ------------------------------------------------------------
RUN wget https://ftp.gromacs.org/gromacs/gromacs-2022.2.tar.gz && \
    tar xfz gromacs-2022.2.tar.gz && \
    cd gromacs-2022.2 && \
    mkdir build && cd build && \
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_DEFAULT_SUFFIX=OFF -DCMAKE_INSTALL_PREFIX=/gromacs && \
    make -j$(nproc) && make install && \
    echo "source /gromacs/bin/GMXRC" >> /etc/bash.bashrc && \
    rm -rf /gromacs-2022.2.tar.gz gromacs-2022.2

ENV PATH="/gromacs/bin:$PATH"

# ------------------------------------------------------------
# Install Martini dependencies
# ------------------------------------------------------------
RUN pip3 install "martinize2==2.6" MDAnalysis mdtraj numpy

# ------------------------------------------------------------
# Install INSANE
# ------------------------------------------------------------
RUN wget -O /usr/local/bin/insane https://raw.githubusercontent.com/Tsjerk/Insane/master/insane.py && \
    chmod +x /usr/local/bin/insane

# ------------------------------------------------------------
# Install Martini 3 bundle
# ------------------------------------------------------------
RUN mkdir -p /opt/martini && \
    wget -O /opt/martini/martini_v300.zip \
      https://cgmartini-library.s3.ca-central-1.amazonaws.com/martini_v300.zip && \
    cd /opt/martini && \
    unzip martini_v300.zip && \
    rm martini_v300.zip

# ------------------------------------------------------------
# Copy your improved run.sh
# ------------------------------------------------------------
COPY run.sh /usr/local/bin/run.sh
RUN chmod +x /usr/local/bin/run.sh

# ------------------------------------------------------------
# Set entrypoint
# ------------------------------------------------------------
ENTRYPOINT ["/usr/local/bin/run.sh"]

