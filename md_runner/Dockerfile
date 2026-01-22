# syntax=docker/dockerfile:1
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GMX_VERSION=2024.4

# Build deps (CPU build, no GPU)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    git \
    python3 \
    python3-pip \
    libfftw3-dev \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

# Download + build + install GROMACS
RUN wget -q "https://ftp.gromacs.org/gromacs/gromacs-${GMX_VERSION}.tar.gz" \
    && tar -xzf "gromacs-${GMX_VERSION}.tar.gz" \
    && mkdir -p "gromacs-${GMX_VERSION}/build" \
    && cd "gromacs-${GMX_VERSION}/build" \
    && cmake .. \
        -DGMX_BUILD_OWN_FFTW=OFF \
        -DGMX_MPI=OFF \
        -DGMX_OPENMP=ON \
        -DGMX_GPU=OFF \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_BUILD_TYPE=Release \
        -GNinja \
    && ninja -j"$(nproc)" \
    && ninja install \
    && rm -rf /tmp/*

# Make gmx available without "source GMXRC"
ENV PATH="/usr/local/gromacs/bin:${PATH}"

RUN gmx --version | head -n 1
CMD ["bash"]

