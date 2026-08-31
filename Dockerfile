# Use NVIDIA CUDA base image
FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

# Set non-interactive mode for apt-get to avoid prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install necessary packages
RUN apt-get update && \
    apt-get clean && \
    apt-get install -y \
    wget \
    python3-pip \
    git \
    libeigen3-dev \
    libboost-all-dev \
    gcc \
    g++ \
    cmake 

# Install CUDA 12.x
RUN echo "Checking and installing CUDA 12.x..."
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
RUN dpkg -i cuda-keyring_1.1-1_all.deb
RUN apt-get update -y
RUN apt-get -y install cuda-toolkit-12-4
ENV PATH=/usr/local/cuda/bin:$PATH
ENV CUDADIR=/usr/local/cuda
ENV CXXFLAGS="-I/usr/local/cuda/include $CXXFLAGS"
ENV LDFLAGS="-L/usr/local/cuda/lib64 $LDFLAGS"

# Set the working directory
WORKDIR /app

# Copy the repository from the build context to the container
RUN mkdir -p /app/stormm
COPY . /app/stormm

# Build STORMM using CMake with CUDA enabled.  The architecture list must be quoted so the shell
# does not treat the semicolon as a command separator.  80 serves data-center Ampere (A100, A30)
# and 89 serves Lovelace (L40S, RTX 40xx).  Each entry emits both SASS and PTX, so the image also
# retains JIT forward-compatibility; adding architectures increases build time and image size
# roughly in proportion to the length of the list.
#
# CUSTOM_NVCC_THREADS lets nvcc run the two per-architecture device passes of a translation unit
# concurrently, which recovers most of the wall-clock cost of the second architecture.
RUN cmake -S stormm -B stormmbuild \
    -DSTORMM_ENABLE_CUDA=YES \
    -DSTORMM_ENABLE_RDKIT=NO \
    -DCUSTOM_GPU_ARCH="80;89" \
    -DCUSTOM_NVCC_THREADS=2

ENV STORMM_HOME=/app/stormm
ENV STORMM_SOURCE=/app/stormm
ENV STORMM_BUILD=/app/stormmbuild
WORKDIR /app/stormmbuild
# Parallelism is capped deliberately rather than left as a bare "make -j".  STORMM's CUDA
# translation units are individually memory-hungry, and Cloud Build's largest available
# --machine-type (e2-highcpu-32) provides only 32 GB of RAM, so an unbounded job count is liable to
# be OOM-killed.  8 make jobs times 2 nvcc threads is roughly 16 concurrent device passes, which
# fits the available memory while still using a fair share of the 32 vCPUs.  Raise this only in
# step with the memory of the machine actually performing the build.
RUN make -j8

# After building all apps
COPY entrypoint /usr/local/bin/entrypoint
RUN chmod +x /usr/local/bin/entrypoint

# Set the default command to run bash
RUN echo "To run this container with GPU support, use the --gpus flag with docker run (e.g. docker run --gpus all stormm-config)"
ENTRYPOINT ["/usr/local/bin/entrypoint"]
