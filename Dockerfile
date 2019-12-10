FROM ubuntu:18.04
LABEL Description "https://github.com/jonahcullen/IndexForTheFuture"

ENV DEBIAN_FRONTEND=noninteractive

# Install the necessary packages ontop of base ubuntu installation
RUN apt-get -y update && apt-get install -y \
    curl \
    lsb-release \
    wget \
    git \
    gcc \
    vim \
    build-essential \
    apt-transport-https \
    python3 \
    python3-dev \
    python3-pip \
    s3cmd \
    zlib1g-dev

RUN cd /root

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f
ENV PATH=/root/miniconda3/bin:${PATH}
RUN conda update -n base conda

# Create a python environment
RUN conda config --add channels conda-forge
RUN conda create -y -n default python=3.6
RUN /bin/bash -c "source activate default"

# Install snakemake
RUN conda install -c bioconda -y snakemake=5.8.1

# Install STAR
RUN conda install -c bioconda star=2.6.1

# Install salmon
RUN conda install -c bioconda salmon=0.13.1

# Install minus80 and locuspocus
RUN pip install minus80
RUN pip install locuspocus


COPY . /root/RNAMapping 
RUN cd /root/RNAMapping
WORKDIR /root/RNAMapping
