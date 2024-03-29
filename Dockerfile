FROM ubuntu:18.04
LABEL Description "https://github.com/jonahcullen/IndexForTheFuture"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get upgrade --yes

# Install zsh
RUN apt-get install zsh --yes
RUN chsh -s $(which zsh)

# Install the necessary packages ontop of base ubuntu installation
RUN apt-get install -y \
    curl \
    lsb-release \
    wget \
    git \
    gcc \
    vim \
    tmux \
    tree \
    htop \
    build-essential \
#    apt-transport-https \
#    python3 \
#    python3-dev \
#    python3-pip \
    s3cmd \
    zlib1g-dev

RUN mkdir -p /home/.local/{bin,src}
WORKDIR /home/.local/src

# Install oh-my-zsh
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN sh miniconda.sh -b -p /home/.conda
ENV PATH=/home/.conda/bin:${PATH}
RUN conda update -n base conda
RUN /home/.conda/bin/conda init zsh

# Install snakemake, STAR, and salmon
RUN conda install -c bioconda -c conda-forge -y \
    snakemake=6.1.1 \
    star=2.7.8 \
    salmon=1.4.0 \
    mashmap=2.0 \
    gffread=0.12.1 \
    bedtools=2.30.0

# Clone SalmonTools and modify path
RUN git clone https://github.com/COMBINE-lab/SalmonTools.git
ENV PATH=/home/.local/src/SalmonTools/scripts:${PATH}

# Install minus80 and locuspocus
RUN pip install minus80==1.0.0 \
    locuspocus==1.0.2
#RUN pip install -e git+git://github.com/LinkageIO/Minus80.git#egg=minus80
#RUN pip install -e git+git://github.com/LinkageIO/LocusPocus.git#egg=locuspocus

#COPY . /home/RNAMapping 
#WORKDIR /home/RNAMapping
#
#ENTRYPOINT ["/usr/bin/zsh"]
