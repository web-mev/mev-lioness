FROM debian:stretch

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    ap-get install -y build-essential \
    wget \
    git \
    python3 \
    python3-pip

# Install all software under /opt/software:
RUN mkdir -p /opt/software 
RUN mkdir -p /opt/software/resources
ADD panda.py /opt/software/ 
ADD lioness.py /opt/software/
ADD merge_lioness.py /opt/software/
ADD determine_scatter.py /opt/software/

# Get the miniconda installer script and run it in silent mode:
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/software/miniconda.sh && \
  bash /opt/software/miniconda.sh -b -p /opt/software/miniconda
ENV PATH="/opt/software/miniconda/bin:${PATH}"

# Install the NetZoo python library
wget https://github.com/netZoo/netZooPy/archive/refs/tags/0.9.3.tar.gz \
    -P /opt/software
tar xzvf /opt/software/0.8.1.tar.gz -C /opt/software
pip3 install -e /ope/software/netZooPy-0.8.1/
pip3 install joblib
