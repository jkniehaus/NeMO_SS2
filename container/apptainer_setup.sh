#!/bin/bash

module add apptainer
# Create sandbox Container
apptainer build --sandbox nemo_env/ docker://ubuntu:22.04

#Create req mount directories for UNC longleaf
mkdir -p nemo_env/nas/longleaf
mkdir -p nemo_env/proj
mkdir -p nemo_env/work
mkdir -p nemo_env/users
mkdir -p nemo_env/overflow
mkdir -p nemo_env/datacommons

#start an interactive container using fakeroot
apptainer shell --writable --fakeroot nemo_env/

#update the system and get your basic tools
apt-get update && apt-get install -y \
        build-essential \
        curl \
        wget \
        git \
        unzip \
        ca-certificates \
        zlib1g-dev \
        xxd \
        cmake \
        libncurses-dev \
        libbz2-dev \
        liblzma-dev \
        software-properties-common \
        locales \
        default-jdk
apt-get update && apt-get install -y locales
# Locale setup to avoid R warnings
locale-gen en_US.UTF-8
update-locale LANG=en_US.UTF-8
#Tool Installation
# STAR
cd /opt
git clone --branch 2.7.11b https://github.com/alexdobin/STAR.git
cd STAR/source
make STAR
cp STAR /usr/local/bin/

# r 4.4.0; this version of ubuntu ships with an earlier R version so we need to install it manually
apt-get install -y software-properties-common
add-apt-repository ppa:marutter/rrutter4.0
apt-get update
apt-get install -y r-base
#get python going
wget https://www.python.org/ftp/python/3.12.2/Python-3.12.2.tgz
tar -xf Python-3.12.2.tgz
cd Python-3.12.2
./configure --enable-optimizations --prefix=/opt/python-3.12.2
make -j$(nproc)
make install
export PATH=/opt/python-3.12.2/bin:$PATH

#Next, install what we use in R:
R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("edgeR")
BiocManager::install("tximport")
install.packages("tidyverse")
quit()

