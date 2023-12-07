#!/bin/bash

# ==================================================================
# Initial setup
# ------------------------------------------------------------------

# Set ENV variables
export APT_INSTALL="apt-get install -y --no-install-recommends"
export PIP_INSTALL="python -m pip --no-cache-dir install --upgrade"
export GIT_CLONE="git clone --depth 10"

# Update apt
sudo apt-get update
sudo apt-get install software-properties-common

# ==================================================================
# Tools
# ------------------------------------------------------------------

DEBIAN_FRONTEND=noninteractive \
sudo $APT_INSTALL \
gcc \
make \
cmake \
gdb \
pkg-config \
apt-transport-https \
build-essential \
apt-utils \
ca-certificates \
wget \
rsync \
git \
vim \
mlocate \
libssl-dev \
curl \
openssh-client \
unzip \
unrar \
zip \
awscli \
csvkit \
emacs \
joe \
jq \
dialog \
man-db \
manpages \
manpages-dev \
manpages-posix \
manpages-posix-dev \
nano \
iputils-ping \
sudo \
ffmpeg \
libsm6 \
libxext6 \
libboost-all-dev


# ==================================================================
# Python
# ------------------------------------------------------------------

#Based on https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa

# # Adding repository for python3.9
# DEBIAN_FRONTEND=noninteractive \
# sudo $APT_INSTALL software-properties-common
# sudo add-apt-repository ppa:deadsnakes/ppa -y
# 
# # Installing python3.9
# DEBIAN_FRONTEND=noninteractive sudo $APT_INSTALL \
# python3.9 \
# python3.9-dev \
# python3.9-venv \
# python3-distutils-extra
# 
# # Add symlink so python and python3 commands use same python3.9 executable
# sudo ln -s /usr/bin/python3.9 /usr/local/bin/python3
# sudo ln -s /usr/bin/python3.9 /usr/local/bin/python
# 
# # Installing pip
# curl -sS https://bootstrap.pypa.io/get-pip.py | python3.9
# export PATH=$PATH:/home/paperspace/.local/bin
# 
# # #
# # # Install other R packages
# # #
# # 
# # conda install --quiet --yes \
#   
# 

# ==================================================================
# Installing CUDA packages (CUDA Toolkit 11.7.1 & CUDNN 8.5.0)
# ------------------------------------------------------------------

# Based on https://developer.nvidia.com/cuda-toolkit-archive
# Based on https://developer.nvidia.com/rdp/cudnn-archive

wget https://developer.download.nvidia.com/compute/cuda/11.7.1/local_installers/cuda_11.7.1_515.65.01_linux.run
sudo sh cuda_11.7.1_515.65.01_linux.run --silent --toolkit
export PATH=$PATH:/usr/local/cuda-11.7/bin
export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64
rm cuda_11.7.1_515.65.01_linux.run

wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/3bf863cc.pub
sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/ /"
sudo $APT_INSTALL libcudnn8=8.5.0.*-1+cuda11.7
sudo $APT_INSTALL libcudnn8-dev=8.5.0.*-1+cuda11.7


# wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
# sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
# wget https://developer.download.nvidia.com/compute/cuda/12.1.1/local_installers/cuda-repo-ubuntu2204-12-1-local_12.1.1-530.30.02-1_amd64.deb
# sudo dpkg -i cuda-repo-ubuntu2204-12-1-local_12.1.1-530.30.02-1_amd64.deb
# sudo cp /var/cuda-repo-ubuntu2204-12-1-local/cuda-*-keyring.gpg /usr/share/keyrings/
# sudo apt-get update
# sudo apt-get -y install cuda


# sudo $GIT_CLONE https://github.com/Kitware/CMake ~/cmake
# cd ~/cmake
# sudo ./bootstrap
# sudo make -j"$(nproc)" install


# ==================================================================
# Node.js and Jupyter Notebook Extensions
# ------------------------------------------------------------------

sudo curl -sL https://deb.nodesource.com/setup_16.x | sudo bash
sudo $APT_INSTALL nodejs
$PIP_INSTALL jupyter_contrib_nbextensions jupyterlab-git
DEBIAN_FRONTEND=noninteractive jupyter contrib nbextension install --user


# ==================================================================
# Config & Cleanup
# ------------------------------------------------------------------

echo "export PATH=${PATH}" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> ~/.bashrc
