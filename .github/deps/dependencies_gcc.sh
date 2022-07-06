#!/usr/bin/env bash

set -eu -o pipefail

sudo add-apt-repository ppa:ubuntu-toolchain-r/test

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
    build-essential csh \
    g++-11 gfortran-11 \
    libopenmpi-dev \
    openmpi-bin \
    libhdf5-openmpi-dev
