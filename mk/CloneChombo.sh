#!/usr/bin/env bash

echo "Getting Chambo "
export SUHMO_HOME=${PWD}/..

mkdir build
git clone --branch feature_SUHMO https://github.com/EnnaDelfen/Chombo_3.2.git build/Chombo_3.2
export CHOMBO_HOME=${PWD}/build/Chombo_3.2/lib
cd ${CHOMBO_HOME}/mk
echo "CXX=gcc-11" >> local/Make.def.GH
echo "MPICXX= mpicxx" >> local/Make.def.GH
echo "FC=gfortran-11" >> local/Make.def.GH
echo "USE_HDF=TRUE" >> local/Make.def.GH
echo "HDFINCFLAGS=-I/usr/include/hdf5/openmpi" >> local/Make.def.GH
echo "HDFLIBFLAGS=-L/usr/lib/x86_64-linux-gnu/ -lhdf5_openmpi -lhdf5_openmpi_hl -lz" >> local/Make.def.GH
echo "HDFMPIINCFLAGS=-I/usr/include/hdf5/openmpi" >> local/Make.def.GH
echo "HDFMPILIBFLAGS=-L/usr/lib/x86_64-linux-gnu/ -lhdf5_openmpi -lhdf5_openmpi_hl -lz" >> local/Make.def.GH
rm -rf Make.defs.local
ln -s local/Make.def.GH Make.defs.local
cat Make.defs.local
