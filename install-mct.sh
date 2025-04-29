#!/bin/bash

ml purge 
ml restore training-coawst
cd ~/COAWST_3.7
mkdir -p MCT_LIB
cd Lib/MCT

export NETCDF=/opt/ohpc/pub/apps/netcdf-4.7.3
# ./configure --prefix=/home/tr16/COAWST_3.8/MCT_LIB F90=gfortran FC=mpif90 FCFLAGS=-fallow-argument-mismatch # gnu9 harus ada fallow argument mismatch
./configure --prefix=/home/tr16/COAWST_3.7/MCT_LIB F90=gfortran FC=mpif90
make
make install

#cd /home/tr16/COAWST_3.7/Lib/SCRIP_COAWST
#make
#make install
