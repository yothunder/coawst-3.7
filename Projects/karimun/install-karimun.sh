#!/bin/bash

ml purge
ml restore training-coawst
export NETCDF=/opt/ohpc/pub/apps/netcdf-4.7.3
export WRF_EM_CORE=1
export NETCDF4=1
export HDF5=/opt/ohpc/pub/apps/hdf5-1.10.4
export jasper=/opt/ohpc/pub/apps/jasper-1.900.1
export JASPERLIB=/opt/ohpc/pub/apps/jasper-1.900.1/lib
export JASPERINC=/opt/ohpc/pub/apps/jasper-1.900.1/include
export WRF_DA_CORE=0
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export MCT_INCDIR=${HOME}/COAWST_3.7/MCT_LIB/include
export MCT_LIBDIR=${HOME}/COAWST_3.7/MCT_LIB/lib
ln -sf ../../WRF/run/RRTM* .
ln -sf ../../WRF/run/ozone* .
ln -sf ../../WRF/run/*.TBL .
sh coawst.bash -j 16