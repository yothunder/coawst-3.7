#!/bin/bash

#SBATCH -J C-SANDY
#SBATCH -o c-sandy.log
#SBATCH --nodes 1
#SBATCH --ntasks 3
#SBATCH --time 01:00:00

ml purge
ml restore training-coawst
export WRF_EM_CORE=1
export NETCDF=/opt/ohpc/pub/apps/netcdf-4.7.3/
export NETCDF4=1
export HDF5=/opt/ohpc/pub/apps/hdf5-1.10.4/
export jasper=/opt/ohpc/pub/apps/jasper-1.900.1/
export JASPERLIB=${jasper}/lib
export JASPERINC=${jasper}/include
export WRF_DA_CORE=0
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
export MCT_INCDIR=/home/tr16/COAWST_3.7/MCT_LIB/include
export MCT_LIBDIR=/home/tr16/COAWST_3.7/MCT_LIB/lib

cd /home/tr16/COAWST_3.7/Projects/Sandy
ln -sf /home/tr16/COAWST_3.7/*TBL .
ln -sf /home/tr16/COAWST_3.7/RRTM* .
mpirun -np 3 ./coawstM coupling_sandy.in 
