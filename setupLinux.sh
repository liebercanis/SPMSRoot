#!bin/bash
# HDF5: C++ interface
export HDF5CPP=/usr/lib64
export HDF5ROOT=$HOME/legend/SPMSRoot
# HDF2ROOT environment:
export LD_LIBRARY_PATH=$HDF5CPP/lib:$HDF5ROOT/obj:$LD_LIBRARY_PATH
#
printenv  HDF5ROOT
printenv  LD_LIBRARY_PATH
