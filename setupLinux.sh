#!bin/bash
# HDF5: C++ interface
export HDF5CPP=/usr/lib64
export LD_LIBRARY_PATH=$HDF5CPP/lib:$LD_LIBRARY_PATH
# HDF2ROOT environment:
export HDF5ROOT=$HOME/legend/HDF5root
export LD_LIBRARY_PATH=$HDF5ROOT/obj:$HDF5ROOT/lib:$LD_LIBRARY_PATH
#
printenv  HDF5ROOT
printenv  LD_LIBRARY_PATH
