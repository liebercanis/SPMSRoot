#!bin/bash
# HDF5: C++ interface
export HDF5CPP=/usr/local/Cellar/hdf5/1.12.2
export DYLD_LIBRARY_PATH=$HDF5CPP/lib:$DYLD_LIBRARY_PATH
# HDF2ROOT environment:
export HDF5ROOT=$HOME/legend/SPMSRoot
export DYLD_LIBRARY_PATH=$HDF5ROOT/obj:$DYLD_LIBRARY_PATH
#
printenv  HDF5CPP
printenv  HDF5ROOT
