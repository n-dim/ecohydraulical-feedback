#!/bin/bash


if [ "$1" == "" ]; then
    echo "usage" 
    echo "BatchRun <inputParameterFile> <outputFolder>"
else

  path=$(pwd)/$(dirname $BASH_SOURCE)          # path of this file
  gfortran $path/../model/coupledModel.f90 &&  # compilation
  #./a.out "$1" "$2" &&                         # simulation run
  $path/readAllCSV.sh "$2" "$path"             # convert to .RData

fi
