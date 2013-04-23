#!/bin/bash


if [ "$1" == "" ]; then
    echo "usage" 
    echo "./BatchRun.sh <inputParameterFile> <outputFolder>"
else

  gfortran ../model/ecohydModel.f90 -o ../model/ecohydModel.out &&    # compilation
  ../model/ecohydModel.out "$1" "$2" &&                               # simulation run
  ./readAllCSV.sh "$2"                                                # convert to .RData

fi
