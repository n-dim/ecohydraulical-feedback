#!/bin/bash


if [ "$1" == "" ]; then
    echo "usage" 
    echo "./BatchRun.sh <inputParameterFile> <outputFolder>"
else
  start_time=$(date +%s)

  gfortran ../model/ecohydModel.f90 -o ../model/ecohydModel.out &&    # compilation
  ../model/ecohydModel.out "$1" "$2" &&                               # simulation run
  ./readAllCSV.sh "$2"                                                # convert to .RData
  ./generateOverview.sh "$2"                                          # generate overview pdf file
  
  finish_time=$(date +%s)
  echo "Calculation time: $((finish_time - start_time)) secs = $(((finish_time - start_time) / 60)) minutes"


fi
