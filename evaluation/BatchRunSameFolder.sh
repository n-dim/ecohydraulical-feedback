#!/bin/bash

# usage: 
# BatchRunSameFolder.sh <inputfile> <simulationNumber>
# output goes into a subfolder with name <simulationNumber> in the same directory as the input file

	inputfile="$1"
	folder=$(dirname $inputfile)
	outfolder=$folder/$2
  mkdir $outfolder
  ecohydModel.out $inputfile $outfolder &&                                   # simulation run
  postprocessing.sh $outfolder
  
  #./readAllCSV.sh $outfolder                                                # convert to .RData
	#./generateOverview.sh $outfolder                                          # generate overview pdf file





######################################

