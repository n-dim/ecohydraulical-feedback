#!/usr/bin/Rscript

#takes output folder of simulation as input

#load cmd arguments:
args <- commandArgs(TRUE)

#set working directory to find files to source
if(!is.na(args[2])) {
  path=args[2]
  setwd(path)
}
source("readAllCSV.r")

library("methods")

readAllCSV(args[1])

message("conversion done")