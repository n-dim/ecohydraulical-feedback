#!/usr/bin/Rscript

#takes output folder of simulation as input

#load cmd arguments:
args <- commandArgs(TRUE)

source("readAllCSV.r")

library("methods")

#set working directory to find files to work on
readAllCSV(args[1])

message("conversion done")