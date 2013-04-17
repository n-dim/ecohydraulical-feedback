#!/usr/bin/Rscript

#load cmd arguments:
args <- commandArgs(TRUE)

library("methods")
source("CSVtoRData.r")

CSVtoRData(args[1])

message("conversion done")