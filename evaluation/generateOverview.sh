#!/usr/bin/Rscript

#takes output folder of simulation as input

#load cmd arguments:
args <- commandArgs(TRUE)

source("generateOverview.r")

generateOverview(args[1])

message("overview .pdf file(s) generated")