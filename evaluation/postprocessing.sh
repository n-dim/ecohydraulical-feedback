#!/usr/bin/Rscript

#load cmd arguments:
args <- commandArgs(TRUE)

library("methods")
source("postprocessing.R")

postprocessing(args[1])