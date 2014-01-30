#!/usr/bin/Rscript

suppressPackageStartupMessages({
  library("R.utils")
  library("methods")
})

#load cmd arguments:
args <- commandArgs(asValues=T)
if(args$args=="TRUE") stop("please provide path to folder with simulation output")
args$args <- normalizePath(args$args)
setwd(dirname(args$file))

source("postprocessing.R")
postprocessing(args$args[1])