#!/bin/bash

# usage: 
# parallelBatchRun.sh <folder>
# output goes into a subfolder with name <simulationNumber> in the same directory as the input file

#export SHELL=/bin/sh
ls "$@"/input-sim*.txt -v | parallel --gnu  -j+0 --noswap --joblog joblog.txt "BatchRunSameFolder.sh {} {#} > {//}/log_sim_{#}.txt"


#ls "$@"/input-sim*.txt | parallel --gnu --eta --sshlogin nanu@141.43.4.32,: --trc -j+0 --joblog joblog.txt "./BatchRunSameFolder.sh {} {#} > {//}/log_sim_{#}.txt"
#ls "$@"/input-sim*.txt | parallel --gnu --eta  -j+0 --trc --sshlogin nanu@141.43.4.32,: --joblog joblog.txt "./BatchRunSameFolder.sh {} {#} > {//}/log_sim_{#}.txt"

