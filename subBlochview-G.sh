#!/bin/sh
# A sample of a script to submit an atk job
# Notice that queue "batch" is used. This queue (in perspective) admits as
# many running jobs as the number of available processors. Please, DO NOT 
# submit to this queue ATK related jobs!
# 
#PBS -q batch

# We ask for 4 processors here (see below).
#  #PBS -l nodes=1:ppn=4

#PBS -l walltime=424:00:00

# Output file. In this file the output that's supposed to be on the screen will be written
#PBS -o STG_out

# This combines error and output messages in the same file
#PBS -j oe

# Local e-mail will be sent when the job is submitted, starts to run and 
# is aborted
#PBS -m bea

# Here we go to the directory where qsub was executed
cd $PBS_O_WORKDIR

# Actual command to be executed. 
# Here we specify the k point. Check below that the number of bands is 
# indeed what you need.
kPoint=G

stanalysis.py --bloch --size_up 30 --size_down 30 -k $kPoint > outSt${kPoint}Bloch.log

DFT_file=$(head -n 1 main.txt)
visBloch-red.py -t Bloch_${kPoint}_$DFT_file > outStVis.log