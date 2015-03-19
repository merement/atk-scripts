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
#PBS -o Scan_output

# This combines error and output messages in the same file
#PBS -j oe

# Local e-mail will be sent when the job is submitted, starts to run and 
# is aborted
#PBS -m bea

# Explicitly go to the working directory. 
# PLEASE, CHANGE THIS TO YOUR WORKING DIRECTORY, WHERE YOU NEED YOUR JOB
# TO BE EXECUTED. Otherwise, you may run a wrong script or may try to enter 
# a wrong directory, which may result in an error.

# Here we go to the directory where qsub was executed
cd $PBS_O_WORKDIR

# Actual command to be executed. 

atkpython scaneps.py --upmin 1 --downmin 1 --upmax 5 --downmax 5 > outScanSpec.log
