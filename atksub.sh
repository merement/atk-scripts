# sample of a script to submit an atk job
# Notice that atk queue is used. This queue allows for maximum two running jobs to prevent license related errors.
# 
#PBS -q atk

# We ask for 4 processors here. Probably this number can be raised (see below)
#PBS -l nodes=1:ppn=4

# We ask for 212 hours of the computer time. It's okay and even desired to 
# overestimate. Underestimated time will result in forced job termination.
#PBS -l walltime=212:00:00

# Output file. In this file the output that's supposed to be on the screen will be written
#PBS -o DFT_out

# This combines error and output messages in the same file
#PBS -j oe

# Local e-mail will be sent when the job is submitted, starts to run and is aborted
#PBS -m bea

# Explicitly go to the working directory. 
# PLEASE, CHANGE THIS TO YOUR WORKING DIRECTORY, WHERE YOU NEED YOUR JOB
# TO BE EXECUTED. Otherwise, you may run a wrong script or may try to enter 
# a wrong directory, which may result in an error.

cd $PBS_O_WORKDIR

# Prevent threading

export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

# Actual command to be executed. Notice that atkpython is run using mpiexec on 
# the same number of cpus (-n <number os cpus>) as requested above.
# This uses atkpython licenses and we have plenty of them. 
# If, however, the python script invokes DFT calculations (which is usually
# the case, when atkpython is used), this also checks out one atkmaster 
# license and we have only two of them.
#
# Distributing calculations with the help of mpiexec requires atkslave licenses
# available. We have seven such licenses, so in principle the number of cpus
# can be set as high as 8 (1 atkmaster + 7 atkslave = 8). Whether this, indeed,
# results in faster calculations, however, should be checked in each case
# separately. Experimentations showed that it's not impossible that calculations
# may run actually longer! Keep this circumstance in mind.

# Don't forget to pass to atkpython the actual script you want to execute.
# You may also want to change the name of the file whether the output is
# redirected. Sometimes it's helpful to keep track of such outputs.

/opt/mpich2_1.5/bin/mpiexec -n 4 atkpython singleMo6x6.py > outpbs.log

