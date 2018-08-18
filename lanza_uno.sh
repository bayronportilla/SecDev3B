#!/bin/bash 

#PBS -l walltime=00:01:00
cd $PBS_O_WORKDIR
LD_LIBRARY_PATH=/home/bayron/lib
export LD_LIBRARY_PATH
./main.x 


