#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=72

cd $PBS_O_WORKDIR
touch start
LD_LIBRARY_PATH=/home/bayron/lib
export LD_LIBRARY_PATH
./main.x 
touch end


