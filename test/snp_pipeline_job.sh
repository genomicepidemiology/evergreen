#!/bin/sh

### Number of nodes
#PBS -l nodes=1:ppn=28
### Memory
#PBS -l mem=120gb
### Requesting time
#PBS -l walltime=48:00:00

# load python2 module
module load anaconda2/4.0.0

# start the pipeline
cd <install directory>
./scripts/parallel_snp_pipeline.py -f <list_of_paths_to_isolates>.iso -b $PWD -p -D -L -E
