#!/usr/bin/sh

mkdir logs
mkdir output
mkdir results_db
mkdir scripts
ln -s /home/projects/cge/evergreen/complete_genomes complete_genomes
mkdir hr_database
mkdir hr_database/current
tar -xzvf /home/projects/cge/people/s151038/bacteria_kma_hq99_201711.tar.gz -C hr_database/current
tar -xzvf /home/projects/cge/people/s151038/snp_pipeline_scripts.tar.gz -C scripts
