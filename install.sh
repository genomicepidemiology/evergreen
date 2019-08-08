#!/usr/bin/sh

# Create the subdirectories and databases for analysis
mkdir logs
mkdir output
mkdir results_db
mkdir complete_genomes
# ~7Gb of complete chromosomes from NCBI RefSeq
wget ftp://ftp.cbs.dtu.dk/public//CGE/databases/Evergreen/refseq_complete_chromosomes_151217.tar.gz
tar -xzf refseq_complete_chromosomes_151217.tar.gz -C complete_genomes

# KMA database with default homology reduction settings
mkdir hr_database
mkdir hr_database/current

# Start environment when running the scripts
conda activate evergreen

# run database builder with default settings
build_database_chr.py -r /path/to/install_dir/evergreen/scripts/refseq_bacteria_2017.lst -o $PWD/hr_database/current

# Stop environment when done
conda deactivate
