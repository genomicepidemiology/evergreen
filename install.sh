mkdir logs
mkdir output
mkdir results_db
mkdir scripts
mkdir complete_genomes
wget ftp://ftp.cbs.dtu.dk/public//CGE/databases/Evergreen/complete_genomes_151217.tar.gz
tar -xzvf complete_genomes_151217.tar.gz -C complete_genomes
mkdir hr_database
mkdir hr_database/current
wget ftp://ftp.cbs.dtu.dk/public//CGE/databases/Evergreen/bacteria_kma_hq99_201711.tar.gz
tar -xzvf bacteria_kma_hq99_201711.tar.gz -C hr_database/current
tar -xzvf compare_pipeline_scripts.tar.gz -C scripts