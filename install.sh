mkdir logs
mkdir output
mkdir results_db
mkdir complete_genomes
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/Evergreen/refseq_complete_chromosomes_151217.tar.gz
tar -xzf refseq_complete_chromosomes_151217.tar.gz -C complete_genomes
mkdir hr_database
mkdir hr_database/current
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/Evergreen/bacteria_kma_hq99_201711.tar.gz
tar -xzf bacteria_kma_hq99_201711.tar.gz -C hr_database/current
