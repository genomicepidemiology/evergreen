mkdir logs
mkdir output
mkdir results_db
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/Evergreen/refseq_bacterial_complete_chromosomes_2021.tar.gz
tar -xzf refseq_bacterial_complete_chromosomes_2021.tar.gz
conda activate evergreen
mkdir hr_database
mkdir hr_database/current
$PWD/scripts/build_database_chr.py -r $PWD/scripts/refseq_bacterial_complete_chromosomes_2021.lst -o $PWD/hr_database/current
