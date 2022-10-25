### Pipeline for automated phylogenomic analysis of bacterial isolates (PAPABAC or Evergreen) - COMPARE Data hubs

This pipeline is for the purpose of continuous phylogenetical analysis of bacterial whole-genome sequencing data.  
Isolates are matched and mapped to complete reference genomes (templates). The resulting consensus sequences are the basis of the SNP based phylogenetic trees, that are inferred for each template. These trees could be completely new, if no isolates were matched previously to that subtype, or could contain isolates that were previously added. Therefore ongoing surveillance is performed just by adding new isolates to the system. There is a clustering step during the distance calculation, where isolates with less than 10 SNPs distance are clustered to a 'cluster representative' isolate. These clustered isolates are denoted with an asterisk (\*) in the phylogenetic trees.  

###### Dependencies

Anaconda Python 2.7  
Joblib package 0.13+  
ETE3 package 3.0+  
[KMA 1.3.28](https://bitbucket.org/genomicepidemiology/kma)  
[IQ-tree 1.6](http://www.iqtree.org)  
[Neighbor from the PHYLIP package 3.697](http://evolution.genetics.washington.edu/phylip.html)


###### Installation
```
# Go to wanted location for Evergreen
cd /path/to/install_dir
# Clone COMPARE branch
git clone https://bitbucket.org/genomicepidemiology/evergreen.git -b COMPARE
mv evergreen scripts
```
```
# Install specific kma version
git clone https://bitbucket.org/genomicepidemiology/kma.git -b 1.3.28 --single-branch
cd kma && make
# Add kma folder to PATH
export PATH="${PATH}:${PWD}"
```
```
# Create Anaconda environment
conda env create --file evergreen/scripts/environment.yml
```

###### Folder structure and databases
```
# Create the subdirectories and databases
mkdir logs
mkdir output
mkdir results_db
# ~16Gb
wget https://cge.food.dtu.dk/services/Evergreen/etc/database.tar.gz 
tar -xzf database.tar.gz

# KMA database with default homology reduction settings
mkdir hr_database
mkdir hr_database/current

# Start environment when running the scripts
conda activate evergreen

# run database builder with default settings
$PWD/scripts/build_database_chr.py -r $PWD/scripts/refseq_bacterial_complete_chromosomes_2021.lst -o $PWD/hr_database/current
```

A config file is included for the customisation of temporary directories.

###### Usage

The data __has to persist__ between runs in the base directory, or at least in *results\_db*, as the sqlite databases and consensus sequence files are kept there. Without those, ongoing monitoring is not possible and trees will be only inferred for the most recent isolates.  
The pipeline was designed for multiprocessing and a computer with at least 8 cores are recommended for use. It determines the number of cpu-s and adjusts the number of parallel processes accordingly.  
Maximum likelihood method works only on less than 300 non-redundant isolates. Above that only neighbor-joining trees are inferred.

_Input_  
The input file should only contain new isolates. As long as *results\_db* is available with the files from the previous runs, the new isolates will be processed in addition to the previous isolates.  
*-f* option: tab separated file with two columns, called _isolates file_. The first column is the isolate identifier, the second has the path(s) to the fastq files, comma separated for paired reads.  
Example:
>I000001 /path/to/read_01_1.fastq.gz,/path/to/read_01_2.fastq.gz  
>I000002 /path/to/read_02_1.fastq.gz,/path/to/read_02_2.fastq.gz  

Or *-i* option: file with hashtag followed by path to directory containing raw reads, and the common file names following on new lines, which will be the isolate identifier in the system. Not preferred.
Example:  
>\# /path/to/folder1  
>read_01  
>read_02  
>\# /path/to/folder2  
>read_03  

_Command_  
The *snp_pipeline_job.sh* script shows how the pipeline could be used with Torque queuing system.
```
$ scripts/parallel_snp_pipeline.py -h
usage: parallel_snp_pipeline.py [-h] [-b BASE] [-f ISOLATES_FILE]
                                [-i COLLECTION_FILE] [-a] [-p] [-D] [-L] [-E]
                                [-q]

Parallel SNP pipeline

optional arguments:
  -h, --help          show this help message and exit
  -b BASE             Base (install) directory, absolute path
  -f ISOLATES_FILE    File with sample names and paths to raw reads. Files for
                      same isolate given comma separated.
  -i COLLECTION_FILE  Collection file of raw reads
  -a                  Consider all positions in column for Ns, or
  -p                  Consider the sequences pairwise for Ns
  -D                  Distance based phylogenic tree
  -L                  Maximum likelihood based phylogenic tree
  -E                  ENA compatible output (tsv + tar.gz)
  -q                  Quiet
```
Example of use with pairwise distance calculation method, which is suited for large numbers of isolates of the same subtype, and both neighbor-joining and maximum likelihood method tree inference.
```
/path/to/install_dir/scripts/parallel_snp_pipeline.py -f <list_of_paths_to_isolates>.iso \
-b /path/to/install_dir -p -D -L -E
```

_Output_  
The default output is a list of templates and corresponding newick trees that were inferred in the current run.  
The output with the *-E* option is a tab separated file (.tsv) listing all isolates in the system, together with their templates, phylogenetic tree(s) and distance matrix. The phylogenetic trees are in newick format. The distance matrix is in phylip format. These files are collected in an archive that has the same filename as the tab separated file.


###### Citation

If you use the software in a scientific publication, please cite the following papers:

Szarvas, J.; Ahrenfeldt, J.; Cisneros, J.L.B.; Thomsen, M.C.F.; Aarestrup, F.M.; Lund, O. Large scale automated phylogenomic analysis of bacterial isolates and the Evergreen Online platform. Commun. Biol. 2020, 3, 137, doi:10.1038/s42003-020-0869-5.

Szarvas, J.; Bartels, M.D.; Westh, H.; Lund, O. Rapid Open-Source SNP-Based Clustering Offers an Alternative to Core Genome MLST for Outbreak Tracing in a Hospital Setting. Front. Microbiol.  2021, 12, 703, doi: 10.3389/fmicb.2021.636608
