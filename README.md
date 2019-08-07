### Evergreen SNP phylogenetic pipeline

The Evergreen SNP phylogenic pipeline is for the purpose of continuous phylogenetical analysis of bacterial whole-genome sequencing data.  
Isolates are matched and mapped to complete reference genomes (templates). The resulting consensus sequences are the basis of the SNP based phylogenetic trees, that are inferred for each template. These trees could be completely new, if no isolates were matched previously to that subtype, or could contain isolates that were previously added. Therefore ongoing surveillance is performed just by adding new isolates to the system. There is a clustering step during the distance calculation, where isolates with less than 10 SNPs distance are clustered to a 'cluster representative' isolate. These clustered isolates are denoted with an asterisk (\*) in the phylogenetic trees.  
[Preprint on BioRxiv](http://biorxiv.org/cgi/content/short/540138v1)

###### Dependencies

Anaconda Python 2.7  
Joblib package 0.13+  
ETE3 package 3.0+  
[KMA 1.1.7](https://bitbucket.org/genomicepidemiology/kma)  
[IQ-tree 1.6](http://www.iqtree.org)  
[Neighbor from the PHYLIP package 3.697](http://evolution.genetics.washington.edu/phylip.html)


###### Installation
```
# Go to install location
cd /path/to/install_dir
# Clone master branch
git clone https://bitbucket.org/genomicepidemiology/evergreen.git
# Add folder to PATH
export PATH="${PATH}:${PWD}/evergreen/scripts"

# Install specific kma version
git clone https://bitbucket.org/genomicepidemiology/kma.git -b 1.1.7 --single-branch
cd kma && make
# Add kma folder to PATH
export PATH="${PATH}:${PWD}"
```

```
# Create Anaconda environment
conda env create --file evergreen/scripts/environment.yml
```
```
# Create the subdirectories and databases for analysis
cd /path/to/analysis_dir
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
```

###### Usage

The data __has to persist__ between runs in the base directory, or at least in *results\_db*, as the sqlite databases and consensus sequence files are kept there. Without those, ongoing monitoring is not possible and trees will be only inferred for the most recent isolates.  
The pipeline was designed for multiprocessing and a computer with at least 8 cores are recommended for use. It determines the number of cpu-s and adjusts the number of parallel processes accordingly.  
Maximum likelihood method works only on less than 300 non-redundant isolates. Above that only neighbor-joining trees are inferred.

_KMA database and complete genomes folder_  
The classification database should be under hr_database/current, and that and the complete_genomes directory should be in the same analysis directory. (Symlinking the directories from one central place to different analysis folders is possible.)

_Config.py_  
Edit the file to define where the temporary folders should be

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
                                [-i COLLECTION_FILE] [-a] [-p] [-D] [-L] [-q]

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
  -q                  Quiet
```
Example of use with pairwise distance calculation method, which is suited for large numbers of isolates of the same subtype, and both neighbor-joining and maximum likelihood method tree inference.
```
/path/to/install_dir/scripts/parallel_snp_pipeline.py -f <list_of_paths_to_isolates>.iso \
-b /path/to/install_dir -p -D -L
```

_Output_  
The default output is a list of templates and corresponding newick trees that were inferred in the current run.  


###### Test data
Download test isolates from [Ahrenfeldt 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3407-6) and [Timme 2017](https://peerj.com/articles/3893/), and create iso file:
```
cd /path/to/install_dir/evergreen/test
wget ftp://ftp.cbs.dtu.dk/public//CGE/databases/Evergreen/evergreen_test_isolates.tar.gz
tar -xzf evergreen_test_isolates.tar.gz
# Create iso file
./proto_iso.sh
```

Run the test analysis
```
cd /path/to/analysis_dir
conda activate evergreen
parallel_snp_pipeline.py -f /path/to/install_dir/evergreen/test/test_1.iso -b $PWD -p -D -L
parallel_snp_pipeline.py -f /path/to/install_dir/evergreen/test/test_1.iso -b $PWD -p -D -L
conda deactivate
```
###### References
- Ahrenfeldt, J. et al. Bacterial whole genome-based phylogeny: construction of a new benchmarking dataset and assessment of some existing methods. BMC Genomics 18, 19 (2017).
- Clausen, P. T. L. C., Aarestrup, F. M. & Lund, O. Rapid and precise alignment of raw reads against redundant databases with KMA. BMC Bioinformatics 19, 307 (2018).
- Huerta-Cepas, J. et al. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Mol. Biol. Evol. 33, 1635–1638 (2016).
- Nguyen, L.-T., Schmidt, H. A., von Haeseler, A. & Minh, B. Q. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol. Biol. Evol. 32, 268–74 (2015).
- Timme, R. E. et al. Benchmark datasets for phylogenomic pipeline validation, applications for foodborne pathogen surveillance. PeerJ 5, e3893 (2017).

###### License  
Copyright (c) 2017-2019, Judit Szarvas, Technical University of Denmark, unless otherwise noted. All rights reserved.  

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  
See the License for the specific language governing permissions and
limitations under the License.
