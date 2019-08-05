#!/usr/bin/env python2.7

from __future__ import print_function
import sys
import os
from datetime import datetime
from distutils.spawn import find_executable
import subprocess
import shlex
import argparse
try:
    import cPickle as pickle
except ImportError:
    import pickle

parser = argparse.ArgumentParser(
    description='Creates KMA database and the necessary pickles for reference sequences (without plasmids)\n from NCBI for the batch Evergreen run with KMA.')
parser.add_argument(
    '-r',
    dest="ref_list",
    default=None,
    help='Path to the list containing paths to the reference sequences, in priority order')
parser.add_argument(
    '-o',
    dest="out_dir",
    default=None,
    help='Output path for KMA database and pickles')
parser.add_argument(
    '-pre',
    dest="pre",
    default=None,
    help='Optional KMA database prefix, if exists it\'s not re-created')
parser.add_argument(
    '-k',
    dest="ksize",
    type=int,
    default=13,
    help='K-mer size, k-mer prefix not included')
parser.add_argument(
    '-t',
    dest="id",
    type=float,
    default=99.0,
    help='Homology reduction threshold, sequence space, percentage')
parser.add_argument(
    '-sparse',
    dest="sparse",
    default="ATG",
    help='Sparse prefix for downsampling k-mers')
parser.add_argument(
    '-append',
    dest="append",
    action="store_true",
    help='Add to existing KMA database in output directory')
args = parser.parse_args()

def exiting(message):
    print(message, file=sys.stderr)
    sys.exit(1)

def logging(message):
    print("# {}".format(message), file=sys.stderr)

def compress_folder_name(organism_name):
    prevc = ""
    compressed = []
    for c in organism_name:
        if c != "_" or c != prevc:
            compressed.append(c)
        prevc = c
    return ''.join(compressed)

## MAIN
## Check the PATH
if find_executable("kma_index") is None:
    exiting("KMA not in path")

if not os.path.exists(args.ref_list):
    exiting("Input reference list is needed")

if args.out_dir is None or not os.path.exists(args.out_dir):
    exiting("No output path given.")

## Create pickles
"""
NC_013791.2 Bacillus pseudofirmus OF4, complete genome
NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
NC_007164.1 Corynebacterium jeikeium K411 complete genome
NC_002162.1 Ureaplasma parvum serovar 3 str. ATCC 700970, complete genome
NC_004088.1 Yersinia pestis KIM10+, complete genome
NC_002620.2 Chlamydia muridarum Nigg, complete genome
NC_002488.3 Xylella fastidiosa 9a5c, complete genome
NC_002505.1 Vibrio cholerae O1 biovar El Tor str. N16961 chromosome I, complete sequence
NC_002506.1 Vibrio cholerae O1 biovar El Tor str. N16961 chromosome II, complete sequence
NC_002516.2 Pseudomonas aeruginosa PAO1 chromosome, complete genome
"""

tax_db = {}
genome_db = {}
foldername_trans=''.join(chr(c) if chr(c).isalnum() else '_' for c in range(256))

# append option
folder_pic_file = os.path.join(args.out_dir,"bacteria.folder.pic")
GFA_pic_file = os.path.join(args.out_dir,"bacteria.fsa_name.pic")
if args.append:
    if os.path.exists(folder_pic_file) and os.path.getsize(folder_pic_file) > 0:
        with open(folder_pic_file, "rb") as inpic:
            tax_db = pickle.load(inpic)

    if os.path.exists(GFA_pic_file) and os.path.getsize(GFA_pic_file) > 0:
        with open(GFA_pic_file, "rb") as inpic:
            genome_db = pickle.load(inpic)

    logging("Number of references before: {}".format(len(tax_db)))


# Extract fasta headers, create KMA ready input
kma_input_lst = "{}.kma_lst".format(args.ref_list)
try:
    lp = open(args.ref_list, "r")
    op = open(kma_input_lst, "w")
except IOError as e:
    exiting(str(e))

# cmd stumps
tmp_files = []
grep_chk = shlex.split('grep "^>"')
replace_cmd = shlex.split('perl -ne \'if (/^>/ and $. != 1){$_ = "NNNNNNNNNNNN\n"}; print $_;\'')
for line in lp:
    fn = os.path.realpath(line.strip())
    if not line.isspace() and os.path.exists(fn):
        filename = os.path.basename(fn)
        # scaffold chromosomes if necessary
        cmd = grep_chk + [fn]
        try:
            out = subprocess.check_output(cmd)
        except subprocess.CalledProcessError:
            logging("Not found: {}".format(fn))
        else:
            headers = out.strip().split("\n")
            entry_num = len(headers)
            # acc, descr_words
            tmp = headers[0][1:].split()
            if tmp[0] not in tax_db:
                if entry_num == 1:
                    print(fn, file=op)
                else:
                    logging("Multiple entries merged: {}".format(fn))
                    if "{}.kma".format(fn) not in tmp_files:
                        tmp_files.append("{}.kma".format(fn))
                        pmd = replace_cmd + [fn]
                        kma_entry = subprocess.check_output(pmd)
                        with open("{}.kma".format(fn), "w") as outfile:
                            outfile.write(kma_entry)
                        print("{}.kma".format(fn), file=op)


                # Create file name pickle
                genome_db[tmp[0]] = filename

                # Create folder names from fasta headers
                # strip the end of description
                orgn = []
                for w in tmp[1:]:
                    if w in ["complete", "chromosome", "genome", "DNA"]:
                        break
                    else:
                        orgn.append(w.strip())
                if not orgn: #starts with the filter words
                    for w in tmp[1:]:
                        if w not in ["complete", "chromosome", "genome", "sequence", "DNA", "1", "2", "I", "II"]:
                            orgn.append(w)
                # remove comma
                if orgn[-1][-1] == ",":
                    orgn[-1] = orgn[-1][:-1]

                orgname = " ".join(orgn)
                # create folder name by adding accession to the end
                foldername = compress_folder_name("{0}_{1}".format(orgname, tmp[0]).translate(foldername_trans))
                tax_db[tmp[0]] = foldername
            else:
                logging("Double accession skipped: {}".format(tmp[0]))
lp.close()
op.close()

# dump files
with open(folder_pic_file, "w") as outpic:
   pickle.dump(tax_db, outpic)

with open(GFA_pic_file, "w") as outpic:
   pickle.dump(genome_db, outpic)

logging("Pickles dumped to {}".format(args.out_dir))

# Create KMA index
# calculate k-mer threshold, ie 90% ~ 18.53, 85% ~ 7.43
kmer_id =  round(((args.id/100)**args.ksize)*100, 2)

## Create KMA database for -Sparse method, homology reducing the sequences
kma_output = ""
if args.pre is not None:
    kma_output = os.path.join(args.out_dir, args.pre)
else:
    kma_output = os.path.join(args.out_dir, "{}_k{}_hr{}_{}".format(os.path.basename(args.ref_list).rsplit(".", 1)[0], args.ksize, args.id, args.sparse))

# KMA index is in decimals
hr_thr = kmer_id
if hr_thr > 1.0:
    hr_thr = hr_thr / 100
hr_part = " -ht {0:.4f} -hq {0:.4f} -and".format(hr_thr)

logging("KMA options:")
logging("k-mer size: {}".format(args.ksize))
logging("Sparse prefix: {}".format(args.sparse))
logging("HR threshold: {} ({})".format(args.id, kmer_id))
logging("Input batch file: {}".format(args.ref_list))
logging("Output path: {}".format(kma_output))

# if not present or zero sized
kma_log_file = ""
if not os.path.exists("{}.comp.b".format(kma_output)) or not os.path.getsize("{}.comp.b".format(kma_output)) > 0:
    index_cmd = "kma_index -batch {0} -o {1} -k {2} -NI -Sparse {3}".format(kma_input_lst, kma_output, args.ksize, args.sparse)
    index_cmd += hr_part

    kma_log_file = "{}.log".format(kma_output)
    with open(kma_log_file, "w") as ofile:
        p = subprocess.call(shlex.split(index_cmd), stdout=ofile, stderr=subprocess.STDOUT)

    logging("KMA indexing done")
# if adding to existing database
elif args.append:
    index_cmd = "kma_index -batch {0} -t_db {1} -k {2} -NI -Sparse {3}".format(kma_input_lst, kma_output, args.ksize, args.sparse)
    index_cmd += hr_part

    kma_log_file = "{}.{}.log".format(kma_output, datetime.today().strftime("%Y-%m-%d"))
    with open(kma_log_file, "w") as ofile:
        p = subprocess.call(shlex.split(index_cmd), stdout=ofile, stderr=subprocess.STDOUT)

    logging("KMA database appended")
else:
    logging("KMA database existing")

# clean up
for tmpfile in tmp_files:
    try:
        os.unlink(tmpfile)
    except OSError as e:
        logging(str(e))

if not os.path.exists("{}.comp.b".format(kma_output)) or not os.path.getsize("{}.comp.b".format(kma_output)) > 0:
    exiting("KMA indexing unsuccessful")
