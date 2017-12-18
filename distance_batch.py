#!/usr/bin/env python2.7
from __future__ import print_function
import sys, time
import os
import math
import numpy as np
import array
import argparse
from operator import itemgetter
import gzip
import pickle
from joblib import Parallel, delayed
import sqlite3

# quick hack
base_path = os.path.dirname(sys.argv[0]).rsplit("/",1)[0]
MAIN_SQL_DB = os.path.join(base_path, "results_db/evergreen.db")
NO_JOBS = 8
THRESHOLD = 10

############# ITERATORS #############
def SeqsFromFile(filename, exit_on_err=False):
   '''Extract sequences from a file

   Name:
      SeqsFromFile
   Author(s):
      Martin CF Thomsen
   Date:
      18 Jul 2013
   Description:
      Iterator which extract sequence data from the input file
   Args:
      filename: string which contain a path to the input file
   Supported Formats:
      fasta, fastq

   USAGE:
   >>> import os, sys
   >>> # Create fasta test file
   >>> file_content = ('>head1 desc1\nthis_is_seq_1\n>head2 desc2\n'
                       'this_is_seq_2\n>head3 desc3\nthis_is_seq_3\n')
   >>> with open('test.fsa', 'w') as f: f.write(file_content)
   >>> # Parse and print the fasta file
   >>> for seq, name, desc in SeqsFromFile('test.fsa'):
   ...    print ">%s %s\n%s"%(name, desc, seq)
   ...
   >head1 desc1
   this_is_seq_1
   >head2 desc2
   this_is_seq_2
   >head3 desc3
   this_is_seq_3
   '''
   # VALIDATE INPUT
   if not isinstance(filename, str):
      msg = 'Filename has to be a string.'
      if exit_on_err: sys.exit('Error: '+msg)
      else: raise IOError, msg
   if not os.path.exists(filename):
      msg = 'File "%s" does not exist.'%filename
      if exit_on_err: sys.exit('Error: '+msg)
      else: raise IOError, msg

   # EXTRACT DATA
   with open_(filename,"r") as f:
      queryseqsegments = []
      seq, name, desc = '', '', ''
      line = ''
      nextline = f.next
      addsegment = queryseqsegments.append
      for line in f:
         if len(line.strip()) == 0: continue
         #print("%s\n"%line, file=sys.stderr)
         fields=line.strip().split()
         if line[0] == ">":
            # FASTA HEADER FOUND
            if queryseqsegments != []:
               # YIELD SEQUENCE AND RESET
               seq = ''.join(queryseqsegments)
               yield (seq, name, desc)
               seq, name, desc = '', '', ''
               del queryseqsegments[:]
            name = fields[0][1:]
            desc = ' '.join(fields[1:])

         elif line[0] == "@":
            # FASTQ HEADER FOUND
            name = fields[0][1:]
            desc = ' '.join(fields[1:])
            try:
               # EXTRACT FASTQ SEQUENCE
               line = nextline()
               seq  = line.strip().split()[0]
               # SKIP SECOND HEADER LINE AND QUALITY SCORES
               line = nextline()
               line = nextline() # Qualities
            except:
               break
            else:
               # YIELD SEQUENCE AND RESET
               yield (seq, name, desc)
               seq, name, desc = '', '', ''

         elif len(fields[0])>0:
            # EXTRACT FASTA SEQUENCE
            addsegment(fields[0])

      # CHECK FOR LAST FASTA SEQUENCE
      if queryseqsegments != []:
         # YIELD SEQUENCE
         seq = ''.join(queryseqsegments)
         yield (seq, name, desc)


############# FUNCTIONS #############
def open_(filename, mode=None, compresslevel=9):
   """Switch for both open() and gzip.open().

   Determines if the file is normal or gzipped by looking at the file
   extension.

   The filename argument is required; mode defaults to 'rb' for gzip and 'r'
   for normal and compresslevel defaults to 9 for gzip.

   """
   if filename[-3:] == '.gz':
      if mode is None: mode = 'rb'
      return gzip.open(filename, mode, compresslevel)
   else:
      if mode is None: mode = 'r'
      return open(filename, mode)

def timing(message):
    if not args.quiet:
        t1 = time.time()
        print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)

def encode_seq_pw(strain):
   encodedinput = [np.zeros((clen), dtype=np.int8) for clen in clens]
   non_nuc_mask = [np.zeros((clen), dtype=np.int8) for clen in clens]
   for j in xrange(len(clens)):
      for i in xrange(clens[j]):
         try:
            encodedinput[j][i] = nuc2num[strain[j][i]]
         except KeyError:
            non_nuc_mask[j][i] = 1
   return encodedinput, non_nuc_mask

def encode_seq_all(strain):
   encodedinput = [np.zeros((clen), dtype=np.int8) for clen in clens]
   non_nuc_mask = [np.ones((clen), dtype=np.bool) for clen in clens]
   for j in xrange(len(clens)):
      for i in xrange(clens[j]):
         try:
            encodedinput[j][i] = nuc2num[strain[j][i]]
         except KeyError:
            non_nuc_mask[j][i] = False
   return encodedinput, non_nuc_mask

# def dist_calc_all(mat1i, iso):
#    return (np.not_equal(mat1i, iso * np.ones((len(mat1i), len(iso)), np.int8))).sum(1)

def dist_calc_all(mat1i, mat2i, i, j):
   return (np.not_equal(mat1i, mat2i[j] * np.ones((slens[0], clens[i]), np.int8))).sum(1)

def dist_calc_all_old(mat1i, i, j):
   return (np.not_equal(mat1i, mat1i[j] * np.ones((slens[0], clens[i]), np.int8))).sum(1)

def dist_calc_pw(inputhrmati, inputnwmati, maskhri, maskni, i, j):
   # distance between _new_ strain and all old
   #ones = np.ones((slens[0], clens[i]), dtype=np.int8)
   dist1 = (np.not_equal(inputhrmati, inputnwmati[j] * np.ones((slens[0], clens[i]), dtype=np.int8))).sum(1)
   # number of non-corresponding N-s
   dist2 = (np.not_equal(maskhri, maskni[j] * np.ones((slens[0], clens[i]), dtype=np.int8))).sum(1)
   return np.array(dist1 - dist2)

def dist_calc_pw_old(inputhrmati, maski, i, j):
   # distance between one _old_ strain and all old !!!!
   # ones = np.ones((slens[0], clens[i]), dtype=np.int8)
   dist1 = (np.not_equal(inputhrmati, inputhrmati[j] * np.ones((slens[0], clens[i]), dtype=np.int8))).sum(1)
   # number of non-corresponding N-s
   dist2 = (np.not_equal(maski, maski[j] * np.ones((slens[0], clens[i]), dtype=np.int8))).sum(1)
   return np.array(dist1 - dist2)

def remain_len(i, j, ones, mask):
   # calculate the remaining length for the two seqs for this chromosome
   # sum([clens[chr] - np.logical_or( arr1[chr], arr2[chr]).sum(1) for chr in [0,1]]) for all chroms
   return (clens[i] - np.logical_or(mask[0][i], mask[1][i][j] * ones).sum(1))

def remain_len_old(i, j, ones, mask):
   # calculate the remaining length for the two seqs for this chromosome
   # sum([clens[chr] - np.logical_or( arr1[chr], arr2[chr]).sum(1) for chr in [0,1]]) for all chroms
   return (clens[i] - np.logical_or(mask[i], mask[i][j] * ones).sum(1))


def n_content(sequences):
   """Test the sequence for the number of N bases."""
   t_len = 0
   n_count = 0.0
   for chrom in sequences:
      t_len += len(chrom)
      n_count += chrom.count("N")
   if ((n_count / t_len) > 0.100):
      return True
   return False

def seq_to_homol(cluster):
    """Convert cluster dict to dict of lists."""
    redict = {} # homol : [seqs]
    for value in set(cluster.values()):
        redict[value] = [key for key in cluster if cluster.get(key) == value]
    return redict

def add_clusters(ref, isolates):
    # clusters = { hr-name0 : [ (isolate0, distance), (isolate1, distance), ... ], ... }
    if ref[0] == 0:
        if ref[1] == 0: # template
            homolname = "template"
        else:
            homolname = oldseqs[ref[1]].split(".")[0]
    else:
        homolname = newseqs[ref[1]].split(".")[0]
    if clusters.get(homolname) is None:
        clusters[homolname] = []
    for ind in isolates:
        isoname = newseqs[ind].split(".")[0]
        if ref[0] == 0:
            dist = dist_old[(ref[1], ind)]
        else:
            dist = dist_new[(ref[1], ind)]
        if (isoname,dist) not in clusters[homolname]: # we dont want to re-add it
            clusters[homolname].append((isoname, dist))
        # prepare db insert
        cluster_insert.append((isoname, homolname, dist))
    return

def create_vcf(ref, seqs):
    """Create VCF file from redundant sequences."""
    #ref: path, vcfname: seqname from path

    date = time.strftime("%Y%m%d")
    version = "VCFv4.2"
    source = sys.argv[0].split("/")[-1]
    if ref[0] == 0:
        reference = inputhrseq[ref[1]]
        reffilename = oldseqs[ref[1]]
    else:
        reference = inputnewseq[ref[1]]
        reffilename = newseqs[ref[1]]
    inputseq = [inputnewseq[x] for x in seqs]
    inputseqname = [newseqs[x].split(".")[0] for x in seqs]
    nchrom = len(reference)
    clens = [len(chrom) for chrom in reference]
    nseq = len(inputseq)

    vcffiles = []
    for i, seq in enumerate(inputseq):
        headertext = []
        headertext.append("##fileformat={0}".format(version))
        headertext.append("##source={0}".format(source))
        headertext.append("##reference={0}".format(os.path.join(args.odir, reffilename)))
        for i, chrom in enumerate(seq):
            headertext.append("##contig=<ID={0},length={1}>".format(i+1, clens[i]))
        headertext.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        vcffiles.append(headertext)

    for j in range(nchrom):
        for k in xrange(clens[j]):
            for i in xrange(0,nseq):
                if reference[j][k] != inputseq[i][j][k]:
                    vcffiles[i].append("{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.".format(j, k, reference[j][k], inputseq[i][j][k]))

    for i, vcfname in enumerate(inputseqname):
        with open(os.path.join(args.odir, vcfname + ".vcf"), "w") as vcffile:
            print("\n".join(vcffiles[i]), file=vcffile)

    return

# Start time to keep track of progress
t0 = time.time()
etta = 0.001

# Parse command line options
parser = argparse.ArgumentParser(
    description='Calculates genetic distance between two set of input sequences, reduces homology')
parser.add_argument("-hr", dest="old", help="read non-homologous sequences from", metavar="INFILE_HR")
parser.add_argument("-n", dest="new", help="read new sequences from", metavar="INFILE_NEW")
parser.add_argument("-o", "--outputdir", dest="odir", help="write to DIR", metavar="DIR")
parser.add_argument("-m", "--outputmat", dest="outputfilenamemat", help="write to OUTMAT", metavar="OUTMAT")
parser.add_argument("-a", "--allcalled", dest="allcalled", action="store_true", help="Only use positions called in all strains")
parser.add_argument("-v", "--vcf", dest="vcf", action="store_true", help="Create VCF files from homologous sequences")
parser.add_argument("-d", "--debug", dest="debug", action="store_true", help="Debug: use .t suffix")
parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Quiet")
args = parser.parse_args()

suffix = ""
if args.debug:
    suffix = ".t"

# File with reads
if args.old is None or args.new is None:
    exiting('No input filelist was provided')
elif not os.path.exists(args.old) or not os.path.exists(args.new):
    exiting('The referred input filelist does not exists')

# path to results dir ../results_db/template/
if args.odir is None:
    exiting('Output directory is needed')

# File for distance matrix output
if args.outputfilenamemat is None:
    if args.allcalled:
        outputmat = os.path.join(args.odir, "dist.all.mat{0}".format(suffix))
    else:
        outputmat = os.path.join(args.odir, "dist.pw.mat{0}".format(suffix))
else:
    # it should be in the output dir
    if os.path.split(args.outputfilenamemat)[0] == args.odir:
        outputmat = args.outputfilenamemat
    else:
        outputmat = os.path.join(args.odir, os.path.split(args.outputfilenamemat)[1])


# New strains
inputnewseq = []
inputnewname = []
inputnewdesc = []
newseqs = []

templ_update = []
runs_update = []
templ = os.path.basename(args.odir)

with open(args.new) as f:
   for l in f:
      l = l.strip()
      if l == "":
        next #avoid empty line
      fp = os.path.join(args.odir, l)
      if os.path.exists(fp) and os.path.getsize(fp) > 0:
         entries = zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)])
         if n_content(entries[0]):
            print("# {0} entry has too many N calls".format(entries[1][0]), file=sys.stderr)
            templ_update.append((0, l.split(".")[0], templ))
            runs_update.append((2, l.split(".")[0]))
            continue
         newseqs.append(l)
         inputnewseq.append(list(entries[0])) #tuple with number of chromosomes
         #inputnewname.append(list(entries[1]))
         #inputnewdesc.append(list(entries[2]))

if templ_update:
    # reach back and update runs table in the main DB to 2 as non-included
    conn = sqlite3.connect(MAIN_SQL_DB)
    conn.execute("PRAGMA foreign_keys = 1")
    conn.commit()
    cur = conn.cursor()
    try:
        cur.executemany('''UPDATE templates SET qc_pass=? WHERE sra_id=? and template=?''', templ_update)
        cur.executemany('''UPDATE runs SET included=? WHERE sra_id=?''', runs_update)
        conn.commit()
        conn.close()
    except sqlite3.Error:
        print("Warning: SQL update failed.", file=sys.stderr)

# Exit if no knew
if len(newseqs) == 0:
    exiting("No new sequences fullfilled the N-content criterium.")

# homology reduced old isolates
inputhrseq = []
inputhrname = []
inputhrdesc = []
oldseqs = [] # hr reduced seqs
# TODO get these from db, where repr_id is NULL
# change so oldseqs, newseqs is just the id, not the filename!

with open(args.old) as f:
   for l in f:
      # First strain should be the the template
      l = l.strip()
      if l == "":
         next #avoid empty line
      fp = os.path.join(args.odir, l)
      if os.path.exists(fp):
         oldseqs.append(l)
         entries = zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)])
         inputhrseq.append(list(entries[0]))
         #inputhrname.append(list(entries[1]))
         #inputhrdesc.append(list(entries[2]))

timing("# Read inputfiles from the file lists.")
# print(oldseqs)
# print(newseqs)

# open database
if args.allcalled:
    db_path = os.path.join(args.odir, "isolates.all.db{}".format(suffix))
else:
    db_path = os.path.join(args.odir, "isolates.pw.db{}".format(suffix))
conn = sqlite3.connect(db_path)
conn.execute("PRAGMA foreign_keys = 1")
cur = conn.cursor()

# Define nucleotides as numbers
nuc2num = {
   # A  adenosine       C  cytidine          G  guanine
   # T  thymidine       N  A/G/C/T (any)     U  uridine
   # K  G/T (keto)      S  G/C (strong)      Y  T/C (pyrimidine)
   # M  A/C (amino)     W  A/T (weak)        R  G/A (purine)
   # B  G/T/C           D  G/A/T             H  A/C/T
   # V  G/C/A           -  gap of indeterminate length
   #   A C G T
   # A A M R W
   # C M C S Y
   # G R S G K
   # T W Y K T
   'A' : 1,
   'T' : 2,
   'C' : 3,
   'G' : 4,
   'M' : 5,
   'R' : 6,
   'W' : 7,
   'S' : 8,
   'Y' : 9,
   'K' : 10
}

clens = [len(chrom) for chrom in inputhrseq[0]] # Length of chromosomes
tot_len = sum(clens) #total length of sequences, for the pw calc
nchrom = len(clens)
slens = [len(inputhrseq), len(inputnewseq)] # Number of strains in each list
if not args.quiet:
    print("# Length of chromosomes %s" %(clens), file=sys.stdout)
    print("# Number of strains: %s" % (slens), file=sys.stdout)
inputhrseqmat = []
inputnewseqmat = []
non_nuc_masks = []
mask_nw = []
mask_hr = []


if args.allcalled:
    for inputseqmat, inputseq, mask in ([inputhrseqmat, inputhrseq, mask_hr], [inputnewseqmat, inputnewseq, mask_nw]):

        arrays = Parallel(n_jobs=NO_JOBS * 2)(delayed(encode_seq_all)(isolate) for isolate in inputseq)

        inputseqmat.extend([np.asarray([item[0][i] for item in arrays]) for i in range(nchrom)])
        mask.extend([np.asarray([item[1][i] for item in arrays]).all(axis=0) for i in range(nchrom)]) # take AND for mask matrix by evaling if all True
    non_nuc_masks = [np.logical_and(np.asarray(mask_hr[i]), np.asarray(mask_nw[i])) for i in range(nchrom)] # take AND for the hr and nw chrom arrays

else:
    for inputseqmat, inputseq in ([inputhrseqmat, inputhrseq], [inputnewseqmat, inputnewseq]):

        arrays = Parallel(n_jobs=NO_JOBS * 2)(delayed(encode_seq_pw)(isolate) for isolate in inputseq)

        inputseqmat.extend([np.asarray([item[0][i] for item in arrays]) for i in range(nchrom)])
        non_nuc_masks.append([np.asarray([item[1][i] for item in arrays]) for i in range(nchrom)])

# free some memory
if not args.vcf:
    del inputhrseq
    del inputnewseq
    del arrays

timing("# Encoded sequences into numpy array.")

# Calculate pairwise distances
if not args.allcalled:
    if slens[0] > 300 or slens[1] > 200:
        NO_JOBS = int(NO_JOBS / max((float(slens[0])/300),(float(slens[1])/200)))
else:
    if slens[0] > 500 or slens[1] > 200:
        NO_JOBS = int(NO_JOBS / max((float(slens[0])/500),(float(slens[1])/200)))


if args.allcalled:
    # Mask non-nucleotide positions
    mat1 = [inputhrseqmat[i].T[non_nuc_masks[i]].T for i in xrange(nchrom)]
    mat2 = [inputnewseqmat[i].T[non_nuc_masks[i]].T for i in xrange(nchrom)]
    clens = [len(chromosome[0]) for chromosome in mat1]
    if not args.quiet:
        print("# Old seqs reduced matrix shape:", [np.shape(mat1[i]) for i in xrange(nchrom)]) #(1, 2, 4606088)
        print("# New seqs reduced matrix shape:", [np.shape(mat2[i]) for i in xrange(nchrom)])


# Calculate genetic distance
if args.allcalled:
    df = []
    for i in xrange(nchrom):
        dist_arr = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_all)(mat1[i], mat2[i], i, j) for j in xrange(slens[1]))
        df.append(np.array(dist_arr))
    dist_o_n = np.array(df).sum(0) # dist_o_n[new, old]
   #dist = np.asarray([[ 9, 7], [0,  15], [25,  4]]) # for 3 new and 2 old strain
else:
    dist_c = []
    #len_of_pairs_c = []
    batch_size = NO_JOBS
    batch_no = int(math.ceil(slens[1] / (float(batch_size))))
    for i in xrange(nchrom):
        dist_arr = []
        for k in xrange(batch_no):
            end = (k+1) * batch_size
            if end > slens[1]:
                end = slens[1]
            distance_batch = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_pw)(inputhrseqmat[i], inputnewseqmat[i], non_nuc_masks[0][i], non_nuc_masks[1][i], i, j) for j in xrange(k * batch_size, end))
            dist_arr += distance_batch
        dist_c.append(np.array(dist_arr))
    dist_o_n = np.array(dist_c).sum(0)

    ## v1.04 of the code, this worked
    # for i in xrange(nchrom): # go over chroms
    #     # for each new strain
    #     dist_arr = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_pw)(inputhrseqmat[i], inputnewseqmat[i], non_nuc_masks, i, j) for j in xrange(slens[1]))
    #     dist_c.append(dist_arr)
    #     #len_arr = Parallel(n_jobs=NO_JOBS, max_nbytes=None)(delayed(remain_len)(i, j, ones, non_nuc_masks) for j in xrange(slens[1]))
    #     #len_of_pairs_c.append(len_arr)
    # dist_o_n = np.array(dist_c).sum(0)

    #print(dist_o_n.tolist())
    #len_o_n = np.array(len_of_pairs_c).sum(0)

# print("Old to new lens:")
# print(len_o_n)
# print("Old to new dist:")
# print(dist_o_n)
timing("# Calculated pairwise distance from query to all previous mapped seqs")

# Cluster

# find those that cluster to previous strains
# take the first one, not the closest
# clustering is done on total num of snps
clustered_to_old = {} #new_index: old_index
dist_old = {}
for i, fromnew in enumerate(dist_o_n):
    min_dist = THRESHOLD
    for j, toold in enumerate(fromnew):
        if toold < min_dist:
            clustered_to_old[i] = j
            dist_old[(j, i)] = toold
            break


# do hobohm1 on the remaining ones
# take the first hit
reduced = [] # the new non-homologous seqs
dist_new = {}
#len_new = {}
new_clusters = {} # index : similar to
non_clustered = [x for x in xrange(slens[1]) if x not in clustered_to_old]
if len(non_clustered) > 1:
    reduced = [non_clustered[0]]
    for i in non_clustered[1:]:
        mpdist = THRESHOLD
        close = None
        for j in reduced:
        # calc the distances
            if args.allcalled:
                pdist = np.array([np.not_equal(mat2[k][i], mat2[k][j]).sum(0) for k in range(nchrom)]).sum(0)
            else:
                pdist = np.array([(np.not_equal(inputnewseqmat[k][i], inputnewseqmat[k][j]).sum(0) - np.not_equal(non_nuc_masks[1][k][i], non_nuc_masks[1][k][j]).sum(0)) for k in range(nchrom)]).sum(0)
                # total remaining length for these two
                #len_pw = tot_len - sum([np.logical_or(non_nuc_masks[1][k][i], non_nuc_masks[1][k][j]).sum(0) for k in range(nchrom)])

            dist_new[(i,i)] = 0
            dist_new[(j,j)] = 0
            dist_new[(i,j)] = pdist
            dist_new[(j,i)] = pdist
            # lengths
            # len_new[(i,i)] = tot_len
            # len_new[(j,j)] = tot_len
            # len_new[(i,j)] = len_pw
            # len_new[(j,i)] = len_pw

            if pdist < mpdist:
                new_clusters[i] = j
                break
        if new_clusters.get(i) is None:
            reduced.append(i)

elif len(non_clustered) == 1:
    reduced = non_clustered
    dist_new[(reduced[0],reduced[0])] = 0
else:
    pass

timing("# Clustered the new isolates.")

# if more than one old seq then calc distances
dist_pic_filepath = os.path.join(args.odir, "pairwise.dist.pic")
if slens[0] > 1:
    if args.allcalled:
        # this works as well, but it is slower
        df = []
        for i in xrange(nchrom): #mat1 is the reduced old
            dist_arr = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_all_old)(mat1[i], i, j) for j in xrange(slens[0]))
            df.append(np.array(dist_arr))
        dist_o_o = np.array(df).sum(0)


        # # lets try dividing it into batches
        # df = []
        # batch_size = NO_JOBS
        # batch_no = int(math.ceil(slens[0] / (float(batch_size))))
        # for i in xrange(nchrom):
        #     dist_arr = []
        #     for k in range(batch_no):
        #         distance_batch = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_all_old)(mat1, i, j) for j in xrange(k * batch_size, (k+1) * batch_size))
        #         dist_arr += distance_batch
        #     df.append(np.array(dist_arr))
        # dist_o_o = np.array(df).sum(0)

        ## this works as well, but it is slower
        # df = []
        # for i, chrom in enumerate(mat1): #mat1 is the reduced old
        #     dist_arr = Parallel(n_jobs=28, max_nbytes=None)(delayed(dist_calc_all)(mat1[i], iso) for iso in chrom)
        #     df.append(np.array(dist_arr))
        # dist_o_o = np.array(df).sum(0)

        matrix = dist_o_o.tolist()


    else:
        # get the dist between the old ones from pickle
        if os.path.exists(dist_pic_filepath):
            with open(dist_pic_filepath, "r") as picfile:
                matrix = pickle.load(picfile)

        else:
            # unless there is no pickle, so if old seqs change, delete to recalc the whole matrix
            dist_c_old = []
            #len_of_pairs_c_old = []
            for i in range(nchrom): # go over chroms
                # for each old strain
                dist_arr = Parallel(n_jobs=NO_JOBS)(delayed(dist_calc_pw_old)(inputhrseqmat[i], non_nuc_masks[0][i], i, j) for j in xrange(slens[0]))
                dist_c_old.append(dist_arr)
                #len_arr = Parallel(n_jobs=NO_JOBS, max_nbytes=None)(delayed(remain_len_old)(i, j, ones, non_nuc_masks[0]) for j in xrange(slens[0]))
                #len_of_pairs_c_old.append(len_arr)
            dist_o_o = np.array(dist_c_old).sum(0)
            #len_o_o = np.array(len_of_pairs_c_old).sum(0)
            #print(len_o_o)
            # len normalized distances
            # elmeletben oszthatunk nullaval, de pici a valoszinusege, ha nincs vmi elbaszva
            #matrix = (dist_o_o / len_o_o * 1000000).tolist()
            matrix = dist_o_o.tolist()

else: # only one old seq, dist 0
    matrix = [[0]]

#print(matrix)
timing("# Calculated distances between old sequences.")

# combine all non-redundant into a distance matrix
seqnames = ["template"]
seqnames.extend([x.split(".")[0] for x in oldseqs[1:]])


orig_rows = len(matrix)
newm = []
for i in reduced:
    newm.append([dist_new[(i,j)] for j in reduced])
    seqnames.append(newseqs[i].split(".")[0])
    new_row = []
    for r in xrange(orig_rows):
        matrix[r].append(dist_o_n[i][r])
        new_row.append(dist_o_n[i][r])
    matrix.append(new_row)
full_rows = len(matrix)
for i in xrange(orig_rows, full_rows):
    matrix[i].extend(newm[i - orig_rows])

# Normalised pw distance calculation, to 1Mbps
# orig_rows = len(matrix)
# newm = []
# if args.allcalled:
#     for i in reduced:
#         newm.append([dist_new[(i,j)] for j in reduced])
#         seqnames.append(newseqs[i].split(".")[0])
#         new_row = []
#         for r in xrange(orig_rows):
#             matrix[r].append(dist_o_n[i][r])
#             new_row.append(dist_o_n[i][r])
#         matrix.append(new_row)
# else:
#     for i in reduced:
#         # corrected dist to original length
#         newm.append([(float(dist_new[(i,j)]) / len_new[(i,j)] * 1000000) for j in reduced])
#         seqnames.append(newseqs[i].split(".")[0])
#         new_row = []
#         for r in xrange(orig_rows):
#             matrix[r].append(float(dist_o_n[i][r]) / len_o_n[i][r] * 1000000)
#             new_row.append(float(dist_o_n[i][r]) / len_o_n[i][r] * 1000000)
#         matrix.append(new_row)
# full_rows = len(matrix)
# for i in xrange(orig_rows, full_rows):
#     matrix[i].extend(newm[i - orig_rows])

# dump matrix in pickle for pairwise
# TODO this should be replaced by a dict-dict thing for full flexibility ?
# or just delete if anything is deleted from db
if not args.allcalled:
    # hope that the old files remain available and keep the same order in the input non-redundant list
    # if not, then the pickle file should be deleted in order to the old-old dists to be re-calculated
    dist_pic_filepath += suffix
    with open(dist_pic_filepath, "w") as picklingfile:
        pickle.dump(matrix, picklingfile)



# print dist matrix in phylip
seqid2name = {}
seq_id = ""
with open(outputmat, "w") as matfile:
    #print("\t","\t".join(seqnames), file=matfile)
    print("  {0}".format(len(matrix)), file=matfile)
    for r, row in enumerate(matrix):
        seq_id = "I{:06}".format(r+1)
        seqid2name[seq_id] = seqnames[r]
        print("{:<10}".format(seq_id), end = "", file=matfile)
        for e in row[:-1]:
            print('{0:.0f}'.format(e), end = "\t", file=matfile)
        print('{0:.0f}'.format(row[-1]), file=matfile)

with open(os.path.join(args.odir, "seqid2name.pic"), "w") as pf:
    pickle.dump(seqid2name, pf)

timing("# Constructed distance matrix.")

# save hr reduced file paths to file
hrseqs = oldseqs[1:]
for i in reduced:
    hrseqs.append(newseqs[i])
if args.allcalled:
    hrfile = open(os.path.join(args.odir,"non-redundant.all.lst{0}".format(suffix)), "w")
else:
    hrfile = open(os.path.join(args.odir,"non-redundant.pw.lst{0}".format(suffix)), "w")
#print("{0}/pt_{1}.fa".format(args.odir, args.template), file=hrfile)
print(oldseqs[0], file=hrfile) # template is the 1st
print("\n".join(hrseqs), file=hrfile)
hrfile.close()
timing("# Saved new isolates to file.")

# add to sequences db
cur.executemany('''INSERT OR IGNORE INTO sequences (sra_id) VALUES (?)''', zip([x.split('.')[0] for x in hrseqs],))

if args.allcalled:
    clfilename = os.path.join(args.odir, "clusters.all.pic")
else:
    clfilename = os.path.join(args.odir, "clusters.pw.pic")

clusters = {}
cluster_insert = []
if os.path.exists(clfilename):
    with open(clfilename, 'rb') as pf:
        clusters = pickle.load(pf)

# make vcf files from the redundant ones and add clusters
if clustered_to_old:
    homologous = seq_to_homol(clustered_to_old)
    for key in homologous.keys():
        ref = (0, key)
        add_clusters(ref, homologous.get(key))
        if args.vcf:
            create_vcf(ref, homologous.get(key))

if new_clusters:
    homologous = seq_to_homol(new_clusters)
    for key in homologous.keys():
        ref = (1, key)
        add_clusters(ref, homologous.get(key))
        if args.vcf:
            create_vcf(ref, homologous.get(key))
            timing("# Made vcf files.")

clfilename += suffix
if clusters:
    with open(clfilename, 'wb') as pf:
        pickle.dump(clusters, pf)
    cur.executemany('''INSERT OR IGNORE INTO sequences VALUES (?,?,?)''', cluster_insert)
    conn.commit()
    conn.close()

    timing("# Pickled redundant sequences.")

else:
    if os.path.exists(clfilename) and args.debug:
        #delete the suffixed pickle if it's there
        try:
            os.unlink(clfilename)
        except IOError as e:
            print("# Warning: {0} couldn't be deleted, {1}".format(clfilename, e), file=sys.stderr)

# Finish up
timing("# Distance calculation is finished.")
print("DONE", file=sys.stderr)
sys.exit(0)
