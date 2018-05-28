#!/usr/bin/env python2.7
from __future__ import print_function
import sys, time
import os, gc
import math
import numpy as np
import array
import argparse
from operator import itemgetter
import gzip
import tempfile
from joblib import Parallel, delayed, dump, load
import sqlite3
try:
    import cPickle as pickle
except ImportError:
    import pickle
import config

# quick hack
base_path = os.path.dirname(sys.argv[0]).rsplit("/",1)[0]
MAIN_SQL_DB = os.path.join(base_path, "results_db/evergreen.db")
no_jobs = 20
MEM_AVAIL = 80 # Gb
THRESHOLD = 10

"""
Change:
    removed N content test, KMA id% -> not more than 0.15
    -n - and -hr -: dont uses files, relies on the db
"""

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

def parallel_opt(no_jobs, tot_len):
    if slens[1] > 200:
        no_jobs = int( min(no_jobs, math.floor(MEM_AVAIL * 0.8 / (sum(slens[1]) * tot_len / 10.0**9)) - 1))
        no_jobs	= max(1, no_jobs)
    return no_jobs

def read_encode_univ(filename, tot_len):
    fp = os.path.join(args.odir, filename)
    if os.path.exists(fp):
        entries = zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)])
        strain = "".join(entries[0])
        if tot_len is None:
            tot_len = len(strain)

        encodedinput = np.zeros((tot_len), dtype=np.int8)
        for i in xrange(tot_len):
            try:
                encodedinput[i] = nuc2num[strain[i]]
            except KeyError:
                pass
        return encodedinput
    else:
        return None

def dist_calc_pw(s1, s2, i, S, no_compared_strains):
    if slens[0] < i + S:
        S = slens[0] - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.int32)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            #print(k,l)
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0) - np.not_equal(s1[l,]!= 0, s2[k,]!= 0).sum(0)
    return dist_t

def dist_calc_all(s1, s2, i, S, no_compared_strains):
    if slens[0] < i + S:
        S = slens[0] - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.int32)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            #print(k,l)
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0)
    return dist_t

def seq_to_homol(cluster):
    """Convert cluster dict to dict of lists."""
    redict = {} # homol : [seqs]
    for value in set(cluster.values()):
        redict[value] = [key for key in cluster if cluster.get(key) == value]
    return redict

def add_clusters(ref, isolates):
    # clusters = { hr-name0 : [ (isolate0, distance), (isolate1, distance), ... ], ... }
    homolname = None
    if ref[0] == 0:
        # clustered to an old one
        if ref[1] == 0: # template
            homolname = "template"
        else:
            homolname = oldseqs[ref[1]].split(".")[0]

        iso_cur.execute('''SELECT count(*) from sequences where repr_id=?''', (homolname,))
        count = iso_cur.fetchone()
        if count is None:
            cluster_increase.append((templ, homolname, len(isolates) + 1))
        else:
            cluster_increase.append((templ, homolname, len(isolates)))
    else:
        # clustered to new
        homolname = newseqs[ref[1]].split(".")[0]
        cluster_increase.append((templ, homolname, len(isolates) + 1))

    for ind in isolates:
        isoname = newseqs[ind].split(".")[0]
        if ref[0] == 0:
            dist = dist_old[(ref[1], ind)]
        else:
            dist = dist_new[(ref[1], ind)]
        # prepare db insert
        cluster_insert.append((homolname, str(dist), isoname))
    return

# Start time to keep track of progress
t0 = time.time()

# Parse command line options
parser = argparse.ArgumentParser(
    description='Calculates genetic distance between two set of input sequences, reduces homology')
parser.add_argument("-hr", dest="old", help="read non-homologous sequences from", metavar="INFILE_HR")
parser.add_argument("-n", dest="new", help="read new sequences from", metavar="INFILE_NEW")
parser.add_argument("-o", "--outputdir", dest="odir", help="write to DIR", metavar="DIR")
parser.add_argument("-m", "--outputmat", dest="outputfilenamemat", help="write to OUTMAT", metavar="OUTMAT")
parser.add_argument("-a", "--allcalled", dest="allcalled", action="store_true", help="Only use positions called in all strains")
parser.add_argument("-d", "--debug", dest="debug", action="store_true", help="Debug: use .t suffix")
parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Quiet")
args = parser.parse_args()

suffix = ""
if args.debug:
    suffix = ".t"

# File with reads
elif args.old is None or args.new is None:
    exiting('No input filelist was provided')

# path to results dir ../results_db/template/
if args.odir is None:
    exiting('Output directory is needed')

# File for distance matrix output
templ = os.path.basename(args.odir)
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

# open database
# MAIN
conn = sqlite3.connect(MAIN_SQL_DB)
conn.execute("PRAGMA foreign_keys = 1")
conn.commit()
cur = conn.cursor()

# for template
mode = 'pw'
if args.allcalled:
    mode = 'all'

db_path = os.path.join(args.odir, "isolates.{0}.db{1}".format(mode, suffix))
non_red_pic = os.path.join(args.odir, "non-red.{0}.pic{1}".format(mode, suffix))
hr_matrix_npy = os.path.join(args.odir, "hr-matrix.{0}{1}.npy".format(mode, suffix))


iso_conn = sqlite3.connect(db_path)
iso_conn.execute("PRAGMA foreign_keys = 1")
iso_cur = iso_conn.cursor()


# New strains
newseqs = []

if args.new == "-":
    # TODO DONE isolates list from db
    iso_cur.execute('''SELECT sra_id from sequences where repr_id='N';''')
    rows = iso_cur.fetchall()
    if rows is not None:
        for r in rows:
            fp = os.path.join(args.odir, r[0])
            newseqs.append("{}.fa".format(r[0]))
else:
    try:
        f = open(args.new, "r")
    except IOError:
        exiting("List of new isolates not found.")
    for l in f:
        l = l.strip()
        if l == "":
            next #avoid empty line
        fp = os.path.join(args.odir, l)
        if os.path.exists(fp) and os.path.getsize(fp) > 0:
            newseqs.append(l)
    f.close()

# Exit if no new
if not newseqs:
    exiting("No new sequences.")

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


# homology reduced old isolates
oldseqs = [] # hr reduced seqs
# TODO DONE get these from db, where repr_id is NULL
non_red_seqs = []
re_calc = False
# pickled list
if os.path.exists(non_red_pic):
    with open(non_red_pic, "rb") as fp:
        non_red_seqs = pickle.load(fp)

# start with the template
db_seqs = ["pt_{}.fa".format(templ)]
iso_cur.execute('''SELECT sra_id from sequences where repr_id is NULL and sra_id!='template';''')
rows = iso_cur.fetchall()
for r in rows:
    db_seqs.append("{}.fa".format(r[0]))

# if input file, compare to db and seqs pic
# order of importance: inputfile > db > pic
if args.old != '-':
    try:
        f = open(args.old, "r")
    except IOError as e:
        exiting("Redundant file list not found.")
    for l in f:
        l = l.strip()
        oldseqs.append(l)
    f.close()

    # not the same as the saved matrix rows
    if non_red_seqs != oldseqs:
        re_calc = True
    else:
        # an isolate was deleted
        if len(non_red_seqs) != len(db_seqs):
            re_calc = True
else:
    # no input file, but there are the pickled ones
    if not non_red_seqs:
        # first time template
        oldseqs = db_seqs
        re_calc = True
    else:
        oldseqs = non_red_seqs

# but if it doesnt match the db
if not re_calc:
    for iso in oldseqs:
        if iso not in db_seqs:
            re_calc = True
            oldseqs = db_seqs
            break

timing("# Read inputfiles from the file lists.")

#TODO DONE
slens = [len(oldseqs), len(newseqs)] # Number of strains in each list
inputhrseqmat = None
inputnewseqmat = None
tot_len = None

# load and encode new isolates
# might not need to be parallel (overhead)
if slens[1] > 30:
    arrays = Parallel(n_jobs=no_jobs)(delayed(read_encode_univ)(isolatefile, None) for isolatefile in newseqs)
    # dump as a memmap and concat
    tot_len = np.shape(arrays[0])[0]
    for i, arr in enumerate(arrays):
        if i == 0:
            inputnewseqmat = np.zeros((slens[1], tot_len), dtype=np.int8)
            inputnewseqmat[0,:] = arr[:]
        else:
            inputnewseqmat[i,:] = arr[:]
    del arrays[:]
else:
    for i, isolatefile in enumerate(newseqs):
        if i == 0:
            tmp_np = read_encode_univ(isolatefile, tot_len)
            tot_len = np.shape(tmp_np)[0]
            inputnewseqmat = np.zeros((slens[1], tot_len), dtype=np.int8)
            inputnewseqmat[0,:] = tmp_np[:]
        else:
            inputnewseqmat[i,:] = read_encode_univ(isolatefile, tot_len)[:]

# save as npy
np.save(os.path.join(args.odir, "new-matrix.npy"), inputnewseqmat, allow_pickle=True, fix_imports=True)

# load / encode old isolates
if not re_calc and os.path.exists(hr_matrix_npy):
    # load into memory
    inputhrseqmat = np.load(hr_matrix_npy, mmap_mode=None, allow_pickle=True, fix_imports=True)
else:
    # read in and encode the seqs in parallel
    arrays = Parallel(n_jobs=no_jobs)(delayed(read_encode_univ)(isolatefile, tot_len) for isolatefile in oldseqs)
    # construct in memory
    iso_not_found = []
    for i, arr in enumerate(arrays):
        if i == 0:
            # do it in the memory
            inputhrseqmat = np.zeros(shape=(slens[0],tot_len), dtype=np.int8)
        if arr is not None:
            # if the file existed
            inputhrseqmat[i,:] = arr[:]
        else:
            iso_not_found.append(i)

    # house-keeping, removing the sequence where the file was not found
    if iso_not_found:
        np.delete(inputhrseqmat, iso_not_found, axis=0)
        for i in iso_not_found[::-1]:
            del oldseqs[i]
    # save as npy
    np.save(hr_matrix_npy, inputhrseqmat, allow_pickle=True, fix_imports=True)

    del arrays[:]

# sanity check
if len(oldseqs) != np.shape(inputhrseqmat)[0]:
    try:
        os.unlink(hr_matrix_npy)
    except OSError:
        pass
    exiting("Number of hr seqs doesnt agree with matrix shape.")

timing("# Loaded sequence matrices into memory.")

# trunc the allcalled matrices
m = None
if args.allcalled:
    # create the 'all' mask from the arrays
    m = np.logical_and((inputhrseqmat != 0).all(axis=0), (inputnewseqmat != 0).all(axis=0))
    inputhrseqmat = inputhrseqmat.T[m].T
    inputnewseqmat = inputnewseqmat.T[m].T

# conserved positions for both type of distance calculation
if slens[0] > 1 and slens[1] > 1:
    m = np.logical_not(np.logical_and(np.all(inputhrseqmat == inputhrseqmat[0,:], axis = 0),np.all(inputnewseqmat == inputnewseqmat[0,:], axis = 0)))
    inputhrseqmat = inputhrseqmat.T[m].T
    inputnewseqmat = inputnewseqmat.T[m].T

del m

# update lengths and isolate counts
slens = [len(oldseqs), len(newseqs)]
tot_len = np.shape(inputnewseqmat)[1]


timing("# Removed non-informative positions from matrices.")
# print(oldseqs)
# print(newseqs)

if not args.quiet:
    print("# Total length: %s" %(tot_len), file=sys.stdout)
    print("# Number of strains: %s" % (slens), file=sys.stdout)

# print(np.info(inputhrseqmat))
# print(np.info(inputnewseqmat))

# TODO DONE
# dump hr array to /dev/shm
temp_folder = tempfile.mkdtemp(prefix='ever_joblib_', dir=config.NUMPY_MEMMAP_DIR)
hr_memmap_fn = os.path.join(temp_folder, 'hr_matrix.npy')
dump(inputhrseqmat, hr_memmap_fn)

# delete from memory, do garbage collection
del inputhrseqmat
gc.collect()

# load as a memmap
hrseqmat_memmap = load(hr_memmap_fn, mmap_mode='r')
# print(np.info(hrseqmat_memmap))
timing("# Dumped hr matrix to disk")

# calculate genetic distance between old and new isolates
no_jobs = min(slens[0], no_jobs)
no_jobs = parallel_opt(no_jobs, tot_len)
batch_size = int(math.ceil(slens[0]/float(no_jobs)))
dist_calc_func = dist_calc_pw
if args.allcalled:
    dist_calc_func = dist_calc_all

dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_func)(hrseqmat_memmap, inputnewseqmat, i, batch_size, slens[1]) for i in xrange(0,slens[0],batch_size))

# put it together, do a np.where on it
dist_nw_o = dist_arr[0]
for np_arr in dist_arr[1:]:
    dist_nw_o = np.concatenate((dist_nw_o,np_arr), axis=1)
    #dist = np.asarray([[ 9, 7], [0,  15], [25,  4]]) # for 3 new and 2 old strain

del dist_arr[:]

timing("# Calculated pairwise distance from query to all previous mapped seqs")

# print(dist_nw_o.tolist())
# collect the old - new clusters
clustered_to_old = {}
dist_old = {}
min_dist = THRESHOLD
nw, hr = np.where(dist_nw_o < min_dist)
pre = -1
for i, hr_index in enumerate(hr):
    # just take the first hit from the hr-s
    if nw[i] != pre:
        clustered_to_old[nw[i]] = hr_index
        dist_old[(hr_index, nw[i])] = dist_nw_o[nw[i]][hr_index]
        pre = nw[i]

# print(clustered_to_old)
# print(dist_old)

# do hobohm1 on the remaining ones
# take the first hit
reduced = [] # the new non-homologous seqs
dist_new = {}
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
                pdist = np.not_equal(inputnewseqmat[i,], inputnewseqmat[j,]).sum(0)
            else:
                pdist = np.not_equal(inputnewseqmat[i,], inputnewseqmat[j,]).sum(0) - np.not_equal(inputnewseqmat[i,]!= 0, inputnewseqmat[j,]!= 0).sum(0)

            dist_new[(i,i)] = 0
            dist_new[(j,j)] = 0
            dist_new[(i,j)] = pdist
            dist_new[(j,i)] = pdist

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
    # get the dist between the old ones from pickle
    if not args.allcalled and not re_calc and os.path.exists(dist_pic_filepath):
        with open(dist_pic_filepath, "rb") as picfile:
            matrix = pickle.load(picfile)

    else:
        dist_arr = []
        dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_func)(hrseqmat_memmap, hrseqmat_memmap, i, batch_size, slens[0]) for i in xrange(0,slens[0],batch_size))

        # put it together, do a np.where on it
        dist_o_o = dist_arr[0]
        for np_arr in dist_arr[1:]:
            dist_o_o = np.concatenate((dist_o_o,np_arr), axis=1)

        matrix = dist_o_o.tolist()

else: # only one old seq, dist 0
    matrix = [[0]]

# print(matrix)
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
        matrix[r].append(dist_nw_o[i][r])
        new_row.append(dist_nw_o[i][r])
    matrix.append(new_row)
full_rows = len(matrix)
for i in xrange(orig_rows, full_rows):
    matrix[i].extend(newm[i - orig_rows])

# dump matrix in pickle for pairwise
if not args.allcalled:
    # if order of files changes, then the old-old dists are re-calculated
    dist_pic_filepath += suffix
    with open(dist_pic_filepath, "wb") as picklingfile:
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

with open(os.path.join(args.odir, "seqid2name.{}.pic".format(mode)), "w") as pf:
    pickle.dump(seqid2name, pf)

timing("# Constructed distance matrix.")

# TODO try add all new sequences to db
if args.new != "-":
    iso_cur.executemany('''INSERT OR IGNORE INTO sequences (sra_id, repr_id) VALUES (?, 'N')''', zip([x.split('.')[0] for x in newseqs],))
    iso_conn.commit()

# save hr reduced file paths to pickle
hrseqs_update = []

for i in reduced:
    hrseqs_update.append(newseqs[i].split('.')[0])
    oldseqs.append(newseqs[i])
with open(non_red_pic, "wb") as fp:
    pickle.dump(oldseqs, fp)

if hrseqs_update:
    iso_cur.executemany('''UPDATE sequences set repr_id=NULL where sra_id=? and repr_id='N';''', zip(hrseqs_update,))
    iso_conn.commit()

timing("# Saved new non-redundant isolates to file.")

# re-load full numpy matrices
inputnewseqmat = np.load(os.path.join(args.odir, "new-matrix.npy"), mmap_mode=None, allow_pickle=True, fix_imports=True)

cluster_insert = []
cluster_increase = []
del_iso = []
# make clusters
if clustered_to_old:
    homologous = seq_to_homol(clustered_to_old)
    for key in homologous.keys():
        ref = (0, key)
        add_clusters(ref, homologous.get(key))

    del_iso = clustered_to_old.keys()

if new_clusters:
    homologous = seq_to_homol(new_clusters)
    for key in homologous.keys():
        ref = (1, key)
        add_clusters(ref, homologous.get(key))

    del_iso += new_clusters.keys()

if cluster_insert:

    inputnewseqmat = np.delete(inputnewseqmat, del_iso, axis=0)

    iso_cur.executemany('''UPDATE sequences set repr_id=?, distance=? where sra_id=? and repr_id='N';''', cluster_insert)
    iso_conn.commit()
    iso_conn.close()

    timing("# Saved clusters to the db.")

    if not args.allcalled and not args.debug:
        # add cluster size changes to db
        # make table for the cluster_sizes
        cur.execute('''CREATE TABLE IF NOT EXISTS clusters
            (db_id INTEGER PRIMARY KEY,
            template TEXT,
            repr_id TEXT,
            change INTEGER,
            ctime DATETIME DEFAULT CURRENT_DATE,
            UNIQUE(template, repr_id, ctime)
        )''')
        conn.commit()

        try:
            cur.executemany('''INSERT OR IGNORE INTO clusters (template, repr_id, change) VALUES (?,?,?)''', cluster_increase)
            conn.commit()
        except sqlite3.Error as e:
            print("Warning: cluster size SQL update failed.", str(e), file=sys.stderr)

# Finish up
conn.close()

if inputnewseqmat.shape[0]:
    inputhrseqmat = np.load(hr_matrix_npy, mmap_mode=None, allow_pickle=True, fix_imports=True)
    inputhrseqmat = np.concatenate((inputhrseqmat, inputnewseqmat), axis = 0)
    np.save(hr_matrix_npy, inputhrseqmat, allow_pickle=True, fix_imports=True)

    del inputhrseqmat

    timing("# Updated the redundant sequences npy.")

# cleanup
del inputnewseqmat
gc.collect()

try:
    os.unlink(os.path.join(args.odir, "new-matrix.npy"))
    os.unlink(hr_memmap_fn)
    os.rmdir(temp_folder)
except OSError:
    import shutil
    shutil.rmtree(temp_folder)


timing("# Distance calculation is finished.")
print("DONE", file=sys.stderr)
sys.exit(0)
