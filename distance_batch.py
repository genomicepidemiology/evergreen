#!/usr/bin/env python2.7
from __future__ import print_function, division
import sys, time
import os, gc
import math
import numpy as np
import array
import argparse
import glob
from operator import itemgetter
import gzip
import tempfile
from multiprocessing import cpu_count
from joblib import Parallel, delayed, load, dump
import sqlite3
try:
    import cPickle as pickle
except ImportError:
    import pickle
import config

cpu_count = cpu_count() // 4

"""
Change: from v2.11
    - isolate cllustered to closest or first closest match (if ties)
    - if too many columns would be removed for pristine, calc known_frac on hr+nw
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

def read_encode_univ(filename, tot_len, odir):
    fp = os.path.join(odir, filename)
    if os.path.exists(fp):
        entries = zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)])
        strain = "".join(entries[0])
        if tot_len is None:
            tot_len = len(strain)

        encodedinput = np.zeros((tot_len), dtype=np.int8)

        n_count = strain.count("N")
        if ((n_count / tot_len) > 0.100):
            print("# {0} entry has too many N calls".format(entries[1][0]), file=sys.stderr)
            # mark it empty, max num (signed int)
            encodedinput[0] = 127
        else:
            for i in xrange(tot_len):
                try:
                    encodedinput[i] = nuc2num[strain[i]]
                except KeyError:
                    pass
        return encodedinput
    else:
        return [127]

def old_arrays_to_volumes(arraymatrix):
    '''Organize old arrays into volumes on disk'''
    no_vol = int(math.ceil(np.shape(arraymatrix)[0]/config.VOL_SIZE))
    for i in range(no_vol):
        matrix_npy_vol = hr_matrix_npy.format(i+1)
        np.save(matrix_npy_vol, arraymatrix[i * config.VOL_SIZE:(i+1)*config.VOL_SIZE], allow_pickle=True, fix_imports=True)

def dist_calc_pw(s1, s2, i, S, no_compared_strains):
    if vol_len < i + S:
        S = vol_len - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.uint16)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            #print(k,l)
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0) - np.not_equal(s1[l,]!= 0, s2[k,]!= 0).sum(0)
    return dist_t

def dist_calc_all(s1, s2, i, S, no_compared_strains):
    if vol_len < i + S:
        S = vol_len - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.uint16)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0)
    return dist_t

def dist_calc_square(s1, i, S, no_compared_strains):
    if no_compared_strains < i + S:
        S = no_compared_strains - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.uint16)
    for j in xrange(S):
        l = i+j
        # want to calculate only below the diagonal (lower triangle of distance matrix)
        for k in xrange(l+1,no_compared_strains):
            if not pristine:
                dist_t[k,j] = np.not_equal(s1[l,], s1[k,]).sum(0) - np.not_equal(s1[l,]!= 0, s1[k,]!= 0).sum(0)
            else:
                dist_t[k,j] = np.not_equal(s1[l,], s1[k,]).sum(0)
    return dist_t

def dist_calc_lower(s1, s2, i, S, no_compared_strains, total_batched):
    if total_batched < i + S:
        S = total_batched - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.uint16)
    for j in xrange(S):
        l = i+j
        for k in xrange(no_compared_strains):
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0) - np.not_equal(s1[l,]!= 0, s2[k,]!= 0).sum(0)
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
        count = int(iso_cur.fetchone()[0])
        if count == 0:
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
parser.add_argument("-a", "--allcalled", dest="allcalled", action="store_true", help="Only use positions called in all strains < 100, then probabilistic ")
parser.add_argument("-d", "--debug", dest="debug", action="store_true", help="Debug: use .t suffix")
parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Quiet")
args = parser.parse_args()

suffix = ""
if args.debug:
    suffix = ".t"

if args.allcalled:
    method = "all"
else:
    method = "pw"

db_path = os.path.join(args.odir, "isolates.{}.db{}".format(method, suffix))
non_red_pic = os.path.join(args.odir, "non-red.{}.pic{}".format(method, suffix))
hr_matrix_npy = os.path.join(args.odir, "hr-matrix.{}{}.{}.npy".format(method, suffix, "{}"))
dist_pic_filepath = os.path.join(args.odir, "pairwise.dist.{}.pic".format(method))


base_path = os.path.dirname(os.path.realpath(args.odir))
main_sql_db = os.path.join(base_path, "evergreen.db")

# adjust number of jobs
if os.environ.get('PBS_NP') is not None:
    # we are on a moab HPC cluster
    cpu_count = int(os.environ.get('PBS_NP'))//4
no_jobs = cpu_count

# File with reads
if args.old is None or args.new is None:
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
conn = sqlite3.connect(main_sql_db)
conn.execute("PRAGMA foreign_keys = 1")
conn.commit()
cur = conn.cursor()

# TEMPLATE
iso_conn = sqlite3.connect(db_path)
iso_conn.execute("PRAGMA foreign_keys = 1")
iso_cur = iso_conn.cursor()


# New strains
newseqs = []
templ_update = []
runs_update = []
seqs_update = []

if args.new == "-":
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
    else:
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

# Homology reduced old isolates
oldseqs = [] # hr reduced seqs
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

## Encode new and load old isolates
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

inputhrseqmat = None
inputnewseqmat = None
tot_len = None
slens = [len(oldseqs), len(newseqs)] # Number of strains in each list

pristine = False
if args.allcalled and slens[0] <= config.PHOLD:
    pristine = True
# load and encode new isolates
# might not need to be parallel (overhead)
if slens[1] > 30:
    arrays = Parallel(n_jobs=no_jobs)(delayed(read_encode_univ)(isolatefile, None, args.odir) for isolatefile in newseqs)
    # dump as a memmap and concat
    tot_len = np.shape(arrays[0])[0]
    for i, arr in enumerate(arrays):
        if i == 0:
            inputnewseqmat = np.zeros((slens[1], tot_len), dtype=np.int8)
            inputnewseqmat[0,:] = arr[:]
        else:
            inputnewseqmat[i,:] = arr[:]
        if arr[0] == 127:
            sra_id = newseqs[i].split(".")[0]
            templ_update.append((0, sra_id, templ))
            runs_update.append((2, sra_id))
            seqs_update.append((sra_id,))
    del arrays[:]
else:
    for i, isolatefile in enumerate(newseqs):
        if i == 0:
            tmp_np = read_encode_univ(isolatefile, tot_len, args.odir)
            tot_len = np.shape(tmp_np)[0]
            inputnewseqmat = np.zeros((slens[1], tot_len), dtype=np.int8)
            inputnewseqmat[0,:] = tmp_np[:]
        else:
            inputnewseqmat[i,:] = read_encode_univ(isolatefile, tot_len, args.odir)[:]
        if inputnewseqmat[i,0] == 127:
            sra_id = newseqs[i].split(".")[0]
            templ_update.append((0, sra_id, templ))
            runs_update.append((2, sra_id))
            seqs_update.append((sra_id,))

# TODO DONE remove > 10% N content
if templ_update:
    # reach back and update runs table in the main DB to 2 as non-included
    try:
        cur.executemany('''UPDATE templates SET qc_pass=? WHERE sra_id=? and template=?''', templ_update)
        cur.executemany('''UPDATE runs SET included=? WHERE sra_id=?''', runs_update)
        conn.commit()
        if args.new == "-":
            iso_cur.executemany('''DELETE FROM sequences WHERE sra_id=? and repr_id='N';''', seqs_update)
            iso_conn.commit()
    except sqlite3.Error:
        print("Warning: SQL update failed.", file=sys.stderr)

    bad_iso = []
    for rec in seqs_update:
        ind = newseqs.index("{}.fa".format(rec[0]))
        if ind > -1:
            bad_iso.append(ind)
    # delete from matrix
    inputnewseqmat = np.delete(inputnewseqmat, bad_iso, axis=0)
    for i in sorted(bad_iso, reverse=True):
        del newseqs[i]
# save as npy
np.save(os.path.join(args.odir, "new-matrix.npy"), inputnewseqmat, allow_pickle=True, fix_imports=True)

# load / encode old isolates
hr_matrix_npy_vol = hr_matrix_npy.format("1")

# if volumes don't exists or we want to remake them
if re_calc or not os.path.exists(hr_matrix_npy_vol):
    no_jobs = min(cpu_count, slens[0])
    arrays = Parallel(n_jobs=no_jobs)(delayed(read_encode_univ)(isolatefile, tot_len, args.odir) for isolatefile in oldseqs)
    # construct in memory
    iso_not_found = []
    for i, arr in enumerate(arrays):
        if i == 0:
            # do it in the memory
            inputhrseqmat = np.zeros(shape=(slens[0],tot_len), dtype=np.int8)
        if arr[0] != 127:
            # if the file existed
            inputhrseqmat[i,:] = arr[:]
        else:
            iso_not_found.append(i)

    # house-keeping, removing the sequence where the file was not found
    if iso_not_found:
        inputhrseqmat = np.delete(inputhrseqmat, iso_not_found, axis=0)
        for i in iso_not_found[::-1]:
            del oldseqs[i]

    old_arrays_to_volumes(inputhrseqmat)

    # clean up
    del arrays[:]
    del inputhrseqmat
    gc.collect()

    timing("# Volumized arrays.")

# update set sizes, volume number
slens = [len(oldseqs), len(newseqs)] # Number of strains in each list
no_vol = int(math.ceil(slens[0]/config.VOL_SIZE))
if no_vol != len(glob.glob(hr_matrix_npy.format("*"))):
    exiting("Number of volumes dont match the expected {}".format(no_vol))

# load 1st volume
inputhrseqmat = np.load(hr_matrix_npy_vol, mmap_mode=None, allow_pickle=True, fix_imports=True)
vol_len = np.shape(inputhrseqmat)[0]

## Remove uncertain regions
# determine masks, probabilistic (< 100 is all (also includes the new isolates)), remove: False
std_dev = None
mean_a = None
m = None
if args.allcalled:
    if pristine:
        m = np.logical_and((inputhrseqmat != 0).all(axis=0), (inputnewseqmat != 0).all(axis=0))
        if (( 1 - (m.sum() / tot_len)) > config.DYHOLD):
            pristine = False
            known_frac = np.round_((((inputhrseqmat != 0).sum(0)) + ((inputnewseqmat != 0).sum(0))) / float(slens[0]+slens[1]), 5)
            m = (known_frac >= config.CLEAN)
    else:
        known_frac = np.round_((inputhrseqmat != 0).sum(0) / vol_len, 5)
        std_dev = np.std(known_frac, dtype=np.float64)
        mean_a = np.mean(known_frac, dtype=np.float64)
        #threshold = mean_a - std_dev
        threshold = config.CLEAN
        m = (known_frac >= threshold)
    inputhrseqmat = inputhrseqmat.T[m].T
    inputnewseqmat = inputnewseqmat.T[m].T

# conserved positions for both type of distance calculation, if there is only one volume (ie limited set size)
cons_m = None
if no_vol == 1 and slens[0] > 1 and slens[1] > 1:
    cons_m = np.logical_not(np.logical_and(np.all(inputhrseqmat == inputhrseqmat[0,:], axis = 0),np.all(inputnewseqmat == inputnewseqmat[0,:], axis = 0)))
    inputhrseqmat = inputhrseqmat.T[cons_m].T
    inputnewseqmat = inputnewseqmat.T[cons_m].T
del cons_m

# update sequence length to reflect masking
tot_len = np.shape(inputnewseqmat)[1]

timing("# Removed non-informative positions from matrices.")

if not args.quiet:
    print("# Known fraction mean, SD: {} {}".format(mean_a, std_dev), file=sys.stdout)
    print("# Minimum known bases per position: {}".format(config.CLEAN), file=sys.stdout)
    print("# Total length: {}".format(tot_len), file=sys.stdout)
    print("# Number of strains: {}".format(slens), file=sys.stdout)

# print(np.info(inputhrseqmat))
# print(np.info(inputnewseqmat))

# which function to use for distance calculation
dist_calc_func = dist_calc_pw
if pristine:
    dist_calc_func = dist_calc_all
# dump hr array to /dev/shm
temp_folder = tempfile.mkdtemp(prefix='ever_joblib_', dir=config.NUMPY_MEMMAP_DIR)
hr_memmap_fn = os.path.join(temp_folder, 'hr_matrix.{}.mmap')

dist_nw_o = np.zeros(shape=(slens[1],slens[0]), dtype=np.uint16)
#offset for complete distance matrix
o_x = 0
for vol in range(no_vol):
    if vol != 0:
        # load the volume and mask it
        hr_matrix_npy_vol = hr_matrix_npy.format(vol+1)
        inputhrseqmat = np.load(hr_matrix_npy_vol, mmap_mode=None, allow_pickle=True, fix_imports=True)
        vol_len = np.shape(inputhrseqmat)[0]
        if m is not None:
            inputhrseqmat = inputhrseqmat.T[m].T

    # dump hr array to /dev/shm
    dump(inputhrseqmat, hr_memmap_fn.format(vol+1))

    # delete from memory, do garbage collection
    del inputhrseqmat
    gc.collect()

    hrseqmat_memmap = load(hr_memmap_fn.format(vol+1), mmap_mode='r')

    # print(np.info(hrseqmat_memmap))
    timing("# Dumped {} volume to disk".format(vol+1))

    # calculate genetic distance between old and new isolates
    no_jobs = min(vol_len, cpu_count)
    batch_size = int(math.ceil(vol_len/float(no_jobs)))

    dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_func)(hrseqmat_memmap, inputnewseqmat, i, batch_size, slens[1]) for i in xrange(0,vol_len,batch_size))

    # put it together, slower but less memory usage
    # dist_arr = [[np.array(shape=(slens[1], batch_size), dtype=np.uint16)],[]]
    #dist = np.asarray([[ 9, 7], [0,  15], [25,  4]]) # for 3 new and 2 old strain
    for np_arr in dist_arr:
        # no of samples in actual batch
        arr_y = np.shape(np_arr)[1]
        arr_x = np.shape(np_arr)[0] # == slens[1]
        for k in xrange(arr_x):
            for l in xrange(arr_y):
                dist_nw_o[k,l+o_x] = np_arr[k,l]
        o_x += arr_y

    # clean up memory and disk
    os.unlink(hr_memmap_fn.format(vol+1))
    del hrseqmat_memmap

del dist_arr[:]

timing("# Calculated pairwise distance from query to all previous mapped seqs")


## Collect the old - new clusters
clustered_to_old = {}
dist_old = {}
nw, hr = np.where(dist_nw_o < config.THRESHOLD)
pre = -1
for i, hr_index in enumerate(hr):
    # take lowest dist or first from equally best
    if nw[i] != pre or dist_nw_o[nw[i]][hr_index] < dist_old[(clustered_to_old[nw[i]], nw[i])]:
        clustered_to_old[nw[i]] = hr_index
        dist_old[(hr_index, nw[i])] = dist_nw_o[nw[i]][hr_index]
        pre = nw[i]

# print(clustered_to_old)
# print(dist_old)

# Insert the found pairs, remove the new clustered ones
cluster_insert = []
cluster_increase = []
del_iso = []
if clustered_to_old:
    homologous = seq_to_homol(clustered_to_old)
    for key in homologous.keys():
        ref = (0, key)
        add_clusters(ref, homologous.get(key))

    del_iso = clustered_to_old.keys()

    # reduce the size of the new sequence matrix and the distance matrix nw-o
    inputnewseqmat = np.delete(inputnewseqmat, del_iso, axis=0)
    dist_nw_o = np.delete(dist_nw_o, del_iso, axis=0)

# do hobohm1 on the remaining ones
## Calculate new-new distances and hobohm1 on the remaining ones
reduced = [] # the new non-homologous seqs
dist_new = {} # distance of the clustered sequences
new_clusters = {} # index : similar to
non_clustered = [x for x in xrange(slens[1]) if x not in del_iso]
upd_new_len = len(non_clustered)
if upd_new_len != np.shape(inputnewseqmat)[0]:
    exiting("Oops.")

dist_nw_nw = None
# only one remaining
if upd_new_len == 1:
    reduced = [non_clustered[0]]
    dist_nw_nw = np.array([[0]], dtype=np.int8)
elif upd_new_len > 1:
    if upd_new_len > 30 and cpu_count >= 4:
        # dump nw array to /dev/shm
        nw_memmap_fn = os.path.join(temp_folder, 'nw_matrix.mmap')
        dump(inputnewseqmat, nw_memmap_fn)

        # delete from memory, do garbage collection
        del inputnewseqmat
        gc.collect()

        nwseqmat_memmap = load(nw_memmap_fn, mmap_mode='r')

        no_jobs = min(upd_new_len, cpu_count)
        batch_size = int(math.ceil(upd_new_len/float(no_jobs)))

        dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_square)(nwseqmat_memmap, i, batch_size, upd_new_len) for i in xrange(0,upd_new_len,batch_size))

        # put it together, fewer so concatenate
        dist_nw_nw = dist_arr[0]
        for np_arr in dist_arr[1:]:
            dist_nw_nw = np.concatenate((dist_nw_nw,np_arr), axis=1)

    else:
        dist_nw_nw = np.zeros(shape=(upd_new_len,upd_new_len), dtype=np.uint16)
        for i in range(upd_new_len):
            # lower triangle
            for j in range(0,i):
                    dist_nw_nw[i,j] = np.not_equal(inputnewseqmat[i,],inputnewseqmat[j,]).sum(0) - np.not_equal(inputnewseqmat[i,]!= 0, inputnewseqmat[j,]!= 0).sum(0)

        # del sequence matrix when done
        del inputnewseqmat
        gc.collect()

    # find clustered isolates in the spirit of hobohm1 using original indeces
    secondary, primary  = np.where(dist_nw_nw < config.THRESHOLD)
    pre = -1
    redundant = []
    for i, pr_index in enumerate(primary):
        if non_clustered[pr_index] not in new_clusters and secondary[i] > pr_index:
            if secondary[i] != pre or dist_nw_nw[secondary[i]][pr_index] < dist_new[(new_clusters[non_clustered[secondary[i]]], non_clustered[secondary[i]])]:
                new_clusters[non_clustered[secondary[i]]] = non_clustered[pr_index]
                dist_new[(non_clustered[pr_index], non_clustered[secondary[i]])] = dist_nw_nw[secondary[i]][pr_index]
                # current indexing
                redundant.append(secondary[i])
                pre = secondary[i]

    # TODO: DONE when reducing is done remove rows&cols from the dist_nw_nw matrix using current indexing
    if new_clusters:
        dist_nw_nw = np.delete(dist_nw_nw, redundant, axis=0)
        dist_nw_nw = np.delete(dist_nw_nw, redundant, axis=1)
        dist_nw_o = np.delete(dist_nw_o, redundant, axis=0)

        # original indexes that were kept
        reduced = [x for x in non_clustered if x not in new_clusters.keys()]
    else:
        reduced = non_clustered

else:
    pass

timing("# Clustered the new isolates.")

# if more than one old seq then calc distances
if slens[0] > 1:
    # get the dist between the old ones from pickle
    if not pristine and not re_calc and os.path.exists(dist_pic_filepath):
        with open(dist_pic_filepath, "rb") as picfile:
            matrix = pickle.load(picfile)
    else:
        dist_o_o = np.zeros(shape=(slens[0],slens[0]), dtype=np.uint16)

        for vol1 in range(no_vol):
            # load the volume and mask it
            hr_matrix_npy_vol1 = hr_matrix_npy.format(vol1+1)
            inputhrseqmat1 = np.load(hr_matrix_npy_vol1, mmap_mode=None, allow_pickle=True, fix_imports=True)
            #vol_len = np.shape(inputhrseqmat)[0]
            if m is not None:
                inputhrseqmat1 = inputhrseqmat1.T[m].T

            # dump hr array to /dev/shm
            dump(inputhrseqmat1, hr_memmap_fn.format(vol1+1))

            # delete from memory, do garbage collection
            del inputhrseqmat1
            gc.collect()

            # load as memmap
            hrseqmat_memmap_1 = load(hr_memmap_fn.format(vol1+1), mmap_mode='r')

            # calculate genetic distance between old and new isolates
            vol_len_1 = np.shape(hrseqmat_memmap_1)[0]
            no_jobs = min(vol_len_1, cpu_count)
            batch_size = int(math.ceil(vol_len_1/float(no_jobs)))

            for vol2 in range(vol1, no_vol):
                # load the volume and mask it
                hr_matrix_npy_vol2 = hr_matrix_npy.format(vol2+1)
                inputhrseqmat2 = np.load(hr_matrix_npy_vol2, mmap_mode=None, allow_pickle=True, fix_imports=True)
                #vol_len = np.shape(inputhrseqmat)[0]
                if m is not None:
                    inputhrseqmat2 = inputhrseqmat2.T[m].T

                # dump hr array to /dev/shm
                dump(inputhrseqmat2, hr_memmap_fn.format(vol2+1))

                # delete from memory, do garbage collection
                del inputhrseqmat2
                gc.collect()

                hrseqmat_memmap_2 = load(hr_memmap_fn.format(vol2+1), mmap_mode='r')
                vol_len_2 = np.shape(hrseqmat_memmap_2)[0]

                dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_lower)(hrseqmat_memmap_1, hrseqmat_memmap_2, i, batch_size, vol_len_2, vol_len_1) for i in xrange(0,vol_len_1,batch_size))
                # place it into big dist matrix
                o_x = 0
                for np_arr in dist_arr:
                    height, width = np.shape(np_arr)
                    #print(np_arr)
                    for x in range(width):
                        for y in range(height):
                            #print(np_arr[y,x])
                            row = vol2 * config.VOL_SIZE + y
                            col = vol1 * config.VOL_SIZE + o_x + x
                            if row > col:
                                dist_o_o[row,col] = np_arr[y,x]
                    o_x += width

                # cleanup
                del hrseqmat_memmap_2
                if vol2 != vol1:
                    os.unlink(hr_memmap_fn.format(vol2+1))

            del hrseqmat_memmap_1
            os.unlink(hr_memmap_fn.format(vol1+1))

        del dist_arr[:]
        matrix = dist_o_o.tolist()

else: # only one old seq, dist 0
    matrix = [[0]]

# print(matrix)
timing("# Calculated distances between old sequences.")

# combine all non-redundant into a distance matrix
seqnames = ["template"]
seqnames.extend([x.split(".")[0] for x in oldseqs[1:]])

if dist_nw_nw is not None:
    orig_rows = len(matrix)
    newm = dist_nw_nw.tolist()
    for i,nw_index in enumerate(reduced):
        seqnames.append(newseqs[nw_index].split(".")[0])
        new_row = dist_nw_o[i,].tolist()
        for r in xrange(orig_rows):
            #matrix[r].append(dist_nw_o[i][r])
            matrix[r].append(0)
        matrix.append(new_row)
    full_rows = len(matrix)
    for i in xrange(orig_rows, full_rows):
        matrix[i].extend(newm[i - orig_rows])

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

    with open(os.path.join(args.odir, "seqid2name.{}.pic".format(method)), "w") as pf:
        pickle.dump(seqid2name, pf)
    # dump matrix in pickle for pairwise
    #if not pristine:
    if not args.allcalled or slens[0] > config.PHOLD:
        # if order of files changes, then the old-old dists are re-calculated
        dist_pic_filepath += suffix
        with open(dist_pic_filepath, "wb") as picklingfile:
            pickle.dump(matrix, picklingfile)


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

    if not args.debug:
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

# TODO DONE add to the latest volume up to 500 and start a new one if needed
if np.shape(inputnewseqmat)[0]:
    # re-load the last volume
    vol_array = np.load(hr_matrix_npy.format(no_vol), mmap_mode=None, allow_pickle=True, fix_imports=True)
    hr_len = np.shape(vol_array)[0]
    nw_len = np.shape(inputnewseqmat)[0]
    fill = min(config.VOL_SIZE - hr_len, nw_len)
    new_vol = int(math.ceil((nw_len - fill)/config.VOL_SIZE))
    # add to/fill the volume
    vol_array = np.concatenate((vol_array, inputnewseqmat[:fill]), axis = 0)
    np.save(hr_matrix_npy.format(no_vol), vol_array, allow_pickle=True, fix_imports=True)
    # rest are saved into volume(s)
    for i in range(new_vol):
        np.save(hr_matrix_npy.format(no_vol + i + 1), inputnewseqmat[fill+i*config.VOL_SIZE:fill+(i+1)*config.VOL_SIZE], allow_pickle=True, fix_imports=True)

    del vol_array

    timing("# Updated the redundant sequences npy.")

# cleanup
del inputnewseqmat
gc.collect()

try:
    os.unlink(os.path.join(args.odir, "new-matrix.npy"))
    for i in range(no_vol):
        os.unlink(hr_memmap_fn.format(vol+1))
    os.rmdir(temp_folder)
except OSError:
    import shutil
    shutil.rmtree(temp_folder)


timing("# Distance calculation is finished.")
print("DONE", file=sys.stderr)
sys.exit(0)
