#!/usr/bin/env python2.7
from __future__ import print_function
import sys, os, time
import argparse
import shutil
import shlex
import subprocess
import pickle
import sqlite3
from operator import itemgetter
from ete3 import Tree

base_path = os.path.dirname(sys.argv[0]).rsplit("/",1)[0]
IQTREE = os.path.join(base_path, "scripts/iqtree")
DTREE = os.path.join(base_path, "scripts/neighbor")
MAIN_SQL_DB = os.path.join(base_path, "results_db/evergreen.db")



parser = argparse.ArgumentParser(
    description='Wrapper for the tree infering program for the Evergreen pipeline')
parser.add_argument(
   '-b',
   dest="bdir",
   help='Base for resutls_db/template directory, absolute path')
parser.add_argument(
    '-m',
    dest="distmat",
    default=None,
    help='Distance matrix for NJ')
parser.add_argument(
    '-f',
    dest="filelist",
    default=None,
    help='List of sequences')
parser.add_argument(
    '-d',
    dest="distance",
    action="store_true",
    help='Distance based phylogenic tree')
parser.add_argument(
    '-l',
    dest="likelihood",
    action="store_true",
    help='Maximum likelihood based phylogenic tree')
parser.add_argument(
    '-a',
    dest="allcalled",
    action="store_true",
    help='Mode is "all"')
parser.add_argument(
    '-t', '--debug',
    dest="debug",
    action="store_true",
    help='Debug: use .t suffix, print trees to screen')
parser.add_argument(
    '-k', '--keep',
    dest="keep",
    action="store_true",
    help='Keep temporary subdirectory')
parser.add_argument(
    '-q',
    dest="quiet",
    action="store_true",
    help='Quiet')
args = parser.parse_args()

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

def input_validation(filename):
    if filename is None:
        exiting("Input file not given.")
    if not os.path.exists(filename):
        exiting("Input file not given.")
    else:
        return os.path.abspath(filename)

def change2subdir(mode):
    if mode in ["all", "pw"]:
        subdir = os.path.join(bdir, "{0}_{1}_{2}{3}".format(mode, method, ctime, suffix))
    else:
        subdir = os.path.join(bdir, "{0}_{1}{2}".format(method, ctime, suffix))

    subdir = subdir.replace(".", "_")
    try:
        os.mkdir(subdir)
        os.chdir(subdir)
    except:
        exiting("Couldn't make {0}".format(subdir))
    return subdir

def decorate_tree(treepath, mode):
    try:
        tree = Tree(treepath)
    except NewickError as e:
        exiting("Couldn't open {0}".format(treepath))

    # go through the tree and check for clusters in db
    full_tree = tree.copy()
    for node in tree.traverse():
        if node.name:
            cur.execute('''SELECT sra_id from sequences where repr_id=?''', (node.name,))
            records = cur.fetchall()
            if records:
                homNode = full_tree.search_nodes(name=node.name)[0]
                oriNode = homNode.add_child(name=node.name, dist=0)
                for rec in records:
                    isoNode = homNode.add_child(name="*{}".format(rec[0]), dist=0)

    outpath = os.path.join(bdir, "_tree", "{0}_{1}_{2}_{3}.newick{4}".format(template_name, mode, method, ctime, suffix))
    full_tree.write(format=0, outfile=outpath)
    if args.debug:
        print(full_tree)
    return

def decode_nw_file(nw_file, mode):
    # load the pickle with id: seqname pairs
    seqid2name = {}
    with open(os.path.join(args.bdir, "seqid2name.{}.pic".format(mode)), "r") as pf:
        seqid2name = pickle.load(pf)

    newick = []
    with open(nw_file, "r") as nw_f:
        for line in nw_f:
            newick.append(line)

    newick_str = "".join(newick)
    for key in seqid2name.keys():
        newick_str = newick_str.replace(key, seqid2name.get(key))

    with open(nw_file, "w") as nw_f:
        nw_f.write(newick_str)

    return
def print2msa(msafile, fsafilename, seqname):
    print(">{0}".format(seqname), file=msafile)
    entries = zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fsafilename)])
    print("".join(list(entries[0])), file=msafile)
    return

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

## Main ##
# Start time to keep track of progress
t0 = time.time()
ctime = int(t0)
etta = 0.001


# Valid base dir
if args.bdir is not None:
    bdir = os.path.realpath(args.bdir)
if args.bdir is None or not os.path.isdir(bdir):
    exiting("Base path is required.")
try:
    oldpath = os.getcwd()
    os.chdir(bdir)
except:
    exiting("Couldn't change to {0}".format(bdir))

suffix = ""
if args.debug:
    suffix = ".t"

mode = 'pw'
if args.allcalled:
    mode = 'all'

treefilename = None

# open database
db_path = os.path.join(bdir, "isolates.{0}.db{1}".format(mode, suffix))
#print(db_path)
conn = sqlite3.connect(db_path)
conn.execute("PRAGMA foreign_keys = 1")
cur = conn.cursor()

if args.distance:
    method = "dist"
    matpath = input_validation(args.distmat)
    wdir = change2subdir(mode)
    matpath = os.path.relpath(matpath)
    treefilename = "outtree"

    inputstr = "{0}\nY\n".format(matpath)
    proc = subprocess.Popen(DTREE, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.stdin.write(inputstr)
    # wait for it. it gets stuck in simple wait() if fails
    (stdoutdata, stderrdata) = proc.communicate()
    if proc.returncode:
        exiting("Neighbor program failed.")

    with open("output", "w") as ofile:
        print(stdoutdata, file=ofile)

    decode_nw_file(os.path.join(wdir, "outtree"), mode)

    timing("# Distance based tree constructed.")

elif args.likelihood:
    method = "ml"
    matpath = input_validation(args.distmat)
    if args.filelist == '-':
        templ = os.path.basename(bdir)
        isolates = ["pt_{}.fa".format(templ)]
        cur.execute('''SELECT sra_id from sequences where repr_id is NULL and sra_id!='template';''')
        rows = cur.fetchall()
        for r in rows:
            isolates.append("{}.fa".format(r[0]))
    else:
        filepath = input_validation(args.filelist)

        try:
            pathfile = open(filepath, "r")
            isolates = pathfile.readlines()
            pathfile.close()
        except IOError as e:
            exiting("File couldn't be opened, {0}.".format(e))

    wdir = change2subdir(mode)
    matpath = os.path.relpath(matpath)
    treefilename = "sequences.fa.treefile"

    # get the NJ tree made
    inputstr = "{0}\nY\n".format(matpath)
    proc = subprocess.Popen(DTREE, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.stdin.write(inputstr)
    # wait for it. it gets stuck in simple wait() if fails
    (stdoutdata, stderrdata) = proc.communicate()
    if proc.returncode:
        exiting("Neighbor program failed.")

    with open("output", "w") as ofile:
        print(stdoutdata, file=ofile)

    decode_nw_file(os.path.join(wdir, "outtree"), mode)

    # make the maxL tree
    try:
        msafile = open("sequences.fa", "w")
    except IOError as e:
        exiting("File couldn't be opened, {0}.".format(e))

    # 1st seq is template -> change name!
    print2msa(msafile, os.path.join(bdir, isolates[0].strip()), "template")
    seq_count = 1
    for iso in isolates[1:]:
        longpath = os.path.join(bdir, iso.strip())
        seqname = iso.split(".")[0]
        if os.path.isfile(longpath):
            print2msa(msafile, longpath, seqname)
            seq_count += 1
            # high five if nothing is wrong
    msafile.close()

    # construct iqtree cmd
    # /data/evergreen/scripts/iqtree-omp
    # -s Escherichia_coli_K_12_substr__MG1655_uid57779_99.fa -st DNA -m GTR
    # -nt AUTO -mem 14Gb -bb 1000
    # dont bootstrap if -fast (not compatible)
    cmd = "{0} -s sequences.fa -st DNA -m GTR -nt AUTO -mem 64Gb -nstop 50 -t outtree".format(IQTREE)
    with open("sequences.fa.out", "w") as ofile:
        exitcode = subprocess.call(shlex.split(cmd), stdout=ofile)
    if exitcode:
        exiting("Iqtree failed.")

    timing("# Maximum likelihood based tree constructed.")

else:
    # neither option was chosen
    exiting("Tree method was not chosen.")

template_name = os.path.basename(bdir)
fulltree_file = os.path.join(bdir, "_tree", "{0}_{1}_{2}_{3}.newick{4}".format(template_name, mode, method, ctime, suffix))
cur.execute('''SELECT count(*) from sequences where repr_id is not Null;''')
numb = cur.fetchone()
if numb and numb[0]:
    decorate_tree(os.path.join(wdir, treefilename), mode)
else:
    try:
        shutil.copy(
        os.path.join(wdir, treefilename),
        fulltree_file
        )
    except:
        exiting("Tree couldn't be copied.")

if args.likelihood:
    try:
        os.unlink(os.path.join(wdir, "sequences.fa"))
    except OSError:
        pass

conn.close()
timing("# Tree decorated.")

# make a tree with metadata
if not args.debug:
    conn = sqlite3.connect(MAIN_SQL_DB)
    conn.execute("PRAGMA foreign_keys = 1")
    conn.commit()
    cur = conn.cursor()

    # save tree to db
    cur.execute('''CREATE TABLE IF NOT EXISTS trees
        (db_id INTEGER PRIMARY KEY,
        template TEXT,
        method TEXT,
        mode TEXT,
        ctime DATETIME DEFAULT CURRENT_TIMESTAMP,
        nw_path TEXT,
        UNIQUE(template, method, mode)
        );''')
    conn.commit()

    cur.execute('''INSERT OR REPLACE INTO trees (template,method,mode,nw_path) VALUES (?,?,?,?)''', (template_name, method, mode, fulltree_file))
    conn.commit()
    conn.close()

os.chdir(oldpath)
if not args.keep:
    try:
        shutil.rmtree(wdir)
    except:
        pass

print("Done.", file=sys.stderr)
sys.exit(0)
