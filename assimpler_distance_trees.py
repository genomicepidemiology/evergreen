#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os, time
import argparse
import shutil
import subprocess
import shlex
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import sqlite3

J_LIMIT = int(cpu_count() / 4)

"""
Takes: one file with the collected info about the isolates (from kmer_tax)

Change:
- for COMPARE: delete ml tree from db if > 650 seqs and thus no new ml tree is done
"""

parser = argparse.ArgumentParser(
    description='Wrapper for Assimpler and distance calculations for the Evergreen pipeline')
parser.add_argument(
   '-b',
   dest="base",
   default="/data/evergreen",
   help='Base output directory, absolute path')
parser.add_argument(
    '-i',
    dest="collection_file",
    default=None,
    help='Input file of isolates with the same template')
parser.add_argument(
    '-a',
    dest="allcalled",
    action="store_true",
    help='Consider all positions in column for Ns')
parser.add_argument(
    '-D',
    dest="distance",
    action="store_true",
    help='Distance based phylogenic tree')
parser.add_argument(
    '-L',
    dest="likelihood",
    action="store_true",
    help='Maximum likelihood based phylogenic tree')
parser.add_argument(
    '-d',
    dest="debug",
    action="store_true",
    help='Debug: use .t suffixes')
parser.add_argument(
    '-k',
    dest="keep",
    action="store_true",
    help='Keep temp tree folders')
parser.add_argument(
    '-q',
    dest="quiet",
    action="store_true",
    help='Quiet')
args = parser.parse_args()

def jobstart(command):
    #job = print(command)
    cmd = shlex.split(command)
    job = subprocess.call(cmd, stdout=logfile, stderr=logfile)
    return job

def jobstart_silent(command):
    #job = print(command)
    cmd = shlex.split(command)
    try:
        devnull = open(os.devnull, 'w')
    except IOError as e:
        exiting("Devnull won't open.")
    job = subprocess.call(cmd, stdout=devnull, stderr=devnull)
    devnull.close()
    return job

def exiting(message):
    print(message, file=logfile)
    print("FAIL", file=logfile)
    sys.exit(1)

def timing(message):
    if not args.quiet:
        t1 = time.time()
        print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=logfile)
        # flush it
        logfile.flush()
        os.fsync(logfile.fileno())
    return

def parse_input():
    ''' Parse input file and return a list of lists with info for each acc.'''

    # SRR5297990    /data/evergreen/data/SRR5297990_1.fastq.gz /data/evergreen/data/SRR5297990_2.fastq.gz	Salmonella_enterica_serovar_Enteritidis_P125109_uid59247

    try:
        ifile = open(args.collection_file, "r")
    except IOError as e:
        exiting("Can't open file {0}".format(e))

    accessions = []
    for line in ifile:
        accessions.append(line.strip().split("\t"))
    ifile.close()

    return accessions

## Main
t0 = time.time()

# Check base dir and db
bdir = os.path.realpath(args.base)
if not os.path.isdir(bdir):
    exiting("Base path is required.")

MAP = os.path.join(bdir, "scripts/assimpler_pipeline.py")
DIST = os.path.join(bdir, "scripts/distance_batch.py")
PTREE = os.path.join(bdir, "scripts/tree_cobbler.py")
MAIN_SQL_DB = os.path.join(bdir, "results_db/evergreen.db")

# Process collection file
inputs = []
logfile = None
if args.collection_file is not None:
    logfilename = os.path.join(bdir, "logs/{}.log".format(os.path.split(args.collection_file)[-1].split(".")[0]))
    logfile = open(logfilename, "w")
    inputs = parse_input()
else:
    exiting("Please specify input file.")

# target directory from the template in the kmerfinder results
template = inputs[0][2]
wdir = os.path.join(bdir, "results_db", template)

if args.allcalled:
    hrfilename = os.path.join(wdir, "non-redundant.all.lst")
    db_path = os.path.join(wdir, "isolates.all.db")
else:
    hrfilename = os.path.join(wdir, "non-redundant.pw.lst")
    db_path = os.path.join(wdir, "isolates.pw.db")

# db management
suffix = ""
if args.debug and os.path.exists(db_path):
    suffix = ".t"
    try:
        shutil.copy(db_path, db_path + suffix)
        db_path = db_path + suffix
    except OSError:
        exiting("Couldnt copy the temp database.")

iso_conn = sqlite3.connect(db_path)
iso_conn.execute("PRAGMA foreign_keys = 1")
iso_cur = iso_conn.cursor()

iso_cur.execute('''CREATE TABLE IF NOT EXISTS sequences
    (sra_id TEXT PRIMARY KEY,
    repr_id TEXT DEFAULT NULL,
    distance INTEGER DEFAULT NULL)''')

# if table is newly created, then the reference sequence should be placed in it
iso_cur.execute('''SELECT count(*) FROM sequences''')
if not iso_cur.fetchone()[0]:
    iso_cur.execute('''INSERT INTO sequences (sra_id) VALUES (?)''', ('template',))
    iso_conn.commit()

# only the template is in the folder
if not os.path.exists(hrfilename):
    with open(hrfilename, "w") as hrfile:
        print("pt_{0}.fa".format(template), end= "", file=hrfile)

timing("# Preparation done for: {} ".format(template))

###
# /data/evergreen/scripts/assimpler_pipeline.py
# -t template_file,
# -i ' '.join(inputFiles),
# -c assimpler_consensus_file,
# -n SRAid
# -z 1.96
# -k 17
# -l 50
###

mapper_cmds = []
cons_files = []
template_file = os.path.join(wdir, "pt_{0}.fa".format(template))
for acc in inputs:
    consensus_file = os.path.join(wdir, "{0}.fa".format(acc[0]))
    cons_files.append("{0}.fa".format(acc[0]))
    # only map if no consensus file exists or size 0
    if not os.path.exists(consensus_file) or not os.path.getsize(consensus_file) > 0:
        cmd = "{0} -i {1} -t {2} -c {3} -n {4} -z 1.96 -k 17 -l 50".format(
              MAP, acc[1], template_file, consensus_file, acc[0]
        )
        mapper_cmds.append(cmd)

if args.debug:
    if mapper_cmds:
        print("# 1st mapping command:", mapper_cmds[0])
    print("# number of mapping commands:", len(mapper_cmds))

# Start assimpler mapper
if mapper_cmds:
    jobs = Parallel(n_jobs=J_LIMIT)(delayed(jobstart_silent)(cmd) for cmd in mapper_cmds)

    # if mapping unsuccessful the file just won't load in DIST
    if sum(jobs) != 0:
        print("Warning: A subprocess was unsuccessful.", file=sys.stderr)

timing("# Mapping done.")

# write cons files list to cf.txt
newfilename = "{0}/new_isolates.lst".format(wdir)
with open(newfilename, "w") as ofile:
    print("\n".join(cons_files), file=ofile)

# TODO DONE insert or ignore those isolates with repr_id = N to the db that have os.path.getsize(consensus_file) > 0
new_cons = []
templ_update = []
runs_update = []
for filenm in cons_files:
    full_cons_path = os.path.join(wdir, filenm)
    if os.path.exists(full_cons_path) and os.path.getsize(full_cons_path) > 0:
        new_cons.append((filenm[:-3], 'N'))
    else:
        templ_update.append((0, filenm[:-3], template))
        runs_update.append((2, filenm[:-3]))

if new_cons:
    iso_cur.executemany('''INSERT OR IGNORE INTO sequences (sra_id, repr_id) VALUES (?,?);''', new_cons)
    iso_conn.commit()

if templ_update:
    # open database
    # MAIN
    conn = sqlite3.connect(MAIN_SQL_DB)
    conn.execute("PRAGMA foreign_keys = 1")
    conn.commit()
    cur = conn.cursor()
    # reach back and update runs table in the main DB to 2 as non-included
    try:
        cur.executemany('''UPDATE templates SET qc_pass=? WHERE sra_id=? and template=?''', templ_update)
        cur.executemany('''UPDATE runs SET included=? WHERE sra_id=?''', runs_update)
        conn.commit()
    except sqlite3.Error:
        print("Warning: SQL update failed.", file=sys.stderr)
    conn.close()

# TODO DONE change to '-' instead of filelist

# call to the distance and hr reduction script
# /data/evergreen/scripts/distance_batch \
# -hr /data/evergreen/results_db/Salmonella_enterica_serovar_Enteritidis_P125109_uid59247/non-redundant.all.lst \
# -n /data/evergreen/results_db/Salmonella_enterica_serovar_Enteritidis_P125109_uid59247/new_isolates.lst \
# -o /data/evergreen/results_db/Salmonella_enterica_serovar_Enteritidis_P125109_uid59247 \
# -m /data/evergreen/results_db/Salmonella_enterica_serovar_Enteritidis_P125109_uid59247/dist.mat \
# -a --vcf -d
if args.allcalled:
    cmd = "{0} -hr {1} -n - -o {3} -a".format(DIST, "-", newfilename, wdir)
else:
    cmd = "{0} -hr {1} -n - -o {3}".format(DIST, "-", newfilename, wdir)

if args.debug:
    cmd += " -d"
    print("# distance command: ", cmd)

# call shlex.split(cmd)
ecode = jobstart(cmd)
if ecode:
    exiting("Error: distance calculation is unsuccessful.")

# check if there is a _trees folder and make one if not
if not os.path.exists(os.path.join(wdir, "_tree")):
    try:
        os.mkdir(os.path.join(wdir, "_tree"))
    except:
        exiting("Error: can't create _tree directory.")

timing("# Distance calculation is done.")

# use .t file if debug
mode = 'pw'
if args.allcalled:
    hrfilename = os.path.join(wdir, "non-redundant.all.lst{0}".format(suffix))
    matfilename = os.path.join(wdir, "dist.all.mat{0}".format(suffix))
    mode = 'all'
else:
    hrfilename = os.path.join(wdir, "non-redundant.pw.lst{0}".format(suffix))
    matfilename = os.path.join(wdir, "dist.pw.mat{0}".format(suffix))



# only if more than 3 seqs in non-redundant.*.lst !
# TODO DONE reimplement this as a db check
iso_cur.execute('''SELECT count(*) from sequences where repr_id is Null;''')
non_redundant_no = iso_cur.fetchone()[0]
iso_conn.close()
if non_redundant_no > 2:
    add_opt = ""
    if args.debug:
        add_opt = "-t"
    if args.keep:
        add_opt += " -k"

    if args.allcalled:
        add_opt += " -a"


    # allow both methods to be run: fork a distance and a ML tree maker
    procs = []
    if args.distance:
        cmd = "{0} -b {1} -m {2} {3} -d".format(PTREE, wdir, matfilename, add_opt)
        if args.debug:
            print("# Tree command: ", cmd)
        procs.append(subprocess.Popen(shlex.split(cmd), stdout=logfile, stderr=logfile))

    if args.likelihood and non_redundant_no > 3 and non_redundant_no <= 300:
        cmd = "{0} -b {1} -f {2} {3} -l -m {4}".format(PTREE, wdir, "-", add_opt, matfilename)
        if args.debug:
            print("# Tree command: ", cmd)
        procs.append(subprocess.Popen(shlex.split(cmd), stdout=logfile, stderr=logfile))
    elif args.likelihood and non_redundant_no > 300:
        # too many isolates for ml, thus the ml tree in the db is deleted
        conn = sqlite3.connect(MAIN_SQL_DB)
        conn.execute("PRAGMA foreign_keys = 1")
        conn.commit()
        cur = conn.cursor()
        # reach back and update runs table in the main DB to 2 as non-included
        try:
            cur.execute('''DELETE or IGNORE from trees where template=? and method='ml' and mode=?;''', (template, mode))
            conn.commit()
        except sqlite3.Error:
            print("Warning: SQL delete failed.", file=sys.stderr)
        conn.close()

    if sum([proc.wait() for proc in procs]) != 0:
        exiting("Error: tree inference is unsuccessful.")

    timing("# Tree inference is done.")
else:
    timing("# Tree inference not possible.")

print("DONE", file=logfile)
logfile.close()
sys.exit(0)
