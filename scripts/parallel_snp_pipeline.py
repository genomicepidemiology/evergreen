#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os, time, math
import argparse
import subprocess
import shlex
import shutil
import glob
from distutils.spawn import find_executable
from random import shuffle
import multiprocessing
import sqlite3

J_LIMIT = 4
KT = "kmer_tax.py"
ADT = "assimpler_distance_trees.py"

parser = argparse.ArgumentParser(
    description='Parallel SNP pipeline')
parser.add_argument(
   '-b',
   dest="base",
   help='Base output directory, absolute path')
parser.add_argument(
    '-f',
    dest="isolates_file",
    default=None,
    help='File with sample names and paths to raw runs. Files for same isolate given comma separated.')
parser.add_argument(
    '-i',
    dest="collection_file",
    default=None,
    help='Collection file of raw reads')
parser.add_argument(
    '-a',
    dest="allcalled",
    action="store_true",
    help='Consider all positions in column for Ns')
parser.add_argument(
    '-p',
    dest="pairwise",
    action="store_true",
    help='Consider the sequences pairwise for Ns')
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
    '-q',
    dest="quiet",
    action="store_true",
    help='Quiet')
args = parser.parse_args()

def jobstart(command):
    #job = print(command)
    cmd = shlex.split(command)
    job = subprocess.call(cmd)
    return job

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)

def timing(message):
    if not args.quiet:
        t1 = time.time()
        print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def callwhere(cmd):
    if args.debug:
        print(cmd)

def adt_call_p(tmpl_file):
    print(tmpl_file)

def adt_call(tmpl_file):

    adt_cmd = "{0} -b {1} -k".format(ADT, args.base)
    if args.allcalled and not args.pairwise:
        adt_cmd += ' -a'
    if args.likelihood:
        adt_cmd += ' -L'
    if args.distance:
        adt_cmd += ' -D'
    if args.debug:
        adt_cmd += ' -d'
    if args.quiet:
        adt_cmd += ' -q'

    adt_cmd += " -i {0}".format(os.path.join(wdir, tmpl_file))

    print("Subprocess:", adt_cmd)
    cmd = shlex.split(adt_cmd)
    subprocess.call(cmd)
    return

def fit_folder(folder):
    """Avoid having more than 10000 files in one folder."""
    full_folder = os.path.join(bdir, folder)
    i = 1
    while len(os.listdir(full_folder)) > 10000:
        i += 1
        full_folder = os.path.join(bdir, "{0}-{1}".format(folder, i))
        if not os.path.exists(full_folder):
            os.mkdir(full_folder)
    return full_folder

## Main
t0 = time.time()
pid = os.getpid()
todaysdate = time.strftime("%d%m%Y")
# unbuffered output
os.environ['PYTHONUNBUFFERED'] = '1'

# check for binaries
if find_executable("kma") is None:
    exiting("KMA not in path")
if find_executable("iqtree") is None:
    exiting("Iqtree not in path")
if find_executable("neighbor") is None:
    exiting("Neighbor not in path")

debug_opt = ""
if args.debug:
    debug_opt = "-d"

# Check base dir
bdir = os.path.realpath(args.base)
if not os.path.isdir(bdir):
    exiting("Base path is required.")

# put scripts folder into path
if find_executable(KT) is None or find_executable(ADT) is None:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] = "{}:{}".format(script_dir,os.environ['PATH'])

wdir = fit_folder("output")

if not args.allcalled and not args.pairwise:
    exiting("Either -a or -p option is needed!")

# db paths
pattern = os.path.join(bdir, "hr_database/current/*.name*")
kma_db = glob.glob(pattern)[0].split(".name")[0]
db_paths = "-db {0} -f_db {1} -fa_db {2}".format(kma_db,
  os.path.join(bdir, "hr_database/current/bacteria.folder.pic"),
  os.path.join(bdir, "hr_database/current/bacteria.fsa_name.pic")
)
MAIN_SQL_DB = os.path.join(bdir, "results_db/evergreen.db")

if args.isolates_file is not None:

    if not os.path.exists(args.isolates_file):
        exiting("List of raw reads required.")

    opt_cmd = "" + debug_opt
    kt_cmd = "{0} -l {1} -o {2} -wdir {2} {3} {4}".format(KT, args.isolates_file, wdir, db_paths, opt_cmd)

elif args.collection_file is not None:
    collection = args.collection_file
    if not os.path.exists(collection):
        exiting("Collections file is unavailable, {0}.".format(collection))

    #
    # /data/evergreen/scripts/kmer_tax
    # -h, --help            show this help message and exit
    # -i COLLECTION_FILE    Input file of collected runs
    # -f RUNS_FILENAME      Comma separated paths to raw runs
    # -o ODIR               Output directory
    # -db DATABASE          Database, absolute path
    # -f_db FOLDER_DATABASE
    #                       Acc number to folder pickle
    # -fa_db FSA_DATABASE   Acc number to filename pickle
    # -wdir WDIR            Working directory, not required
    # -d, --debug           Debug
    # -q                    Quiet

    opt_cmd = "" + debug_opt
    kt_cmd = "{0} -i {1} -o {2} -wdir {2} {3} {4}".format(KT, collection, wdir, db_paths, opt_cmd)

else:
    exiting("Input files not given.")

callwhere(kt_cmd)
retcode = jobstart(kt_cmd)
if retcode:
    exiting("Kmer_tax exited with error.")

start_dir = os.environ.get('PBS_O_WORKDIR')
os.chdir(wdir)
prefix = "kmer_{0}_{1}".format(todaysdate, pid)
tmpl_files = glob.glob("{}_*.tmpl".format(prefix))
job_no = len(tmpl_files)
if not job_no:
    exiting("Warning: no templates were found for prefix: {0}\n{1}".format(prefix, collection))

if start_dir is not None:
    os.chdir(start_dir)

J_LIMIT = min(job_no, J_LIMIT)

if __name__ == '__main__':

    p = multiprocessing.Pool(J_LIMIT)
    p.map(adt_call, tmpl_files)

# update the db to see what was included in the trees
conn = sqlite3.connect(MAIN_SQL_DB)
conn.execute("PRAGMA foreign_keys = 1")
cur = conn.cursor()

cur.execute('''SELECT count(*) from sqlite_master where type='table' and name='templates';''')
if cur.fetchone()[0] == 1:
    included_update = []
    cur.execute('''SELECT sra_id from runs where included=2;''')
    svar = cur.fetchall()
    for row in svar:
        # if not on at least one tree then not included
        cur.execute('''SELECT sum(qc_pass) from templates where sra_id=?''', (row[0],))
        templ_no = cur.fetchone()[0]
        if templ_no is not None and templ_no != 0:
            included_update.append((1, row[0]))
        else:
            included_update.append((0, row[0]))

    if included_update:
        cur.executemany('''UPDATE runs SET included=? WHERE sra_id=?''', included_update)
        conn.commit()

# get a list of trees made in this run
cur.execute('''SELECT count(*) from sqlite_master where type='table' and name='trees';''')
if cur.fetchone()[0] == 1:
    cur.execute('''SELECT template,ctime,nw_path from trees where ctime > datetime(?, 'unixepoch');''', (t0,))
    print("\n=== Inferred trees in this run ===\n")
    print("template", "time", "path", sep="\t")
    for row in cur.fetchall():
        print("\t".join(row))


print("DONE", file=sys.stderr)
sys.exit(0)
