#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os, time, math
import argparse
import subprocess
import shlex
import glob
import shutil
from random import shuffle
import multiprocessing
import sqlite3
try:
    import cPickle as pickle
except ImportError:
    import pickle

'''
Debug option is removed. Production version.
'''

J_LIMIT = 4
KT = "scripts/kmer_tax.py"
ADT = "scripts/assimpler_distance_trees.py"

parser = argparse.ArgumentParser(
    description='Parallel SNP pipeline')
parser.add_argument(
   '-b',
   dest="base",
   help='Base (install) directory, absolute path')
parser.add_argument(
    '-f',
    dest="isolates_file",
    default=None,
    help='File with sample names and paths to raw reads. Files for same isolate given comma separated.')
parser.add_argument(
    '-i',
    dest="collection_file",
    default=None,
    help='Collection file of raw reads')
parser.add_argument(
    '-a',
    dest="allcalled",
    action="store_true",
    help='Consider all positions in column for Ns, or')
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
    '-E',
    dest="ena",
    action="store_true",
    help='ENA compatible output (tsv + tar.gz)')
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

def adt_call(tmpl_file):

    # don't keep the tree temp dir
    # adt_cmd = "{0} -b {1} -k".format(os.path.join(args.base, ADT), args.base)
    adt_cmd = "{0} -b {1}".format(os.path.join(args.base, ADT), args.base)
    if args.allcalled and not args.pairwise:
        adt_cmd += ' -a'
    if args.likelihood:
        adt_cmd += ' -L'
    if args.distance:
        adt_cmd += ' -D'
    if args.quiet:
        adt_cmd += ' -q'

    adt_cmd += " -i {0}".format(os.path.join(wdir, tmpl_file))

    # print("Subprocess:", adt_cmd)
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

def isolates_to_tsv(template, outfile, newicks = "", treetime = "", matrix_filename = ""):
    iso_conn = sqlite3.connect(db_path.format(template))
    iso_cur = iso_conn.cursor()

    # get non-redundant isolates
    iso_cur.execute('''SELECT sra_id from sequences where repr_id is NULL;''')
    nonred_rows = iso_cur.fetchall()
    for nonred in nonred_rows:
        repr_iso = nonred[0]
        cl_dist = ""
        if repr_iso:
            # get cluster for that non-redundant
            iso_cur.execute('''SELECT sra_id,distance from sequences WHERE repr_id=? ORDER BY sra_id;''', (repr_iso,))
            clust_rows = iso_cur.fetchall()
            # there migth not be one
            if clust_rows:
                if repr_iso != "template":
                    # print repr_iso as the cluster rep for the cluster rep for easier filtering
                    print(repr_iso, repr_iso, "", template, treetime, newicks, matrix_filename, sep="\t", file=outfile)
                for iso in clust_rows:
                    print(iso[0], repr_iso, iso[1], template, treetime, newicks, matrix_filename, sep="\t", file=outfile)
            else:
                if repr_iso != "template":
                    # singletons are not clusters
                    print(repr_iso, "None", "", template, treetime, newicks, matrix_filename, sep="\t", file=outfile)

    iso_conn.close()

def decode_dist_matrix(seq2name_file, dist_mat_file):
    # load the pickle with id: seqname pairs
    seqid2name = {}
    with open(seq2name_file, "r") as pf:
        seqid2name = pickle.load(pf)

    mat_cont = []
    with open(dist_mat_file, "r") as mat_f:
        for line in mat_f:
            first_word = line.split(" ")[0]
            try:
                mat_cont.append(line.replace(first_word, seqid2name[first_word]))
            except KeyError:
                mat_cont.append(line)

    return "".join(mat_cont)

## Main
t0 = time.time()
pid = os.getpid()
todaysdate = time.strftime("%d%m%Y")
# unbuffered output
os.environ['PYTHONUNBUFFERED'] = '1'

# Check base dir
bdir = os.path.realpath(args.base)
if not os.path.isdir(bdir):
    exiting("Base path is required.")

wdir = fit_folder("output")

if not args.allcalled and not args.pairwise:
    exiting("Either -a or -p option is needed!")

# db paths
db_paths = "-db {0} -f_db {1} -fa_db {2}".format(os.path.join(bdir, "hr_database/current/bacterial_compl_genomes_hq99_k13_ATG"),
  os.path.join(bdir, "hr_database/current/bacteria.folder.pic"),
  os.path.join(bdir, "hr_database/current/bacteria.fsa_name.pic")
)
MAIN_SQL_DB = os.path.join(bdir, "results_db/evergreen.db")

if args.isolates_file is not None:

    if not os.path.exists(args.isolates_file):
        exiting("List of raw reads required.")

    opt_cmd = ""
    kt_cmd = "{0} -l {1} -o {2} -wdir {2} {3} {4}".format(os.path.join(bdir, KT), args.isolates_file, wdir, db_paths, opt_cmd)

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

    opt_cmd = ""
    kt_cmd = "{0} -i {1} -o {2} -wdir {2} {3} {4}".format(os.path.join(bdir, KT), collection, wdir, db_paths, opt_cmd)

else:
    exiting("Input files not given.")


retcode = jobstart(kt_cmd)
if retcode:
    exiting("Kmer_tax exited with error.")

start_dir = bdir
if os.environ.get('PBS_O_WORKDIR') is not None:
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
    nonred = cur.fetchall()
    for row in nonred:
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


if args.ena:

    mode = "pw"
    if args.allcalled:
        mode = "all"
    db_path = os.path.join(bdir, "results_db", "{}", "isolates.{}.db".format(mode))
    distmat_path = os.path.join(bdir, "results_db", "{}", "dist.{}.mat".format(mode))
    seq2name_path = os.path.join(bdir, "results_db", "{}", "seqid2name.{}.pic".format(mode))

    run_output_path = os.path.join(bdir, "run_{0}_{1}".format(todaysdate, pid))
    try:
        os.mkdir(run_output_path)
        tsv_fp = open(os.path.join(bdir, "run_{0}_{1}.tsv".format(todaysdate, pid)), "w")
    except OSError as e:
        exiting(str(e))

    # seq_id   ctime    template    qc_pass    repr_id
    print("Seq_ID","Repr_ID", "Distance", "Template","Time","Tree(s)", "Matrix", sep="\t", file=tsv_fp)
    # get templates with trees of given mode as they could differ
    cur.execute('''SELECT template FROM trees WHERE mode=? GROUP BY template ORDER BY MAX(ctime) DESC;''', (mode,))
    templates = cur.fetchall()
    # get isolates in the trees
    if templates is not None:
        for row in templates:
            # tree info
            templ = row[0]
            cur.execute('''SELECT ctime,nw_path from trees WHERE template=? and mode=? ORDER BY ctime DESC;''', (templ,mode))
            ctime = ""
            trees = []
            tree_rows = cur.fetchall()
            if tree_rows is not None:
                for tr in tree_rows:
                    #  ctime then tree path
                    if not ctime:
                        ctime = tr[0].split()[0]
                    # append just the filename
                    trees.append(os.path.basename(tr[1]))
                    # copy tree to be zipped
                    try:
                        shutil.copy(
                        tr[1],
                        os.path.join(run_output_path, trees[-1])
                        )
                    except:
                        exiting("Tree couldn't be copied.")

            # TODO DONE decode and copy matrix
            zipped_matrix_filename = ""
            if os.path.exists(distmat_path.format(templ)) and os.path.exists(seq2name_path.format(templ)):
                zipped_matrix_filename = "{0}_{1}.mat".format(templ, mode)
                with open(os.path.join(run_output_path, zipped_matrix_filename), "w") as mat_p:
                    print(decode_dist_matrix(seq2name_path.format(templ), distmat_path.format(templ)), file=mat_p)

            isolates_to_tsv(template = templ, outfile = tsv_fp, newicks = " ".join(trees), treetime = ctime,  matrix_filename = zipped_matrix_filename)

    # templates where < 3 isolates thus no trees
    cur.execute('''SELECT distinct(template) from templates WHERE template NOT IN (SELECT template FROM trees WHERE mode=? GROUP BY template);''', (mode,))
    minor_templates = cur.fetchall()
    if minor_templates is not None:
        for row in minor_templates:
            templ = row[0]

            zipped_matrix_filename = ""
            if os.path.exists(distmat_path.format(templ)) and os.path.exists(seq2name_path.format(templ)):
                zipped_matrix_filename = "{0}_{1}.mat".format(templ, mode)
                with open(os.path.join(run_output_path, zipped_matrix_filename), "w") as mat_p:
                    print(decode_dist_matrix(seq2name_path.format(templ), distmat_path.format(templ)), file=mat_p)

            if os.path.exists(db_path.format(templ)) and os.path.getsize(db_path.format(templ)):
                isolates_to_tsv(template = templ, outfile = tsv_fp,  matrix_filename = zipped_matrix_filename)

    # non-included isolates
    cur.execute('''SELECT sra_id FROM runs WHERE included=0 ORDER BY sra_id ASC;''')
    non_inc = cur.fetchall()
    if non_inc is not None:
        for iso in non_inc:
            print(iso[0], "", "", "", "", "",  "", sep="\t", file=tsv_fp)

    conn.close()
    tsv_fp.close()

    # zip the output dir
    os.chdir(bdir)
    exit_code = subprocess.call(shlex.split("tar -zcf {0}.tar.gz {0}/".format(os.path.basename(run_output_path))))
    if exit_code:
        exiting("Zipping error.")

    try:
        shutil.rmtree(run_output_path)
    except:
        exiting("Run_output folder couldn't be removed")

else:
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
