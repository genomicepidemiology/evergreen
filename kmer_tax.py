#!/usr/bin/env python2.7

from __future__ import print_function
import sys, os, time
import argparse
import tempfile
import pickle
import shutil
import shlex
import subprocess
import sqlite3
from joblib import Parallel, delayed
from multiprocessing import cpu_count
import config

"""
Obs: Offline version, the db is populated here with runs!!
"""


NAP = 30
WAIT = 600
J_LIMIT = cpu_count()
base_path = os.path.dirname(sys.argv[0]).rsplit("/",1)[0]
REF_GEN_DIR = os.path.join(base_path, "complete_genomes") #the ref genomes without the plasmids !!!
MAIN_SQL_DB = os.path.join(base_path, "results_db/evergreen.db")
KMA_SHM = os.path.join(base_path, "scripts/kma_shm")
KMA = os.path.join(base_path, "scripts/kma")

parser = argparse.ArgumentParser(
    description='Wrapper for KMA for the Evergreen pipeline')
parser.add_argument(
    '-i',
    dest="collection_file",
    default=None,
    help='Input file of collected runs')
parser.add_argument(
    '-l',
    dest="runs_filename",
    default=None,
    help='File with sample names and paths to raw runs. Files for same isolate given comma separated.')
parser.add_argument(
    '-o',
    dest="odir",
    default="/data/evergreen/output",
    help='Output directory')
parser.add_argument(
    '-db',
    dest="database",
    default="/data/evergreen/hr_database/current/bacterial_compl_genomes_hq99_k13_ATG",
    help="Database, absolute path")
parser.add_argument(
    '-f_db',
    dest="folder_database",
    default="/data/evergreen/hr_database/current/bacteria.folder.pic",
    help="Acc number to folder pickle")
parser.add_argument(
    '-fa_db',
    dest="fsa_database",
    default="/data/evergreen/hr_database/current/bacteria.fsa_name.pic",
    help="Acc number to filename pickle")
parser.add_argument(
    '-wdir',
    dest="wdir",
    default=None, # /data/evergreen/output
    help='Working directory, not required')
parser.add_argument(
    '-d', '--debug',
    dest="debug",
    action="store_true",
    help='Debug')
parser.add_argument(
    '-q',
    dest="quiet",
    action="store_true",
    help='Quiet')
args = parser.parse_args()

def timing(message):
    if not args.quiet:
        t1 = time.time()
        print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)

def jobstart(command):
    #job = print(command)
    cmd = shlex.split(command)
    job = subprocess.call(cmd) # if shell=True, then it returns when the shell returns!
    return job


def jobstart_safe_n_silent(command):
    """Launch process with output redirected to /dev/null, kill process after WAIT time"""
    try:
        devnull = open(os.devnull, 'w')
    except IOError as e:
        exiting("Devnull won't open.")
    #print(command, file=sys.stderr)
    tx = time.time()
    proc = subprocess.Popen(shlex.split(command), stdout=devnull, stderr=devnull)
    time.sleep(NAP)
    returncode = proc.poll()
    while returncode is None:
        if time.time() - tx > WAIT:
            proc.kill()
            #print("signalled", proc.pid)
            returncode = 5
            continue
        #print(proc.returncode)
        time.sleep(NAP)
        returncode = proc.poll()
    return returncode

def jobstart_silent(command):
    """Launch process with output redirected to /dev/null"""
    try:
        devnull = open(os.devnull, 'w')
    except IOError as e:
        exiting("Devnull won't open.")
    #print(command, file=sys.stderr)
    job = subprocess.call(shlex.split(command), stdout=devnull, stderr=devnull)
    devnull.close()
    return job

def get_template(filename):
    '''
    Parse the outputs from taxonomy and find closest matching templates.

    Based on code for the Evergreen service by Johanne Ahrenfeldt

    '''
    try:
        templates = []
        with open("{0}/{1}".format(wdir, filename), "r") as f:
            # Skip header
            _ = f.readline()
            # Find all related templates (>80% kmer template coverage)
            for l in f:
                #Template	Num	Score	Expected	Template_length	query_coverage	Coverage	Depth	tot_query_coverage	tot_coverage	tot_depth	q_value	p_value
                #CP000857.1 Salmonella enterica subsp. enterica serovar Paratyphi C strain RKS4594, complete genome	4866	4267322	684	155800	   94.19	   99.91	   27.39	   94.19	   99.91	   27.39	4265268.49	1.0e-26
                d = [x.strip() for x in l.split('\t')]
                if len(d) > 9:
                    template_coverage = float(d[9])
                    depth = float(d[10])
                    if (template_coverage > 87.752 and depth >= 11.0):
                        # take all templates above XXkmerID%, the HR on the db should ensure limited similarity btw the templates
                        tmp = d[0].split()
                        species = " ".join(tmp[1:3])
                        # convert the templ_acc to the folder/template name we have in the dir system
                        template = folder_db[tmp[0]]
                        templates.append([filename[:-4], inputs[filename[:-4]], template, species, template_coverage, depth, tmp[0]])
        if not len(templates) > 0:
            print('Warning: No appropriate template could be found for {0}!'.format(filename), file=sys.stderr)
            return
        else:
            return templates
    except IOError as e:
        sys.stderr.write('All or some Kmer Results could not be found for %s!\n%s'%(filename, e))
        return

def prepare_tmpl_file(s, i, templ_acc):
    tfn = "{0}/kmer_{1}_{2}_{3}.tmpl".format(odir, time.strftime("%d%m%Y"), os.getppid(), i)
    tfiles.append(tfn)
    todirs.append(os.path.join(results_dir, s))
    if not os.path.isdir(todirs[-1]):
        try:
            os.mkdir(todirs[-1])

            #copy template sequence to the new directory
            if templ_acc.find(".") == -1:
                # oop it's the folder name already
                tmp = templ_acc.split("_")[-3:]
                templ_acc = "{}_{}.{}".format(tmp[0], tmp[1], tmp[2])

            src_fsa_file = os.path.join(REF_GEN_DIR, fsa_db[templ_acc])
            dest_file = os.path.join(todirs[-1], "pt_{0}.fa".format(s))
            shutil.copy(src_fsa_file, dest_file)

        except OSError as e:
            exiting(str(e))
    return tfn

##MAIN##
t0 = time.time()

# Check db
if not os.access("{0}.sparse.b".format(args.database), os.R_OK):
    exiting("Database unaccessible.")
if not os.access(args.folder_database, os.R_OK):
    exiting("Binary for folders unaccessible.")

tmpdir = ""
if os.environ.get('TMPDIR') is not None:
    tmpdir = os.path.join(os.environ['TMPDIR'], "kma-%s"%(os.getpid()))
    os.mkdir(tmpdir)
elif os.environ.get('PBS_JOBID') is not None:
    tmpdir = os.path.join(config.TMP_FOLDER, "{0}_{1}".format(os.environ['PBS_JOBID'], "kma-%s"%(os.getpid())))
    os.mkdir(tmpdir)
else:
    tmpdir = tempfile.mkdtemp(dir=config.TMP_FOLDER)

# If we want to keep some results, then it's in wdir
wdir = tmpdir
if args.wdir is not None and os.path.isdir(args.wdir):
    wdir = os.path.realpath(args.wdir)
    shutil.rmtree(tmpdir)
if os.path.isdir(args.odir):
    odir = os.path.realpath(args.odir)
results_dir = os.path.realpath(os.path.join(odir, "../results_db")) # results folder rel. path

# open sql db to templates
#print(MAIN_SQL_DB)

#start up the db during first run_id
conn = sqlite3.connect(MAIN_SQL_DB)
conn.execute("PRAGMA foreign_keys = 1")
cur = conn.cursor()

cur.execute('''CREATE TABLE IF NOT EXISTS runs
    (sra_id TEXT PRIMARY KEY,
    folder TEXT,
    file TEXT,
    cdate TEXT DEFAULT CURRENT_DATE,
    included INTEGER DEFAULT NULL);''')
conn.commit()

cur.execute('''CREATE TABLE IF NOT EXISTS templates
    (db_id INTEGER PRIMARY KEY,
    sra_id TEXT,
    template TEXT,
    species TEXT,
    templ_cov REAL,
    depth REAL,
    qc_pass INT DEFAULT 1,
    UNIQUE(sra_id, template),
    FOREIGN KEY(sra_id) REFERENCES runs(sra_id) ON DELETE CASCADE
)''')

# Process collection file
inputs = {} # acc : "/path/to/read_1.fastq.gz /path/to/read_2.fastq.gz"
runs_insert = []
if args.collection_file is not None:
    try:
        ifile = open(args.collection_file, "r")
    except IOError as e:
        exiting("Can't open file, {}".format(e))

    ddir = ""
    extensions = [".fastq.gz", "_1.fastq.gz", "_2.fastq.gz"]
    for line in ifile:
        if line[0] == "#":
            ddir = line[1:].strip()
        else:
            acc = line.strip()
            if acc not in inputs: #there migth be duplicates in the list
                tmp = []
                for ext in extensions:
                    fp = os.path.join(ddir, acc + ext)
                    if os.path.isfile(fp):
                        tmp.append(fp)
                inputs[acc] = " ".join(tmp)
                runs_insert.append((acc, ddir, ",".join(tmp)))
    ifile.close()

elif args.runs_filename is not None:
    with open(args.runs_filename, "r") as fp:
        for line in fp:
            # sample_id    pathtoread1,pathtoread2
            tmp = line.strip().split("\t")
            sid = tmp[0]
            ddir = os.path.dirname(tmp[1].split(",")[0])
            inputs[sid] = tmp[1].replace(","," ") #/path/to/read_1.fastq.gz /path/to/read_2.fastq.gz
            runs_insert.append((sid, ddir, tmp[1]))


else:
    exiting("Paths to input files are missing.")

cur.executemany('''INSERT OR REPLACE INTO runs (sra_id, folder, file) VALUES (?,?,?)''', runs_insert)
cur.execute('''CREATE INDEX IF NOT EXISTS sra_id_index ON runs(sra_id)''')
conn.commit()

timing('# Preparation done.')

###Set up shared db and then run KMA on it
# /home/projects/cge/apps/kma_shm
# -t_db /data/evergreen/hr_database/current/bacterial_compl_genomes_hr99_and_k13_ATG
# -Sparse


# /home/projects/cge/apps/kma
# -i /data/evergreen/data-3/ERR016880_1.fastq.gz /data/evergreen/data-3/ERR016880_2.fastq.gz
# -o /data/evergreen/test/ERR016880
# -t_db /data/evergreen/hr_database/current/bacterial_compl_genomes_hr99_and_k13_ATG
# -Sparse
# -shm
###

old_templs = []
kmf_cmds = []
spa_files = []
acc_found = set()
for acc in inputs.keys():
    # don't run kma again if there is a included result from it already
    # obv. it would not be added from q_p normally
    # check against db
    cur.execute("SELECT template FROM templates WHERE sra_id=? and qc_pass=1", (acc,))
    existing_templ = cur.fetchall()
    if existing_templ:
        acc_found.add(acc)
        for templ in existing_templ:
            old_templs.append([acc, inputs[acc], str(templ[0])])
            #print(old_templs[-1])
    else:
        # if no accepted KMA result yet in the db then do over
        tax_tsv_path = os.path.join(wdir, "{0}.spa".format(acc))

        if not os.path.exists(tax_tsv_path) or not os.path.getsize(tax_tsv_path) > 0:
            km_cmd = "{0} -i {1} -o {2}/{3} -t_db {4} -Sparse -shm".format(
            KMA, inputs[acc], wdir, acc, args.database
            )
            kmf_cmds.append(km_cmd)

        # km_cmd = "{0} -i {1} -o {2}/{3} -t_db {4} -Sparse -shm".format(
        # KMA, inputs[acc], wdir, acc, args.database
        # )
        # kmf_cmds.append(km_cmd)

        spa_files.append("{0}.spa".format(acc))


# Start kma if it needs to be run
if kmf_cmds:
    # load the db into shared memory
    shm_cmd = "{0} -t_db {1} -Sparse".format(KMA_SHM, args.database)
    task = jobstart_silent(shm_cmd)
    if task:
        exiting("KMA db couldn't be loaded into memory.")

    jobs = Parallel(n_jobs=J_LIMIT)(delayed(jobstart_safe_n_silent)(cmd) for cmd in kmf_cmds)
    if sum(jobs) != 0:
        print("Warning: A subprocess was unsuccessful.", file=sys.stderr)

    #delete shared mem db
    shm_cmd += " -destroy"
    task = jobstart_silent(shm_cmd)
    if task:
        exiting("KMA db couldn't be unloaded from memory.")

    timing("# KMA finished running.")

# load the folder databases
folder_db = {}
with open(args.folder_database, "r") as folder_db_fh:
    folder_db = pickle.load(folder_db_fh)

# Parse KMA results
runs_update = []
ut = []
ref_found = 0
if spa_files:
    references = Parallel(n_jobs=J_LIMIT)(delayed(get_template)(fn) for fn in spa_files)

    # Untangle templates
    templ_insert = []
    for isolate in references:
        if isolate is not None: # template(s) were found
            ref_found += 1

            if type(isolate[0]) is str: # just one
                ut.append(isolate)
                # sra_id, templ, species, tot_cov, depth
                templ_insert.append(tuple([isolate[0]]+isolate[2:6]))
                # sra_id, "1" if included
                runs_update.append((1, isolate[0]))
                acc_found.add(isolate[0])
            else:
                ut.extend(isolate)
                for tmpl in isolate:
                    templ_insert.append(tuple([tmpl[0]]+tmpl[2:6]))
                runs_update.append((1, isolate[0][0]))
                acc_found.add(isolate[0][0])



    timing("# KMA done, {0} runs have templates.".format(ref_found))

    # update sql db with the templates and whether the runs have been included
    if templ_insert:
        try:
            cur.executemany('''INSERT OR REPLACE INTO templates (sra_id, template, species, templ_cov, depth) VALUES (?,?,?,?,?)''', templ_insert)
            cur.execute('''CREATE INDEX IF NOT EXISTS db_id_index ON templates(db_id)''')
            conn.commit()
        except sqlite3.Error as e:
            exiting("SQL insert into templates failed: " + str(e))


acc_not_found = set(inputs.keys()) - acc_found
for acc in acc_not_found:
    runs_update.append((0, acc))

if runs_update:
    try:
        cur.executemany('''UPDATE runs SET included=? WHERE sra_id=?''', runs_update)
        conn.commit()
        conn.close()
    except sqlite3.Error:
        print("Warning: SQL update failed.", file=sys.stderr)
        # should not happen

timing("# Db updates/inserts are done.")

# add old templates to the list of templates
ut.extend(old_templs)
if not ut:
    exiting("No templates were found.")


# load reference genomes db
fsa_db = {}
with open(args.fsa_database, "r") as fsa_db_fh:
    fsa_db = pickle.load(fsa_db_fh)

# acc, inputs[acc], template, species, template_coverage, depth, template_acc
templ_sort = sorted(ut, key=lambda x: x[2])
s = templ_sort[0][2]
i = 0

#check if the results_db folder exists, if not, create it and copy the reference
# create input files for adt
tfiles = []
todirs = []
try:
    tfilename = prepare_tmpl_file(s, i, templ_sort[0][-1])
    tfile = open(tfilename, "w")
    for tmpl in templ_sort:
        if tmpl[2] != s:
            s = tmpl[2]
            tfile.close()
            i += 1
            tfilename = prepare_tmpl_file(s, i, tmpl[-1])
            tfile = open(tfilename, "w")
        # acc, inputs, template
        print("\t".join([str(x) for x in tmpl[:3]]), file=tfile)
    tfile.close()
except Exception as e:
    exiting("Error, {}".format(e))

if not args.debug:
    shutil.rmtree(tmpdir, ignore_errors=True)

timing("# Templates are found:\n{0}\n\n{1}".format("\n".join(tfiles), "\n".join(todirs)))


print("DONE", file=sys.stderr)
sys.exit(0)
