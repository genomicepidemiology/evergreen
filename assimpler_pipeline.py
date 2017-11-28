#!/usr/bin/env python2.7

# Import libraries
import sys, time
import os
import gc
from optparse import OptionParser
from operator import itemgetter
import re
from math import sqrt
import gzip

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
   >>> # MULTIFILE LOOP
   >>> for fil in ['test.fsa','test.fsa']:
   ...    for seq, name, desc in SeqsFromFile(fil):
   ...       print ">%s %s\n%s"%(name, desc, seq)
   ...

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
         #sys.stderr.write("%s\n"%line)
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

#
# Functions
#
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

# Construct the reverse complemant from a sequence
def reversecomplement(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if   s == 'A': comp = comp + 'T'
        elif s == 'T': comp = comp + 'A'
        elif s == 'C': comp = comp + 'G'
        elif s == 'G': comp = comp + 'C'
        else:          comp = comp + s
    return comp[::-1]

# Color to base
color2base = {
  "A0" : 'A',
  "A1" : 'C',
  "A2" : 'G',
  "A3" : 'T',
  "T0" : 'T',
  "T1" : 'G',
  "T2" : 'C',
  "T3" : 'A',
  "C0" : 'C',
  "C1" : 'A',
  "C2" : 'T',
  "C3" : 'G',
  "G0" : 'G',
  "G1" : 'T',
  "G2" : 'A',
  "G3" : 'C'
}

tooneletter = {
  #
  # A  adenosine          C  cytidine             G  guanine
  # T  thymidine          N  A/G/C/T (any)        U  uridine
  # K  G/T (keto)         S  G/C (strong)         Y  T/C (pyrimidine)
  # M  A/C (amino)        W  A/T (weak)           R  G/A (purine)
  # B  G/T/C              D  G/A/T                H  A/C/T
  # V  G/C/A              -  gap of indeterminate length
  #   A C G T
  # A A M R W
  # C M C S Y
  # G R S G K
  # T W Y K T
  #
  "AA" : 'A',
  "AC" : 'M',
  "AG" : 'R',
  "AT" : 'W',
  "AN" : 'A',
  "TA" : 'W',
  "TC" : 'Y',
  "TG" : 'K',
  "TT" : 'T',
  "TN" : 'T',
  "CA" : 'M',
  "CC" : 'C',
  "CG" : 'S',
  "CT" : 'Y',
  "CN" : 'C',
  "GA" : 'R',
  "GC" : 'S',
  "GG" : 'G',
  "GT" : 'K',
  "GN" : 'G',
  "NA" : 'A',
  "NC" : 'C',
  "NG" : 'G',
  "NT" : 'T',
  "NN" : 'N'
}

# Start time to keep track of progress
t0 = time.time()
etta = 0.001

# Parse command line options
parser = OptionParser()
parser.add_option("-i", "--inputfiles", dest="inputfiles", help="read from INFILE", metavar="INFILE")
parser.add_option("-t", "--templatefile", dest="templatefilename", help="read from TEMFILE", metavar="TEMFILE")
parser.add_option("-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE")
parser.add_option("-s", "--columnfile", dest="columnfilename", help="write to column file", metavar="COLUMN")
parser.add_option("-c", "--consensusfile", dest="consensusfilename", help="write fasta assembly to CONFILE", metavar="CONFILE")
parser.add_option("-n", "--entryname", dest="entryname", help="Set name in fastafile to ENTRYNAME", metavar="ENTRYNAME")
parser.add_option("-u", "--unusedfile", dest="unusedfilename", help="write to UNUSEDFILE", metavar="UNUSEDFILE")
parser.add_option("-m", "--mincoverage", dest="mincoverage", help="lower bound of coverage for assembly to print out. Dekault is 0 (~1 read )", metavar="MINCOVERAGE")
parser.add_option("-a", "--maxlatcalls", dest="maxaltcalls", help="maximum fraction of alternative calls. Default is 0.1", metavar="MAXALTCALLS")
parser.add_option("-z", "--minzscore", dest="minzscore", help="Minimum Z-score to call bases. Deafault is 3.29 ~ p = 0.001", metavar="MINZSCORE")
parser.add_option("-M", "--mismatchscore", dest="mismatchscore", help="Mismatch score (penalty). Deafault is 3", metavar="KMERLENGTH")
parser.add_option("-k", "--kmerlength", dest="kmerlength", help="k-mer length to use in hash table. default is 17", metavar="MISMATCHSCORE")
parser.add_option("-l", "--minscore", dest="minscore", help="Minscore to accept mapping. Deafault is 50", metavar="MINSCORE")
parser.add_option("-f", "--frogdna", dest="usefrogdna", action="store_true", help="Use template DNA in regions with no reads")
parser.add_option("-e", "--extendtemplate", dest="extendtemplate", action="store_true", help="extend templates if reads overhang in the end/create new templat if no match is found")
parser.add_option("-d", "--diploid", dest="diploid", action="store_true", help="Make base calls for diploid genomes")
(options, args) = parser.parse_args()
#
# Open files
#

# File with templare
if options.templatefilename != None:
  templatefile = open(options.templatefilename,"r")

# File for general output
if options.outputfilename != None:
  outputfile = open(options.outputfilename,"w")
else:
  outputfile = sys.stdout

# File for consensus sequence in fasta format
if options.consensusfilename != None:
  consensusfile = open(options.consensusfilename,"w")
else:
  consensusfile = sys.stdout

# Entryname in fastafile
if options.entryname != None:
  entryname = options.entryname + " "
else:
  entryname = ""

# File for assembly in columns in SAM like format (which is not really sam like)
if options.columnfilename != None:
  columnfile = open(options.columnfilename,"w")

# File for unused reads
if options.mincoverage != None:
  mincoverage = options.mincoverage
else:
  mincoverage = 0

if options.maxaltcalls != None:
  maxaltcalls = float(options.maxaltcalls)
else:
  maxaltcalls = 0.1
mincovfrac = float(1.0-maxaltcalls)

if options.minzscore != None:
  minzscore = float(options.minzscore)
else:
  minzscore = 3.29

# Penalty for mismatch
if options.mismatchscore != None:
  mismatchscore = -1.0*float(options.mismatchscore)
else:
  mismatchscore = -3

# Minscore to accept mapping
if options.minscore != None:
  minscore = 1.0*float(options.minscore)
else:
  minscore = 50

# K-mer length to use for mapping
if options.kmerlength != None:
  kmerlength = int(options.kmerlength)
else:
  kmerlength = 17

# Minimum coverage to print out
if options.unusedfilename != None:
  unusedfile = open(options.unusedfilename,"w")

# Read Template fasta file
templateseq = []
templateseqsegments = []
consensusseq = []
templatename = []
templatedesc = []
Ntemplates=0
i=0
if options.templatefilename != None:
   sys.stdout.write("%s\n" % ("# Reading templatefile"))
   for line in templatefile:
      fields=line.split()
      if len(line)>1:
         if fields[0][0] == ">":
            if (i>0):
               templateseq[-1] = ''.join(templateseqsegments)
            del templateseqsegments
            templateseqsegments = []
            i=0
            templateseq.append("")
            consensusseq.append("")
            templatename.append(fields[0][1:])
            templatedesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
         else:
            templateseqsegments.append("")
            templateseqsegments[i] = fields[0]
            i+=1
   templateseq[-1] = ''.join(templateseqsegments)
del templateseqsegments

# Make nmer index for templatefiles
sys.stdout.write("%s\n" % ("# Making index"))
templateindex = {}
# Length of K-mer
oligolen = kmerlength
# Steps between k-mers in template (large value reduces memory usage)
kmerstep = kmerlength
for i in xrange(0, len(templateseq)):
  seqlen = len(templateseq[i])
  for j in xrange(0, seqlen-oligolen+1,kmerstep):
    #take larger steps to conserve memory
    submer = templateseq[i][j:j+oligolen]
    name = ("%d%s%d" % (i,"_",j))
    if submer in templateindex:
      # below is just to test
      templateindex[submer].append(name)
    else:
       templateindex[submer] = [name,]

# Initialize read depth matrices
sys.stdout.write("%s\n" % ("# Initializing read statistics"))
Nmatch = []
NmatchA = []
NmatchT = []
NmatchG = []
NmatchC = []

# non-numpy version
for i in xrange(0, len(templateseq)):
#  Nmatch.append([0,0,0,0,0]*len(templateseq[i]))
  Nmatch.append([0]*len(templateseq[i]))
  NmatchA.append([0]*len(templateseq[i]))
  NmatchT.append([0]*len(templateseq[i]))
  NmatchG.append([0]*len(templateseq[i]))
  NmatchC.append([0]*len(templateseq[i]))
#
# Scan input file for matches to template
#
# Number of reads read
nreads =0
# Number of lookups of oligos in template databases
oligotests=0
# Number of oligos matching template
oligohits=0
# Number og hits that are significant after extension
fullhits=0
# Score for matching nucleotides
matchscore = 1
# Score for alitning to nucleotides with low read depth
singlescore = 0.1
# Number of unmatched reads
Nunmatchedreads= 0
# Min mength of read for it to be mapped
minreadlength = int(minscore)
# Stepsize for selecting next oligo for searching for match in template
oligosteps = 1
sys.stdout.write("%s\n" % ("# Scanning input file for matches to template"))
t1 = time.time()
for fil in options.inputfiles.split():
   for seq, name, desc in SeqsFromFile(fil):
      nreads += 1
      if nreads % 1000 == 0:
         t2 = time.time()
         sys.stdout.write("\r# %s nreads (%s nreads / s)" % ("{:,}".format(nreads), "{:,}".format(int(nreads / (t2-t1)))))
         sys.stdout.flush()
      # Search for sequence in template
      foundmatch = 0
      if len(seq) >= minreadlength:
         # Filter out short reads
         for qseq in [seq,reversecomplement(seq)]:
            # Try both the sequence and its reverse complement
            qseqlen = len(qseq)
            for i in xrange(0, qseqlen-oligolen,oligosteps):
               # Look for first exact match of a subsequence of length oligolen
               # Taking longer steps here may speed things up (and degrade perormance a bit!)
               if (foundmatch > 0):
                  break
               # Go to reverse complement or next sequence if mactch is found
               submerstart = i
               submer = qseq[submerstart:submerstart+oligolen]
               oligotests += 1
               if submer in templateindex:
                  #  Subsequence is found in template
                  foundmatch += 1
                  oligohits += 1
                  # Make vector of matches
                  matchindexes = templateindex[submer]
                  for matchindex in matchindexes:
                     # split each match into sequence No (j) and nucleotide number in sequence (k)
                     (jj,kk) = matchindex.split("_")
                     j = int(jj)
                     k = int(kk)
                     # Extend match to the left and right to find max scoring segment
                     #extend left
                     l = submerstart-1
                     m = k-1
                     leftscore = 0
                     bestleftscore = 0
                     bestleftqstart = submerstart
                     bestlefttstart = k
                     while (l >= 0 and m >= 0 and m > k-qseqlen):
                        # Search within sequeces and not longer than length of qseqlen away
                        # The last is really only needed if code is generalized to gaped alignments
                        # The following 10 lines is maybe where speed may be gained
                        if (qseq[l:l+1] == templateseq[j][m:m+1]):
                           leftscore += matchscore
                        else:
                           leftscore += mismatchscore
                        if (leftscore > bestleftscore):
                           bestleftscore = leftscore
                           bestleftqstart = l
                           bestlefttstart = m
                        l -= 1
                        m -= 1
                     #extend right - see comments above for extend left
                     l = submerstart+oligolen
                     m = k+oligolen
                     rightscore = 0
                     bestrightscore = rightscore
                     bestrightqend = l
                     bestrighttend = m
                     templateseqlen = len(templateseq[j])
                     while (l< qseqlen  and m < templateseqlen and m <= k+oligolen+qseqlen):
                        disttoend = (templateseqlen -1) - bestrighttend
                        if (qseq[l:l+1] == templateseq[j][m:m+1]):
                           rightscore += matchscore
                        # Allow well scoring alignments to be extended if it can extend long to the right of template
                        elif (disttoend < 5 and totalscore >= minscore and qseqlen-bestrightqend-1 > 5):
                           rightscore += singlescore
                        else:
                           rightscore += mismatchscore
                        if (rightscore > bestrightscore):
                           bestrightscore = rightscore
                           bestrightqend = l
                           bestrighttend = m
                        l += 1
                        m += 1
                     totalscore = bestleftscore+bestrightscore+oligolen*matchscore
                     if (totalscore >= minscore):
                        fullhits += 1
                        for n in xrange (0,bestrighttend-bestlefttstart):
                           Nmatch[j][n+bestlefttstart] += 1
                           nucleotide =  qseq[n+bestleftqstart:n+bestleftqstart+1]
                           if (nucleotide == 'A'):
                              NmatchA[j][n+bestlefttstart] += 1
                           elif (nucleotide == 'T'):
                              NmatchT[j][n+bestlefttstart] += 1
                           elif (nucleotide == 'G'):
                              NmatchG[j][n+bestlefttstart] += 1
                           elif (nucleotide == 'C'):
                              NmatchC[j][n+bestlefttstart] += 1
                        if options.extendtemplate != None:
                           # Extend reads to the right
                           disttoend = (templateseqlen -1) - bestrighttend
                           if disttoend == 0: # or <2 eg
                              # Alignment goes all the way to the end of the template
                              for o in xrange (disttoend,qseqlen-bestrightqend-1):
                                 tpos = templateseqlen+o
                                 qpos = bestrightqend+1+o
                                 templateseq[j] = templateseq[j] + qseq[qpos:qpos+1]
                                 Nmatch[j].append(0)
                                 NmatchA[j].append(0)
                                 NmatchT[j].append(0)
                                 NmatchG[j].append(0)
                                 NmatchC[j].append(0)
                                 # Update read depth statistics - i am not sure the below section works
                                 nucleotide = qseq[qpos:qpos+1]
                                 if (nucleotide == 'A'):
                                    NmatchA[j][tpos] += 1
                                 elif (nucleotide == 'T'):
                                    NmatchT[j][tpos] += 1
                                 elif (nucleotide == 'G'):
                                    NmatchG[j][tpos] += 1
                                 elif (nucleotide == 'C'):
                                    NmatchC[j][tpos] += 1
                                 # add kmer to index sequence (I should acrually go in kmer steps here)
                                 indexseq = qseq[qpos-oligolen+1:qpos+1]
                                 name = ("%d%s%d" % (j,"_",tpos-oligolen+1))
                                 if indexseq in templateindex:
                                    templateindex[indexseq].append(name)
                                 else:
                                    templateindex[indexseq] = [name,]
      if (foundmatch == 0):
         Nunmatchedreads += 1;
         if options.unusedfilename != None:
            unusedfile.write("%s\n" % (seq))
         # Add new template
         if options.extendtemplate != None:
            i = len(templateseq)
            templateseq.append("")
            consensusseq.append("")
            templatename.append("New template")
            templatedesc.append("")
            templateseq[i] = seq
            seqlen = len(seq)
            for j in xrange(0, seqlen-oligolen+1,kmerstep):
               # Take larger steps to conserve memory
               submer = templateseq[i][j:j+oligolen]
               name = ("%d%s%d" % (i,"_",j))
               if submer in templateindex:
                  templateindex[submer].append(name)
               else:
                  templateindex[submer] = [name,]
            Nmatch.append([0]*len(templateseq[i]))
            NmatchA.append([0]*len(templateseq[i]))
            NmatchT.append([0]*len(templateseq[i]))
            NmatchG.append([0]*len(templateseq[i]))
            NmatchC.append([0]*len(templateseq[i]))

# Print template based assembly
sys.stdout.write("\n")
sys.stdout.write("%s\n" % ("# Writing template based assembly"))
# total number of matches summed over all positions
Nmatchsum = 0
Nmatchsumv = []
# Number af positions covered by atleast one read
covered = 0
# Number af positions covered by atleast mincoverage reads
mincovered = 0
# Length of template
templatelentot = 0

mincoveredcon = []
for i in xrange(0, len(templateseq)):
  mincoveredcon.append("0")

# If user have selected sam style outputfile then write header
if options.columnfilename != None:
  columnfile.write("%s\n" % ("# index templatename templateseq consensusseq No_of_matches  No_of_A No_of_T No_of_G No_of_C"))

# Call bases
for i in xrange(0, len(templateseq)):
   # Loop over template sequences
   Nmatchsumv.append(0)
   consensuslist = []
   templatelentot += len(templateseq[i])
   mincoveredcon[i] =0
   for j in xrange (0,len(templateseq[i])):
      # Loop over positions in templatefile
      if options.usefrogdna != None:
         # Use template (frog in Jurassic Parc) DNA
         consensus = templateseq[i][j]
      else:
         # The more conservative choice: if not enough reads call it as an N
         consensus = 'N'
         nextbase = 'N'
      maxcoverage = mincoverage
      if options.diploid == None:
         # Not diploid genome
         if (float(NmatchA[i][j])>maxcoverage):
            maxcoverage = NmatchA[i][j]
            nextbase = 'A'
         if (float(NmatchT[i][j])>maxcoverage):
            maxcoverage = NmatchT[i][j]
            nextbase = 'T'
         if (float(NmatchG[i][j])>maxcoverage):
            maxcoverage = NmatchG[i][j]
            nextbase = 'G'
         if (float(NmatchC[i][j])>maxcoverage):
            maxcoverage = NmatchC[i][j]
            nextbase = 'C'
         tot = float(Nmatch[i][j])
         alt =  tot - maxcoverage
         zscore = (maxcoverage-alt)/sqrt(tot+etta)
         if (maxcoverage > mincoverage and maxcoverage/(tot+etta)>mincovfrac and zscore > minzscore):
            consensus = nextbase
      else:
         # Diploid genome base calls
         firstbase = 'N'
         secondbase = 'N'
         nfirstbase = 0
         nsecondbase = 0
         # Find most common base
         nfirstbase = mincoverage
         if (float(NmatchA[i][j])>nfirstbase):
            nfirstbase = NmatchA[i][j]
            firstbase = 'A'
         if (float(NmatchT[i][j])>nfirstbase):
            nfirstbase = NmatchT[i][j]
            firstbase = 'T'
         if (float(NmatchG[i][j])>nfirstbase):
            nfirstbase = NmatchG[i][j]
            firstbase = 'G'
         if (float(NmatchC[i][j])>nfirstbase):
            nfirstbase = NmatchC[i][j]
            firstbase = 'C'
         # Find second most common base
         nsecondbase = mincoverage
         if (float(NmatchA[i][j])>nsecondbase and firstbase != 'A'):
            nsecondbase = NmatchA[i][j]
            secondbase = 'A'
         if (float(NmatchT[i][j])>nsecondbase and firstbase != 'T'):
            nsecondbase = NmatchT[i][j]
            secondbase = 'T'
         if (float(NmatchG[i][j])>nsecondbase and firstbase != 'G'):
            nsecondbase = NmatchG[i][j]
            secondbase = 'G'
         if (float(NmatchC[i][j])>nsecondbase and firstbase != 'C'):
            nsecondbase = NmatchC[i][j]
            secondbase = 'C'
         # Find significance of base calls
         tot = float(Nmatch[i][j])
         alt =  tot - (nfirstbase+nsecondbase)
         zscore = (nfirstbase-alt)/sqrt(tot-nsecondbase+etta)
         if (nfirstbase > mincoverage and nfirstbase/(tot-nsecondbase+etta)>mincovfrac and zscore > minzscore):
            # Call for most comon base is significant
            zscore = (nsecondbase-alt)/sqrt(tot-nfirstbase+etta)
            if (nsecondbase > mincoverage and nsecondbase/(tot-nfirstbase+etta)>mincovfrac and zscore > minzscore):
               # heterozygote call Call for secondmost comon base is significant
               diconsensus = firstbase + secondbase
            else:
               # homozygore call
               diconsensus = firstbase + firstbase
         else:
            diconsensus = 'NN'
         consensus = tooneletter[diconsensus]
      consensuslist.append(consensus)
      # Write columns in column format
      if options.columnfilename != None:
         columnfile.write("%10d %20s %5d %1s %1s %5d %5d %5d %5d %5d\n" % (
            i+1,
            templatename[i],
            j+1,
            templateseq[i][j],
            consensus, Nmatch[i][j],
            NmatchA[i][j],
            NmatchT[i][j],
            NmatchG[i][j],
            NmatchC[i][j]))
      # Update asembly statistics
      if (Nmatch[i][j]>0):
         Nmatchsumv[i] +=  Nmatch[i][j]
         covered += 1
      if (Nmatch[i][j]>=mincoverage):
         mincovered += 1
         mincoveredcon[i] += 1
   # The following is the fast way to construct sequences according to an internet site, and it really speed things up
   consensusseq[i] = ''.join(consensuslist)
   del consensuslist
   Nmatchsum += Nmatchsumv[i]
   sys.stdout.write("# Average cov template no. %d: %f %s %s\n"%(
      i+1,
      float(Nmatchsumv[i])/(float(len(templateseq[i]))+etta),
      templatename[i],
      templatedesc[i]
      ))
   # Write consensusseq as fastafile
   if options.consensusfilename != None:
      if mincoveredcon[i] > 0:
         start=0
         consensusfile.write(">%sContig_%d %s %s\n" % (entryname,i+1, templatename[i], templatedesc[i]))
         start=0
         while start < len(templateseq[i]):
            consensusfile.write("%s\n" % (consensusseq[i][start:start+60]))
            start +=60

# Print statistics
sys.stdout.write("%s %d\n" % ("# Oligotests:   ", oligotests))
sys.stdout.write("%s %d\n" % ("# Fullhits:     ", fullhits))
sys.stdout.write("%s %d\n" % ("# Nmatchsum:    ", Nmatchsum))
sys.stdout.write("%s %f\n" % ("# Average cov: ", float(Nmatchsum)/float(templatelentot+etta)))
sys.stdout.write("%s %d\n" % ("# templatelentot:  ", templatelentot))
sys.stdout.write("%s %d\n" % ("# covered:      ", covered))
sys.stdout.write("%s %f\n" % ("# frac covered: ", float(covered)/float(templatelentot+etta)))
sys.stdout.write("%s %d\n" % ("# mincovered:   ", mincovered))
sys.stdout.write("%s %f\n" % ("# frac min covered: ", float(mincovered)/float(templatelentot+etta)))
sys.stdout.write("%s %d\n" % ("# N unmap read: ", Nunmatchedreads))

t2 = time.time()
sys.stdout.write("%s %d %s\n" % ("# Finishing. Time used: ", int(t2-t0)," seconds"))
sys.stderr.write("Done\n")
