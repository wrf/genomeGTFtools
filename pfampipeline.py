#!/usr/bin/env python
#
# pfampipeline.py

'''pfampipeline.py  last modified 2016-04-22

pfampipeline.py proteins.fasta

    requires Bio Python library

    requires PFAM-A hmm database, found at:
  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

    optional SignalP (v4.1) stand-alone binary, academic download link at:
  http://www.cbs.dtu.dk/services/SignalP/

    hmmscan (from hmmer package at http://hmmer.org/) must be in PATH
'''

import sys
import os
import subprocess
import argparse
import time

def call_hmmscan(inputfasta, threadcount, pfamhmms):
	outfile = "{}.pfam.tab".format(os.path.splitext(inputfasta)[0])
	hmmscan_args = ["hmmscan","--cpu", threadcount, "--domtblout", outfile, pfamhmms, inputfasta]
	DEVNULL = open(os.devnull, 'w')
	print >> sys.stderr, "# Searching PFAM against {}".format(inputfasta), time.asctime()
	print >> sys.stderr, "Calling:\n{}".format(' '.join(hmmscan_args))
	subprocess.call(hmmscan_args, stdout=DEVNULL)
	return outfile

def call_pfamgff(hmmtable):
	outfile = "{}.gff".format(os.path.splitext(hmmtable)[0])
	pfam2gff_args = ["pfam2gff.py", "-i", hmmtable]
	print >> sys.stderr, "Calling:\n{}".format(' '.join(pfam2gff_args))
	with open(outfile, 'w') as pfamgff:
		subprocess.call(pfam2gff_args, stdout=pfamgff)
	return outfile

def call_pfamcdd(pfamgff, clanlinks, inputprots):
	outfile = pfamgff.replace("pfam","clan")
	pfam2cdd_args = ["pfamgff2clans.py", "-i", pfamgff, "-c", clanlinks, "-s", inputprots]
	print >> sys.stderr, "Calling:\n{}".format(' '.join(pfam2cdd_args))
	with open(outfile, 'w') as pfamcdd:
		subprocess.call(pfam2cdd_args, stdout=pfamcdd)
	return outfile

def call_signalp(signalp, inputfasta, cddgff, dscorecutoff):
	if not os.path.exists(signalp):
		print >> sys.stderr, "# Cannot find {}, skipping...".format(signalp), time.asctime()
		return 1
	signalp_args = [signalp, inputfasta]
	print >> sys.stderr, "Calling:\n{}".format(' '.join(signalp_args))
	spcall = subprocess.Popen(signalp_args, stdout=subprocess.PIPE)
	spstdoutlines = spcall.communicate()[0]
	# name                     Cmax  pos  Ymax  pos  Smax  pos  Smean   D     ?  Dmaxcut    Networks-used
	# scict1.028757.1_0          0.229  41  0.371  41  0.814  34  0.381   0.375 N  0.500      SignalP-TM
	# 0                     1        2     3        4     5        6     7        8        9    10       11
	# ['scict1.028757.1_0', '0.229', '41', '0.371', '41', '0.814', '34', '0.381', '0.375', 'N', '0.500', 'SignalP-TM']
	with open(cddgff, 'a') as cdo:
		for line in spstdoutlines.split("\n"):
			line = line.strip()
			if line and not line[0]=="#":
				lsplits = line.split()
				protid = lsplits[0]
				dscore = float(lsplits[8])
				if dscore >= dscorecutoff: # default for SignalP is 0.5
					cutpos = int(lsplits[4]) - 1
					print >> cdo, "{0}\tSignalP\tsignal_peptide\t1\t{1}\t{2}\t.\t.\tID={0}.sp".format(protid, cutpos, dscore)
	return 0

def call_draw_domains(drawdomains, cddgff):
	drawdomain_args = ["Rscript", drawdomains, cddgff]
	print >> sys.stderr, "Calling:\n{}".format(' '.join(drawdomain_args))
	subprocess.call(drawdomain_args)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input', help="fasta format file of proteins")
	parser.add_argument('-c','--clans', default=os.path.expanduser("~/db/Pfam-A.clans.tsv"), help="PFAM clan information tsv")
	parser.add_argument('-d','--d-score', type=float, default=0.25, help="D-score cutoff for SignalP")
	parser.add_argument('-n','--remove-null', action="store_true", help="remove domains that have no annotation")
	parser.add_argument('-p','--processors', help="number of CPUs for hmmscan [1]", default="1")
	parser.add_argument('-P','--PFAM', default=os.path.expanduser("~/PfamScan/data/Pfam-A.hmm"), help="PFAM-A hmm file")
	parser.add_argument('-R','--rscript', default=os.path.expanduser("~/Dropbox/rscripts/draw_protein_gtf.R"), help="path for Rscript of draw_protein_gtf.R")
	parser.add_argument('-S','--signalp', default=os.path.expanduser("~/signalp-4.1/signalp"), help="path for signalp")
	parser.add_argument('-t','--type', help="gff type [CDD]", default="CDD")
	args = parser.parse_args(argv)

	hmmtblout = call_hmmscan(args.input, args.processors, args.PFAM)
	pfamgff = call_pfamgff(hmmtblout)
	clangff = call_pfamcdd(pfamgff, args.clans, args.input)
	if args.signalp:
		signalp_flag = call_signalp(args.signalp, args.input, clangff, args.d_score)
	### TODO do something with signalp_flag
	call_draw_domains(args.rscript, clangff)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
