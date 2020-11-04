#!/usr/bin/env python
#
# pfampipeline.py

'''pfampipeline.py  last modified 2019-09-26

    USAGE requires only a fasta file of proteins
pfampipeline.py proteins.fasta

    PROGRAM WILL AUTOMATICALLY NAME ALL OUTPUT FILES
    if the starting file were named proteins.fasta, will produce:
    proteins.pfam.tab - hmmscan domain table
    proteins.pfam.gff - PFAM domains converted to protein GFF
    proteins.clan.gff - filtered PFAM domains, renamed by protein clan
    proteins.clan.pdf - visualised protein domain clans as PDF

    requires Bio Python library

    requires PFAM-A hmm database, found at:
  ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

    optional SignalP (v4.1) stand-alone binary, academic download link at:
  http://www.cbs.dtu.dk/services/SignalP/

    hmmscan (from hmmer package at http://hmmer.org/) must be in PATH

    this error usually arises from issues of the header with SignalP
    ValueError: could not convert string to float: N
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
	sys.stderr.write("# Searching PFAM against {}  ".format(inputfasta) + time.asctime() + os.linesep)
	sys.stderr.write("Calling:\n{}\n".format(' '.join(hmmscan_args)) )
	subprocess.call(hmmscan_args, stdout=DEVNULL)
	return outfile

def call_pfamgff(hmmtable, evalue):
	outfile = "{}.gff".format(os.path.splitext(hmmtable)[0])
	pfam2gff_args = ["pfam2gff.py", "-i", hmmtable, "-e", str(evalue)]
	sys.stderr.write("Calling:\n{}\n".format(' '.join(pfam2gff_args)) )
	with open(outfile, 'w') as pfamgff:
		subprocess.call(pfam2gff_args, stdout=pfamgff)
	return outfile

def call_pfamcdd(pfamgff, clanlinks, inputprots):
	outfile = pfamgff.replace("pfam","clan")
	pfam2cdd_args = ["pfamgff2clans.py", "-i", pfamgff, "-c", clanlinks, "-s", inputprots]
	sys.stderr.write("Calling:\n{}\n".format(' '.join(pfam2cdd_args)) )
	with open(outfile, 'w') as pfamcdd:
		subprocess.call(pfam2cdd_args, stdout=pfamcdd)
	return outfile

def call_signalp(signalp, inputfasta, gfftype, cddgff, dscorecutoff):
	if not os.path.exists(signalp):
		sys.stderr.write("# Cannot find {}, skipping...  ".format(signalp) + time.asctime() + os.linesep)
		return 1
	signalp_args = [signalp, inputfasta]
	sys.stderr.write("Calling:\n{}\n".format(' '.join(signalp_args)) )
	spcall = subprocess.Popen(signalp_args, stdout=subprocess.PIPE)
	spstdoutlines = spcall.communicate()[0].decode()
	# name                     Cmax  pos  Ymax  pos  Smax  pos  Smean   D     ?  Dmaxcut    Networks-used
	# scict1.028757.1_0          0.229  41  0.371  41  0.814  34  0.381   0.375 N  0.500      SignalP-TM
	# 0                     1        2     3        4     5        6     7        8        9    10       11
	# ['scict1.028757.1_0', '0.229', '41', '0.371', '41', '0.814', '34', '0.381', '0.375', 'N', '0.500', 'SignalP-TM']
	with open(cddgff, 'a') as cdo: # append additional signalP lines
		for line in spstdoutlines.split("\n"):
			line = line.strip()
			if line and not line[0]=="#":
				lsplits = line.split()
				protid = lsplits[0]
				dscore = float(lsplits[8]) # has problems with headers ###TODO convert to try
				if dscore >= dscorecutoff: # default for SignalP is 0.5
					cutpos = int(lsplits[4]) - 1
					cdo.write( "{0}\tSignalP\t{3}\t1\t{1}\t{2}\t.\t.\tID={0}.sp\n".format(protid, cutpos, dscore, gfftype) )
	return 0

def call_draw_domains(drawdomains, cddgff):
	drawdomain_args = ["Rscript", drawdomains, cddgff]
	sys.stderr.write("Calling:\n{}\n".format(' '.join(drawdomain_args)) )
	subprocess.call(drawdomain_args)

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input', help="fasta format file of proteins")
	parser.add_argument('-c','--clans', default=os.path.expanduser("~/db/Pfam-A.clans.tsv"), help="PFAM clan information tsv")
	parser.add_argument('-d','--d-score', type=float, default=0.25, help="D-score cutoff for SignalP")
	parser.add_argument('-e','--evalue', type=float, default=1e-1, help="evalue cutoff for domain filtering [1e-1]")
	parser.add_argument('-p','--processors', help="number of CPUs for hmmscan [1]", default="1")
	parser.add_argument('-P','--PFAM', default=os.path.expanduser("~/db/Pfam-A.hmm"), help="path to PFAM-A hmm file")
	parser.add_argument('-R','--rscript', default=os.path.expanduser("~/git/genomeGTFtools/draw_protein_gtf.R"), help="path for Rscript of draw_protein_gtf.R")
	parser.add_argument('-S','--signalp', default=os.path.expanduser("~/signalp-4.1/signalp"), help="path for signalp")
	parser.add_argument('-t','--type', help="gff type for SignalP [signal_peptide]", default="signal_peptide")
	args = parser.parse_args(argv)

	if not os.path.isfile(args.input):
		sys.exit("ERROR: CANNOT FIND INPUT FILE {}".format(args.input))
	if not os.path.isfile(args.clans):
		sys.exit("ERROR: CANNOT FIND CLAN LINKS FILE {}".format(args.clans))

	hmmtblout = call_hmmscan(args.input, args.processors, args.PFAM)
	pfamgff = call_pfamgff(hmmtblout, args.evalue)
	clangff = call_pfamcdd(pfamgff, args.clans, args.input)
	if args.signalp and os.path.exists(os.path.expanduser(args.signalp)):
		signalp_flag = call_signalp(args.signalp, args.input, args.type, clangff, args.d_score)
	### TODO do something with signalp_flag
	call_draw_domains(args.rscript, clangff)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
