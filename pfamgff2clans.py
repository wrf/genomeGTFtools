#!/usr/bin/env python
#
# pfamgff2clans.py v1.0 created 2016-04-19

'''
pfamgff2clans.py  last modified 2016-04-20

pfamgff2clans.py -i proteins.pfam.gtf -c Pfam-A.clans.tsv > proteins.clan.gtf

    GENERATE PFAM GFF BY:
pfam2gff.py -i proteins.pfam.tab > proteins.pfam.gtf

    download PFAM clan information from:
ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
'''

import sys
import time
import argparse
import re
from collections import defaultdict
from Bio import SeqIO

def parse_clan_links(clanlinks):
	'''read in PFAM ID to PFAM clan tsv and make a dict where keys are PFAM accessions and values are cl accessions'''
    # The columns are: Pfam accession, clan accession, clan ID, Pfam ID, Pfam description.
	pfamtoclan = {}
	pfamannotation = {}
	print >> sys.stderr, "# Parsing clan links from {}".format(clanlinks), time.asctime()
	for line in open(clanlinks, 'r').readlines():
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			pfamacc = lsplits[0] # not always pfam
			clanacc = lsplits[1]
			clanname = lsplits[2]
			pfamname = lsplits[3]
			if clanacc: # some entries are blank, so skip them
				pfamtoclan[pfamacc] = clanacc
			if pfamname.startswith("DUF") and clanname:
				pfamannotation[pfamacc] = clanname
			else:
				pfamannotation[pfamacc] = pfamname
	print >> sys.stderr, "# Found {} clan links".format(len(pfamtoclan)), time.asctime()
	return pfamtoclan, pfamannotation

def parse_pfam_gtf(pfamgtf, overlaplimit, verbose=False):
	'''read PFAM GTF, merge identical annotations, and print the domain-merged GTF'''
	gtfbyprot = defaultdict(list) # keys are protein IDs and 
	domcount = 0
	print >> sys.stderr, "# Parsing GTF from {}".format(pfamgtf), time.asctime()
	for line in open(pfamgtf,'r').readlines():
		line = line.strip()
		if line and not line[0]=="#":
			domcount += 1
			lsplits = line.split("\t")
			protid = lsplits[0]
			domstart = int(lsplits[3])
			domend = int(lsplits[4])
			domlength = domend - domstart + 1
			qscore = float(lsplits[5])
			for i,protstats in enumerate(gtfbyprot[protid]):
				sstart, send, sscore = protstats[3:6] # meaning 3,4,5
				if sstart > domend or domstart > send: # means zero overlap
					continue
				else: # some overlap possible
					overlap = min(send, domend) - max(sstart, domstart) + 1
					if verbose:
						print >> sys.stderr, "{} {} overlap from ({},{}) to ({},{})".format(protid, overlap, domstart, domend, sstart, send)
					slength = send - sstart + 1
					qoverlap = overlap * 1.0 / domlength
					soverlap = overlap * 1.0 / slength
					if qoverlap >= overlaplimit: # default is 0.67 overlap, maybe needs to be lower
						if qscore < sscore: # worse domain hit, break
							break
					if soverlap >= overlaplimit:
						if qscore >= sscore:
							gtfbyprot[protid].pop(i)
			else: # if no break, add to list
				lsplits[3] = domstart # write back integer and floats
				lsplits[4] = domend
				lsplits[5] = qscore
				gtfbyprot[protid].append(lsplits)
	print >> sys.stderr, "# Found {} domains for {} proteins".format(domcount, len(gtfbyprot) ), time.asctime()
	return gtfbyprot

def convert_domains(domainsbyprot, programname, outputtype, wayout, pfamtoclandict, annotdict, fastalendict=None):
	print >> sys.stderr, "# Coverting domains to clans", time.asctime()
	writecount = 0
	for protid, domlist in domainsbyprot.iteritems():
		if fastalendict: # if length is available
			# print one entry for each protein
			# this could also be id: SO:0000104 polypeptide
			print >> wayout, "{0}\t{1}\tprotein\t1\t{2}\t.\t.\t.\tID={0}".format(protid, programname, fastalendict[protid])
		domaincounter = defaultdict(int)
		for domainstats in domlist:
			attributes = domainstats[8]
			hitid = re.search('ID=([\w.|-]+)', attributes).group(1)
			pfamid, targetname, domnumber = hitid.split('.') # ID should appear as ID={6}.{7}.{8};
			cldomain = pfamtoclandict.get(pfamid,pfamid) # if no clan is found, use the PFAM ID
			domaincounter[cldomain] += 1
			writecount += 1
			domainstats[8] = "ID={}.{}.{}".format(cldomain, annotdict.get(pfamid,"None"), domaincounter[cldomain])
			print >> wayout, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(*domainstats)
	print >> sys.stderr, "# Wrote {} domains".format(writecount), time.asctime()
	# NO RETURN

def get_prot_lengths(sequences):
	'''from a fasta file, return a dictionary where protein ID is the key and length is the value'''
	seqlendict = {}
	print >> sys.stderr, "# Parsing proteins from {}".format(sequences), time.asctime()
	for seqrec in SeqIO.parse(sequences,'fasta'):
		seqlendict[seqrec.id] = len(seqrec.seq)
	return seqlendict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="PFAM protein gff")
	parser.add_argument('-c','--clans', help="PFAM clan information tsv")
	parser.add_argument('-o','--overlap', type=float, default=0.67, help="minimum overlap to try and merge [0.67]")
	parser.add_argument('-p','--program', help="program for 2nd column in output [hmmscan]", default="hmmscan")
	parser.add_argument('-s','--sequences', help="fasta format file of protein sequences")
	parser.add_argument('-t','--type', help="gff type [PFAM]", default="PFAM")
	args = parser.parse_args(argv)

	pfamtoclandict, pfamannot = parse_clan_links(args.clans)

	seqlens = get_prot_lengths(args.sequences) if args.sequences else None

	domainsbyprot = parse_pfam_gtf(args.input, args.overlap)
	convert_domains(domainsbyprot, args.program, args.type, wayout, pfamtoclandict, pfamannot, seqlens)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
