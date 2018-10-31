#!/usr/bin/env python

'''mcscan_to_gff.py  last modified 2018-10-31
convert mcscanX synteny output to GFF-like file

mcscan_to_gff.py -c g1_v_g2.collinearity > g1_v_g2.collinearity.gff
'''

import sys
import argparse

def read_scaffold_keys(renamefile):
	'''read short scaffold keys into dictionary, where key is short key, and value is full name'''
	scaffold_names = {} # key is short key (meaning can be indexed), value is long name
	# should be tab-delimited, in order of
	# short01   full_contig_name_01
	print >> sys.stderr, "# reading scaffold names from {}".format(renamefile)
	for line in open(renamefile,'r'):
		line = line.strip()
		if line:
			shortname, longname = line.split("\t")
			scaffold_names[shortname] = longname
	print >> sys.stderr, "# found {} short names".format(len(scaffold_names))
	return scaffold_names

def read_short_positions(shortposfile):
	'''read positions from psuedo-GFF, return dict where key is gene name and value is interval'''
	gene_intervals = {} # key is gene name, value is interval, scaffold is ignored
	# ta1	ta_g00001	11252	15952
	print >> sys.stderr, "# reading scaffold names from {}".format(shortposfile)
	for line in open(shortposfile,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			geneid = lsplits[1]
			boundaries = [ int(lsplits[2]), int(lsplits[3]) ]
			gene_intervals[geneid] = boundaries
	print >> sys.stderr, "# found {} gene positions".format(len(gene_intervals))
	return gene_intervals

def read_gene_positions(gfffile):
	'''read gene positions from GFF file, return dict where key is gene name and value is interval'''
	gene_intervals = {} # key is gene name, value is interval, scaffold is ignored
	## TODO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-c','--collinearity', help="collinearity output file of mcscanX")
	parser.add_argument('--transcript', default="transcript", help="optional name for transcript features, default is transcript")
	parser.add_argument('--short-names', help="optional short name conversion vector")
	parser.add_argument('--short-positions', help="gene positions, in short tab format")
	args = parser.parse_args(argv)

	scaffold_names = read_scaffold_keys(args.short_names) if args.short_names else None
	genepositions = read_short_positions(args.short_positions)

	linecounter = 0
	blockcounter = 0
	genecounter = 0
	lastalignment = ""
	genelist = [] # list of lists, each is a GFF line
	print >> sys.stderr, "# Converting {} to GFF".format(args.collinearity)
	for line in open(args.collinearity,'r'):
		line = line.strip()
		if line:
			linecounter += 1
			if linecounter < 12: # first 11 lines are summary stats
				continue
			if line[0:3]=="## ":
				# ## Alignment 0: score=11767.0 e_value=0 N=237 hh1&ta7 plus
				# if finding a new alignment, print the previous one
				if lastalignment:
					# print parent feature
					outline = "{0}\tmcscanX\tmatch\t{1}\t{2}\t{3}\t{4}\t.\tID=alignment_{5};Name=alignment_{5}_to_{6};Target={6} {7} {8};length={9};evalue={10}".format( queryscaffold, parentstart, parentend, fullscore, fullstrand, alignnum, subscaffold, subscafstart, subscafend, fulllength, fullevalue)
					print >> wayout, outline
					for gene in genelist:
						genecounter += 1
						outline = "\t".join(map(str,gene))
						print >> wayout, outline
				blockcounter += 1
				# reset parameters
				genelist = []
				parentstart, parentend = 0,0
				subscafstart, subscafend = 0,0
				# begin processing the line
				lsplits = line.split()
				alignnum = lsplits[2].replace(":","")
				lastalignment = alignnum
				#fullscore = re.search('score=([\d.]+)', line).group(1)
				fullscore = lsplits[3].split("=")[1]
				fullevalue = lsplits[4].split("=")[1]
				fulllength = lsplits[5].split("=")[1]
				queryscaffold, subscaffold = lsplits[6].split("&")
				# if short names are given, meaning dict is not None
				# then try to rename all scaffolds
				if scaffold_names:
					queryscaffold = scaffold_names.get(queryscaffold, queryscaffold)
					subscaffold = scaffold_names.get(subscaffold, subscaffold)
				fullstrand = "+" if lsplits[7]=="plus" else "-"
			else:
				#   0-  0:	hh_g00284	ta_g06734	      0
				lsplits = line.split("\t")
				genenum = lsplits[0].split("-")[1].strip().replace(":","")
				querygene = lsplits[1]
				subjectgene = lsplits[2]
				evalue = lsplits[3].strip()
				# adjust range of query genome
				querystart, queryend = genepositions.get(querygene)
				if parentstart==0 and parentend==0:
					parentstart = querystart
					parentend = queryend
				elif querystart > parentend:
					parentend = queryend
				elif queryend < parentstart:
					parentstart = querystart
				# adjust range of target genome
				subjectstart, subjectend = genepositions.get(subjectgene)
				if subscafstart==0 and subscafend==0:
					subscafstart = subjectstart
					subscafend = subjectend
				elif subjectstart > subscafend:
					subscafend = subjectend
				elif subjectend < parentstart:
					subscafstart = subjectstart
				# make list, to convert into GFF string
				features = [queryscaffold, "mcscanX", "match_part", querystart, queryend, evalue, fullstrand, ".", "ID=alignment_{0}.{1};Parent=alignment_{0};Target={2} {3} {4}".format(alignnum, genenum, subjectgene, subjectstart, subjectend) ]
				genelist.append(features)
	else: # print final gene
		outline = "{0}\tmcscanX\tmatch\t{1}\t{2}\t{3}\t{4}\t.\tID=alignment_{5};Name=alignment_{5}_to_{6};Target={6} {7} {8};length={9};evalue={10}".format( queryscaffold, parentstart, parentend, fullscore, fullstrand, alignnum, subscaffold, subscafstart, subscafend, fulllength, fullevalue)
		print >> wayout, outline
		for gene in genelist:
			genecounter += 1
			outline = "\t".join(map(str,gene))
			print >> wayout, outline
	print >> sys.stderr, "# read {} lines for {} genes in {} blocks".format(linecounter, genecounter, blockcounter)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
