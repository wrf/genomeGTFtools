#!/usr/bin/env python
#
# alignmentpos2gff.py v1.0 created 2017-09-19
#
# Sequence Ontology terms from:
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo

'''
alignmentpos2gff.py  last modified 2017-09-19

    EXAMPLE USAGE:
alignmentpos2gff.py -a all_prots.aln -s 100,150 > all_prots.sites.gff

    to exclude all amino acids except a certain set, say FYW
alignmentpos2gff.py -a all_prots.aln -s 100,150 -f FYW > all_prots.FYW_only.gff

    or append to the protein GFF files from pfam2gff.py or pfampipeline.py

alignmentpos2gff.py -a all_prots.aln -s 100,150 >> all_prots.clan.gff
'''

import sys
import time
import argparse
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-a','--alignment', help="alignment in fasta format")
	parser.add_argument('-f','--filter', default="ACDEFGHIKLMNPQRSTVWY", help="keep only certain amino acids [default: all]")
	parser.add_argument('-s','--sites', help="comma separate list of sites, like '75,79,150'")
	parser.add_argument('-i','--keep-index', action="store_true", help="sites given are already indexed for python (meaning minus 1)")
	args = parser.parse_args(argv)

	sites = [int(s) for s in args.sites.split(",")]
	if not args.keep_index:
		sites = [i-1 for i in sites]

	for seqrec in SeqIO.parse(args.alignment, "fasta"):
		#seqsites = [seqrec.seq[s] for s in sites]
		#seqpos = [len(str(x1.seq[:s+1]).replace("-","")) for s in sites]
		for ss in sites:
			residueatpos = seqrec.seq[ss]
			nogappos = len(str(seqrec.seq[:ss+1]).replace("-",""))
			if residueatpos in args.filter: # by default this will ignore any gaps
				# Q16665	UniProtKB	Modified residue	564	564	.	.	.	Note=4-hydroxyproline
				outline = "{0}\talignment\tmodified residue\t{1}\t{1}\t1\t.\t.\tNote={2}".format(seqrec.id, nogappos, residueatpos)
				print >> wayout, outline

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
