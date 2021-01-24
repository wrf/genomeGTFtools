#!/usr/bin/env python
#
# extract_coordinates.py

'''extract_coordinates.py  last modified 2021-01-24

    extract relevant regions from a GFF, to make simple figures
    resembling a view in a genome browser

    for example
    to extract coordinates from scaffold ML0011
    along the window of 63k-88k (includes 3 genes)

extract_coordinates.py -s ML0011 -b 63000 -e 88000 -g ML2.2.gff3 > ml0011_63k-88k_annot_short.tab

    -g can take multiple files to extract from different tracks at the same locus

    for bacterial genomes, where each gene is a single CDS, use -p

extract_coordinates.py -g GCF_000011805.1_ASM1180v1_genomic.gff -s NC_006841.2 -b 1041700 -e 1051700 -p > lux_locus_short_annot.tab
'''

import sys
import argparse
import gzip
from collections import defaultdict
#import matplotlib.pyplot as plt
#from matplotlib.patches import Rectangle
#from matplotlib.patches import Arrow
#from matplotlib.patches import Polygon

def extract_features(gtffile, target_scaffold, target_start, target_end, keep_only_full=False, gene_level_only=False, use_cds=False, is_bacterial=False):
	if gtffile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Parsing gff from {} as gzipped\n".format(gtffile) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Parsing gff from {}\n".format(gtffile) )

	linecounter = 0
	genecounter = 0
	scaffold_tracker = {} # key is scaffold, value is True
	# iterate through file
	for line in opentype(gtffile,'rt'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			feature = lsplits[2]
			scaffold = lsplits[0]
			scaffold_tracker[scaffold] = True
			if scaffold != target_scaffold:
				continue
			fstart = int(lsplits[3])
			fend = int(lsplits[4])
			if fend < target_start or fstart > target_end:
				# out of relevant bounds, so skip
				continue
			if keep_only_full:
				if fstart < target_start or fend > target_end:
					# either the start or end of the feature is out of bounds, so skip
					continue
			strand = lsplits[6]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
			if feature=="transcript" or feature=="mRNA":
				if feature=="transcript":
					feature = "mRNA"
				if gene_level_only:
					feature = "gene" # will be unique treated by R script
				gene_id = attrd.get("ID")
				outline = "{}\t{}\t{}\t{}\t{}\n".format( gene_id, feature, fstart, fend, strand )
				sys.stdout.write( outline )
				genecounter += 1
			elif feature=="exon" and not gene_level_only:
				gene_id = attrd.get("Parent")
				outline = "{}\t{}\t{}\t{}\t{}\n".format( gene_id, feature, fstart, fend, strand )
				sys.stdout.write( outline )
			elif feature=="CDS":
				if is_bacterial:
					gene_id = attrd.get("Parent")
					if gene_id is None: # some formats will use ID instead of Parent for CDS features
						gene_id = attrd.get("ID")
					outline = "{}\tgene\t{}\t{}\t{}\n".format( gene_id, fstart, fend, strand )
					sys.stdout.write( outline )
					genecounter += 1
				elif use_cds:
					gene_id = attrd.get("Parent")
					outline = "{}\texon\t{}\t{}\t{}\n".format( gene_id, fstart, fend, strand )
					sys.stdout.write( outline )
	sys.stderr.write("# Read {} lines, kept {} transcripts/genes\n".format(linecounter, genecounter) )
	if genecounter==0:
		if target_scaffold in scaffold_tracker:
			sys.stderr.write("# NO GENES KEPT, check -b or -e\n")
		else:
			sys.stderr.write("# NO GENES KEPT, SCAFFOLD:  {}  NOT FOUND, check -s\n".format(target_scaffold) )

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-g','--gff-files', nargs='*', help="GFF-type feature files, can be .gz")
	parser.add_argument('-b', '--begin', type=int, metavar='N', default=1, help="beginning scaffold position (starting with 1)")
	parser.add_argument('-e', '--end', type=int, metavar='N', default=1000000000, help="ending scaffold position")
	#parser.add_argument('-o','--output', help="name of output file (as .png)")
	parser.add_argument('-s','--scaffold', help="name of target scaffold")
	parser.add_argument('-f','--full-only', action="store_true", help="only keep features completely within the selected window")
	parser.add_argument('-c','--use_cds', action="store_true", help="use CDS features when drawing genes if exons are not given, not used with -G or -p")
	parser.add_argument('-G','--genes-only', action="store_true", help="draw only genes, ignore exon features")
	parser.add_argument('-p','--prokaryote', action="store_true", help="assume annotation is for a prokaryote, use gene or CDS features, and no exons or mRNA are specified")
	parser.add_argument('-v','--verbose', action="store_true", help="extra output")
	args = parser.parse_args(argv)

	sys.stderr.write("# Extracting features on {} from position {} to {}\n".format(args.scaffold, args.begin, args.end) )

	axis_line = "{}\taxis\t{}\t{}\t.\n".format( args.scaffold, args.begin, args.end )
	sys.stdout.write( axis_line )

	for gff_file in args.gff_files:
		extract_features( gff_file, args.scaffold, args.begin, args.end, args.full_only, args.genes_only, args.use_cds, args.prokaryote )



if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
