#!/usr/bin/env python
#
# microsynteny_to_ortho_pairs.py v1 2023-07-31

"""microsynteny_to_ortho_pairs.py  last modified 2023-07-31

    generate a fasta file for each pair of putative orthologs based on microsynteny.py

microsynteny_to_ortho_pairs.py -m spA_against_spB_microsynteny.tab -f spA_prots.fa spB_prots.fa

"""

import sys
import os
import argparse
import gzip
from Bio import SeqIO

def safe_mkdir(path):
	if os.path.isfile(path):
		raise ValueError("{0} is a regular file, a directory is required".format(path))
	elif os.path.isdir(path):
		print("# The directory {0} already exists, and will be used".format(path), file=sys.stderr)
	else:
		print("# Making the directory {0}".format(path), file=sys.stderr)
		os.mkdir(path)

def get_protein_dict(seq_file_list):
	"""get list of fasta files and return a SeqIO dict containing all proteins"""
	dbdictionary = dict()
	for fastadb in seq_file_list:
		print("# Reading proteins from {}".format( fastadb ), file=sys.stderr)
		dbdictionary.update( SeqIO.to_dict(SeqIO.parse(fastadb,"fasta")) )
	print("# Counted {} proteins from all files".format( len(dbdictionary) ), file=sys.stderr)
	return dbdictionary

def parse_microsynteny( microsyn_file, prot_dict, cluster_dir, do_print_doubles ):
	"""from the microsynteny table, create fasta files of each syntenic pair"""

	# 0             1           2       3               4       5     6   7               8       9     10  11
	#NC_058106.1  NC_056533.1  blk-1  XP_044205213.1  215858  215890  -  XP_042271582.1  36087   36119   -  526.0
	#NC_058106.1  NC_056533.1  blk-1  XP_044207100.1  224447  224464  +  XP_042260048.1  44421   44423   +  1561.0
	#NC_058106.1  NC_056533.1  blk-1  XP_044212745.1  359975  360078  +  XP_042273658.1  168171  168274  +  940.0
	#NC_058106.1  NC_056533.1  blk-1  XP_044228147.1  482097  482173  -  XP_042271438.1  418155  418231  -  627.0

	gene_pair_count = 0
	block_counts = {} # key is block ID, value is count
	used_query_genes = {} # key is query ID, value is count
	used_subject_genes = {} # key is subject ID, value is count
	if microsyn_file.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		_opentype = gzip.open
		print("# Starting parsing on {} as gzipped".format( microsyn_file ), file=sys.stderr)
	else: # otherwise assume normal open for fasta format
		_opentype = open
		print("# Starting parsing on {}".format( microsyn_file ), file=sys.stderr)
	for line in _opentype(microsyn_file,'rt'):
		line = line.strip()
		if line: # skip blanks, there should not be any
			gene_pair_count += 1
			lsplits = line.split('\t') # should split to 12 columns
			q_prot_id = lsplits[3]
			s_prot_id = lsplits[7]
			used_query_genes[q_prot_id] = used_query_genes.get(q_prot_id,0) + 1
			used_subject_genes[s_prot_id] = used_subject_genes.get(s_prot_id,0) + 1

			q_scaffold = lsplits[0]
			s_scaffold = lsplits[1]
			block_num = lsplits[2]
			block_counts[block_num] = block_counts.get(block_num,0) + 1

			records = [q_prot_id, s_prot_id]
			# each file appears as:
			# NC_058106.1_NC_056533.1_blk1_g1_pair.fasta
			clustoutname = "{}_{}_{}_g{}_pair.fasta".format(q_scaffold, s_scaffold, block_num.replace('-',''), block_counts.get(block_num) )

			# write fasta format of each cluster
			clust_fasta_file = os.path.join(cluster_dir, clustoutname)
			with open(clust_fasta_file, 'w') as fo:
				for record in records:
					fo.write( prot_dict[record].format("fasta") )
	print("# Counted {} queries and {} subjects".format( len(used_query_genes) , len(used_subject_genes) ), file=sys.stderr)
	print("# Wrote {} protein pairs".format( gene_pair_count ), file=sys.stderr)
	if do_print_doubles:
		for k,v in used_query_genes.items():
			if v > 1:
				print("q\t{}\t{}".format(k,v) , file=sys.stderr )
		for k,v in used_subject_genes.items():
			if v > 1:
				print("s\t{}\t{}".format(k,v) , file=sys.stderr )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-m','--microsynteny', help="microsynteny table, can be .gz")
	parser.add_argument('-f','--fasta', nargs='*', help="original protein fasta sequences, concatenated or several files, IDs must be unique")
	parser.add_argument('-o','--output-dir', help="directory of output files [./blockpairs]", default='./blockpairs')
	parser.add_argument('-P','--print-doubles', help="print names of any sequence ID occurring more than once", action="store_true")
	args = parser.parse_args(argv)

	prot_dict = get_protein_dict(args.fasta)
	safe_mkdir(args.output_dir)
	parse_microsynteny( args.microsynteny, prot_dict, args.output_dir, args.print_doubles )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
