#!/usr/bin/env python
# v1.0 created 2016-03-16

'''number_contigs_by_length.py    last modified 2019-10-04

number_contigs_by_length.py contigs.fasta > renumbered_contigs.fasta

    by default, prints fasta output as:
scaffold_0123

    broken down as:
    -n   -d        -z                 -l              -o
    name delimiter [omit-zero] number optional-length optional-old-name

    with all options, might appear as:
number_contigs_by_length.py -n avictoria -d "." -l length -o contigs.fasta > renumbered_contigs.fasta

    output would look like:
avictoria.0001.length.248369 scaffold4_cov199

    to output optional conversion vector list, use -c conversions.txt
    which is a tab separated list like:
oldnamecontig123    newnamecontig001

    conversion vector can be used with script - rename_gtf_contigs.py -
'''

import sys
import argparse
import gzip
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', help="contigs file, by default in fasta format, can be .gz")
	parser.add_argument('-c','--conversion', help="optional output file for naming conversions")
	parser.add_argument('-d','--delimiter', default='_', help="delimiter for headers")
	parser.add_argument('-f','--format', default='fasta', help="format of sequences [fasta] or fastq")
	parser.add_argument('-l','--length', nargs='?', const="len", help="include length in fasta ID, printed as [len]")
	parser.add_argument('-n','--name', help="name for output [scaffold]", default="scaffold")
	parser.add_argument('-o','--old', action="store_true", help="keep old name as space separated description")
	parser.add_argument('-R','--only-reorder', action="store_true", help="do not rename any sequences, only reorder by length")
	parser.add_argument('-z','--omit-zero', action="store_true", help="do not add extra zeroes to number")
	args = parser.parse_args(argv)

	# make dict of all contigs
	if args.input_file.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading sequences from {} as gzipped\n".format(args.input_file) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading sequences from {}\n".format(args.input_file) )
	contigdict = SeqIO.to_dict( SeqIO.parse( opentype(args.input_file,'rt'), args.format) )
	# count number of sequences, and return number of padded zeroes
	if not args.omit_zero:
		reclog = str(len(str(len(contigdict)))) # literally counting the number of digits in the length

	# optionally make dict for conversion vector file
	conversiondict = {} if args.conversion else None
	counter = 0

	# print example header
	if args.only_reorder:
		sys.stderr.write("# Ordering sequences by length, without renaming\n")
	else:
		if args.omit_zero:
			contignum = str(counter)
		else:
			contignum = str("{:0"+reclog+"}").format(counter) # format number with zeroes
			outputid = [args.name, contignum]
			if args.length:
				outputid.extend( [ args.length, "000" ] )
		sys.stderr.write("# Ordering sequences by length, renaming as {}\n".format(args.delimiter.join(outputid)) )

	# do the renaming and printing of contigs
	for seqrec in sorted(contigdict.values(), key=lambda x: len(x.seq), reverse=True):
		counter += 1 # first contig is 1, not 0
		if args.only_reorder:
			wayout.write(seqrec.format(args.format))
		else:
			if args.omit_zero:
				contignum = str(counter)
			else:
				contignum = str("{:0"+reclog+"}").format(counter) # format number with zeroes
			outputid = [args.name, contignum]
			if args.length:
				outputid.extend( [ args.length, str(len(seqrec.seq)) ] )
			seqrec.description = seqrec.id if args.old else ""
			if conversiondict is not None:
				conversiondict[str(seqrec.id)] = args.delimiter.join(outputid)
			seqrec.id = args.delimiter.join(outputid)
			wayout.write(seqrec.format(args.format))

	# make conversion vector file
	if args.conversion:
		sys.stderr.write("# Writing conversion file {}\n".format(args.conversion) )
		with open(args.conversion,'w') as cf:
			for k,v in conversiondict.items():
				cf.write("{}\t{}\n".format(k,v) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)

