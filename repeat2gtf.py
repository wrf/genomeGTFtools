#! /usr/bin/env python
#
# repeat2gtf.py v1 created 2016-02-18
# v1.1 2023-11-06 allow gzip input 
#
# for SOFA terms:
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo

"""repeat2gtf.py  v1.1 last modified 2023-11-06
    generates a GFF3 format file of repeats, typically Ns as gaps
  the script only searches the FORWARD strand, meaning would need
  to be run twice for non-palindromic sequences (e.g. CACA vs GTGT)
    Requires Bio library (biopython)

repeat2gtf.py scaffolds.fasta > scaffolds_gaps.gff

    for gzip files, use stdin as

gzip -dc scaffolds.fasta.gz | repeat2gtf.py -

    for non-gap repeats (anything other than N)
    change -t to direct_repeat

    GFF output contains 9 tab-separated columns for:
contig  program  type  start  end  score (length)  strand  phase  attributes
    where attributes are composed of:

  identifier.repeat.unique-number.length

    a typical line will appear as:
scaffold123  assembler  gap  2001  2105  105  .  .  ID=gap.N.1.105

    meaning gap (default name), of N, number 1, length is 105
"""

import sys
import argparse
import time
import os
import re
import gzip
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default = '-', help="fasta format file, can be .gz, or stdin as -")
	parser.add_argument('-a', '--above', type=int, metavar='N', default=2, help="only print repeats/gaps longer than N, in letters [2]")
	parser.add_argument('-b', '--below', type=int, metavar='N', default=1000000000, help="only print sequences repeats/gaps shorter than N, in letters")
	parser.add_argument('--format', metavar='fastq', default='fasta', help="import fastq format sequences")
	parser.add_argument('-i','--identifier', help="tag for ID attribute [gap]", default="gap")
	parser.add_argument('-l','--lowercase', action="store_true", help="search for lowercase letters as well")
	parser.add_argument('-p','--program', help="program for 2nd column in output [assembler]", default="assembler")
	parser.add_argument('-r','--repeat', metavar='N', default='N', help="measure length of the longest polyN, NX repeat, etc.")
	parser.add_argument('-t','--type', help="feature type for 3rd column [gap]", default="gap")
	parser.add_argument('--attribute', help="attribute for 9th column [ID]", default="ID")
	parser.add_argument('-v','--verbose', action="store_true", help="extra output")
	args = parser.parse_args(argv)

	# all integers initialized
	seqcount = 0
	combined_rep_len = 0
	repcounter = 0
	longestrepeat = 0
	longestrepcontig = ""

	# make the one or two regexps
	repeatregex = re.compile("({0})+".format(args.repeat) )
	if args.lowercase: # lowercase version of the same repeat
		lowerrepeat = args.repeat.lower()
		lcregex = re.compile("({0})+".format(lowerrepeat) )
		sys.stderr.write("# Parsing repeats of {} and {} from {}  {}\n".format(args.repeat, lowerrepeat, args.input_file.name, time.asctime() ) )
	else:
		sys.stderr.write("# Parsing repeats of {} from {}  {}\n".format(args.repeat, args.input_file.name, time.asctime() ) )

	# begin iterating through sequences, then search the regular expression
	if args.input_file.name.rsplit('.',1)[-1]=="gz":
		input_handler = gzip.open(args.input_file.name,'rt')
	else:
		input_handler = args.input_file
	for seqrec in SeqIO.parse(input_handler, args.format):
		seqcount += 1
		contig = seqrec.id
		reptracker = {} # key is position, value is list of start, end, length
		for rep in repeatregex.finditer(str(seqrec.seq)): # iterate through normal repeats
			replen = rep.end() - rep.start()
			if replen < args.above or replen > args.below:
				continue
			if replen > longestrepeat:
				longestrepeat = replen
				longestrepcontig = contig
			combined_rep_len += replen
			repstart = rep.start()+1
			repend = rep.end()
			reptracker[repstart] = [repstart, repend, replen, args.repeat]
		if args.lowercase:
			for rep in lcregex.finditer(str(seqrec.seq)): # iterate through lowercase repeats
				replen = rep.end() - rep.start()
				if replen < args.above or replen > args.below:
					continue
				combined_rep_len += replen
				repstart = rep.start()+1
				repend = rep.end()
				reptracker[repstart] = [repstart, repend, replen, lowerrepeat]
		for startpos in sorted(reptracker.keys()):
			repcounter += 1
			sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t.\t.\t{}={}.{}.{}.{}\n".format(contig, args.program, args.type, reptracker[startpos][0], reptracker[startpos][1], reptracker[startpos][2], args.attribute, args.identifier, reptracker[startpos][3], repcounter, reptracker[startpos][2]) )

	sys.stderr.write("# Counted {} sequences  {}\n".format(seqcount, time.asctime() ) )
	try:
		sys.stderr.write("# Counted {} repeats of {} total bases, average {:.2f}\n".format(repcounter, combined_rep_len, combined_rep_len*1.0/repcounter) )
	except ZeroDivisionError: # repcounter remains at 0, no repeats
		sys.stderr.write( "# No repeats found\n" )
	if longestrepeat:
		sys.stderr.write("# Longest repeat was {} bases on {}\n".format(longestrepeat, longestrepcontig) )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
