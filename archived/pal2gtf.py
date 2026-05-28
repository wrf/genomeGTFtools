#!/usr/bin/env python
#
# pal2gtf.py v1.0 2016-01-19
#
# for SOFA terms:
# https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo

'''pal2gtf.py  last modified 2016-01-19

pal2gtf.py mt_dna.pal > mt_dna.gtf

    for gtf attributes, sequences appear as:
  ID=123.palindrome.n3;Note=GAGCCAC
    where:
    123 is the unique number of that repeat for the gtf file
    n3 is the copy number of that specific repeat
    palindrome can be changed to another tag using -i
    if the repeat is detected (e.g. H1) then this name is used instead:
  ID=123.H1.n3;Note=GAGCCAC

    generate .pal format files at:
http://emboss.bioinformatics.nl/cgi-bin/emboss/palindrome
'''

import sys
import argparse
import time

# pal output of EMBOSS explorer looks like:
# Palindromes:
# 33       ctaggagcccctag       46
#          ||||||||||||||
# 60       gatcctcggggatc       47

# repeats from Lavrov et al. 2012. Gene 505:1
repeatindex = {"CCACCTAG":"H1", 
               "CCAGGAT":"H2", "ACCAGGA":"H2", "ACCAGGAT":"H2", 
               "TTAGGCAT":"H3", "TTGGGCAT":"H3", 
               "ATTAATATATCCTATGA":"H4", "ATATCCTATGA":"H4", 
               "GAGCCAC":"H5", "TAGCCAC":"H5", "GAGCCACCT":"H5", 
               "TAGGCATC":"H6", 
               "GACCGAC":"H7", "GGACCGAC":"H7",
               "GGGATATCC":"H8", 
               "TAGGATAT":"H9"}

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input', type = argparse.FileType('rU'), default = '-', help="pal format file")
	parser.add_argument('-i','--identifier', help="tag for ID attribute [palindrome]", default="palindrome")	
	parser.add_argument('-p','--program', help="program for 2nd column in output [EMBOSS]", default="EMBOSS")
	parser.add_argument('-t','--type', help="gff type [inverted_repeat]", default="inverted_repeat")
	parser.add_argument('-v','--verbose', action="store_true", help="extra output")
	args = parser.parse_args(argv)

	startdromes = False # after header information, start getting palindromes
	getend = False # switch to get end position instead of start position
	genename = args.input.name.rsplit(".",1)[0] # by default, if not found in the file, use the file name

	linecounter = 0
	repeatcounter = 0
	repcounter = {}
	print >> sys.stderr, "# Starting parsing on {}".format(args.input.name), time.asctime()
	for line in open(sys.argv[1],'r').readlines():
		line = line.strip()
		if line: # skip empty lines
			linecounter += 1
			if line.find("Palindromes of:")==0: # get scaffold name from "Palindromes of:  mt-DNA_final"
				genename = line.split(":")[1].strip()
			if startdromes:
				lsplits = line.split()
				if getend:
					endpos = lsplits[0]
					if endpos[0]=="|": # skip alignment line
						continue
					repeatcounter += 1
					print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}\t.\t.\tID={}.{}.n{};Note={}".format(genename, args.program, args.type, startpos, endpos, dromelen, repeatcounter, repeatindex.get(palindrome, args.identifier), repcounter[palindrome], palindrome)
					getend = False
				else:
					startpos = lsplits[0]
					palindrome = lsplits[1].upper()
					dromelen = len(palindrome)
					repcounter[palindrome] = repcounter.get(palindrome,0) + 1
					getend = True
			if line.find("Palindromes:")==0: # indicating start of palindrome list
				startdromes = True
	print >> sys.stderr, "# Parsed {} lines".format(linecounter), time.asctime()
	print >> sys.stderr, "# Found {} repeats".format(repeatcounter), time.asctime()
	if args.verbose: # print table of repeats and occurrence
		for k,v in sorted(repcounter.items(), key=lambda x: x[1], reverse=True):
			print >> sys.stderr, k, v, repeatindex.get(k,"ID")

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
