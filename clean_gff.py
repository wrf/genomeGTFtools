#! /usr/bin/env python
# cleangff.py modified from reformatgmapgff.py

'''clean_gff.py  last modified 2019-09-26

  operation will:
    remove comment lines
    fix any erroneous start/end positions, start must always be less than end
    rename non-standard features:
      5'-UTR, 3'-UTR, cds
      five_prime_UTR, three_prime_UTR, CDS

clean_gff.py gmap.gff > new.gff
'''

import sys

if len(sys.argv) < 2:
	sys.stderr.write( __doc__ )
else:
	infile = sys.argv[1]
	mode = "n"
	if len(sys.argv) > 2:
		mode = sys.argv[2].replace("-","")
	sys.stderr.write("# Reformatting {}\n".format(infile) )
	with open(infile) as go:
		for line in go:
			line = line.strip()
			line = line.split("#")[0] # this should remove all in line comments
			if line:
				lsplits = line.split("\t")
				if int(lsplits[3]) > int(lsplits[4]) and lsplits[6]=="-":
					lsplits[3], lsplits[4] = lsplits[4], lsplits[3]
				if lsplits[2]=="5'-UTR":
					lsplits[2] = "five_prime_UTR"
				elif lsplits[2]=="3'-UTR":
					lsplits[2] = "three_prime_UTR"
				elif lsplits[2]=="cds":
					lsplits[2] = "CDS"
				sys.stdout.write("{}\n".format( "\t".join(lsplits) ) )
					# Name=comp12345 is always second position from gmap
				#	trans_name = splits[8].split(";")[1].split("=")[1]
				#	contig = splits[0]
				#	target = "Target={}".format(contig)
				#	tc = [contig, trans_name]
				#	tc.extend(splits[2:])
				#	print >> sys.stdout, "{};{}".format(line.rstrip(), target)
	sys.stderr.write("# Done reformatting {}\n".format(infile) )

