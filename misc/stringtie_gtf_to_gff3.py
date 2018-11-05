#!/usr/bin/env python

'''convert stringtie GTF to standard GFF (ID-Parent format)
stringtie_gtf_to_gff3.py stringtie.gtf > stringtie.gff3

   change transcript features with -t, such as
   -t mRNA
'''

import sys
import argparse

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('gtf', help="stringtie GTF file")
	parser.add_argument('--transcript', default="transcript", help="optional name for transcript features, default is transcript")
	args = parser.parse_args(argv)

	print >> sys.stderr, "# Converting {} to GFF".format(args.gtf)
	for line in open(args.gtf,'r'):
		line = line.strip()
		if line and line[0] != "#":
			lsplits = line.split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split(" ")) for field in attributes.split(";") if field])
			if feature=="transcript": # changed to ID and Name
				lsplits[2] = args.transcript
				transcriptid = attrd.get("transcript_id").replace('"','')
				newattrs = "ID={0};Name={0}".format(transcriptid)
			elif feature=="exon": # change to ID and Parent
				transcriptid = attrd.get("transcript_id").replace('"','')
				exonnum = attrd.get("exon_number").replace('"','')
				newattrs = "ID={0}.exon{1};Parent={0}".format(transcriptid,exonnum)
			lsplits[8] = newattrs
			print >> wayout, "\t".join(lsplits)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
