#!/usr/bin/env python

'''stringtie_gtf_to_gff3.py  last modified 2023-01-27

convert stringtie GTF to standard GFF (ID-Parent format)
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

	sys.stderr.write( "# Converting {} to GFF\n".format(args.gtf) )
	for line in open(args.gtf,'r'):
		line = line.strip()
		if line and line[0] != "#":
			# stringtie format "transcript" features:
			# gene_id "B1_LR.1"; transcript_id "B1_LR.1.1"; cov "5.364951"; FPKM "0.929265"; TPM "2.592418";

			# stringtie format "exon" features:
			# gene_id "B1_LR.1"; transcript_id "B1_LR.1.1"; exon_number "1"; cov "5.649920";
			lsplits = line.split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]
			if attributes.find("ID=")>-1:
				raise ValueError("ERROR: ID tag found, file may already by GFF format:\n{}\n".format(line) )
			attrd = dict([(field.strip().split(" ")) for field in attributes.split(";") if field])
			if feature=="transcript" or feature=="mRNA": # changed to ID and Name
				lsplits[2] = args.transcript
				transcriptid = attrd.get("transcript_id").replace('"','')
				# if using StringTie with FPKM etc, then add transfer those as well
				fpkm = attrd.get("FPKM",None)
				tpm = attrd.get("TPM",None)
				coverage = attrd.get("cov",None)
				newattrs = "ID={0};Name={0};".format(transcriptid)
				if coverage is not None:
					newattrs += "cov={};".format(coverage.replace('"','') )
				if fpkm is not None:
					newattrs += "FPKM={};".format(fpkm.replace('"','') )
				if tpm is not None:
					newattrs += "TPM={};".format(tpm.replace('"','') )
			elif feature=="exon": # change to ID and Parent
				transcriptid = attrd.get("transcript_id").replace('"','')
				exonnum = attrd.get("exon_number","").replace('"','')
				if exonnum:
					newattrs = "ID={0}.exon{1};Parent={0}".format(transcriptid,exonnum)
				else:
					newattrs = "Parent={0}".format(transcriptid)
			lsplits[8] = newattrs
			print( "\t".join(lsplits) , file=wayout )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
