#!/usr/bin/env python
#
# 2020-09-03 previous version
# 2022-05-02 fix for old format without tags for transcript features

'''
  augustus_to_gff3.py  last modified  2022-05-02
  generate true GFF3 from augustus output
  involves creating mRNA features from transcript
  and exon from CDS
  also will formalize some ID-Parent tags

augustus_to_gff3.py augustus.gff > augustus_no_comments.gff3

  some formats include exons already, rather than needing to create them from CDS
  use -x

augustus_to_gff3.py braker2.gff -x > braker2_no_comments.gff3
'''

import sys

if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	if len(sys.argv)==3 and sys.argv[2]=="-x":
		exon_from_CDS = False
	else:
		exon_from_CDS = True

	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line and line[0] != "#":
			lsplits = line.split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]

			if feature=="gene": # attributes should be single item
				newattrs = "ID={0};Name={0}".format(attributes)
			elif feature=="transcript": # changed to ID and Name
				lsplits[2] = "mRNA"
				try: # dictionary update sequence element #0 has length 1; 2 is required
					attrd = dict([(field.strip().split(" ",1)) for field in attributes.split(";") if field])
				except ValueError: # catch for cases of old version with format of g4.t1 as the only line item
					attrd = {"gene_id":attributes.split(".",1)[0] , "transcript_id": attributes.strip() }
				if attrd.get("gene_id", False) and attrd.get("transcript_id", False):
					newattrs = "ID={0};Parent={1}".format(attrd.get("transcript_id").replace('"',''), attrd.get("gene_id").replace('"',''))
				else:
					newattrs = "ID={0};Parent={1};Name={0}".format(attributes, attributes.split(".")[0])
			elif feature=="CDS": # change to ID and Parent
				attrd = dict([(field.strip().split(" ",1)) for field in attributes.split(";") if field])
				transcriptid = attrd.get("transcript_id").replace('"','')
				newattrs = "ID={0}.cds;Parent={0}".format(transcriptid)
			elif feature=="stop_codon" or feature=="start_codon":
				attrd = dict([(field.strip().split(" ",1)) for field in attributes.split(";") if field])
				transcriptid = attrd.get("transcript_id").replace('"','')
				newattrs = "Parent={0}".format(transcriptid)
			else:
				continue
			lsplits[8] = newattrs
			if feature=="CDS" and exon_from_CDS:
				exonsplits = list(lsplits)
				exonsplits[2] = "exon"
				exonsplits[7] = "."
				exonsplits[8] = "ID={0}.exon;Parent={0}".format(transcriptid)
				sys.stdout.write("{}\n".format( "\t".join(exonsplits) ) )

			sys.stdout.write("{}\n".format( "\t".join(lsplits) ) )
