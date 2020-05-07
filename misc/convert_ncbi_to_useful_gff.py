#!/usr/bin/env python
# convert_ncbi_to_useful_gff.py v1 created 2020-05-06

'''convert_ncbi_to_useful_gff.py  last modified 2020-05-06

    convert default NCBI bacterial genome GFF to a more useful format
    input can be .gz

convert_ncbi_to_useful_gff.py GCF_000011805.1_ASM1180v1_genomic.gff > GCF_000011805.1_ASM1180v1_genomic.useful.gff
'''

import sys
import gzip

if len(sys.argv)<2:
	sys.exit(__doc__)
else:
	genebuffer = {} # key is ID, value is line
	pseudo_counter = {} # key is CDS ID, value is True
	gff_file = sys.argv[1]
	if gff_file.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading GFF file {} as gzipped\n".format(gff_file) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading GFF file {}\n".format(gff_file) )
	# parse GFF
	for line in opentype(gff_file,'rt'):
		if line.strip() and line[0]!="#":
			lsplits = line.strip().split("\t")
			feature = lsplits[2]
			attributes = lsplits[8]
			attrd = dict([(field.strip().split("=")) for field in attributes.split(";") if field.count("=")])
			if feature=="region" or feature=="pseudogene": # ignore region features and pseudogenes
				continue
			elif feature=="gene":
				gene_id = attrd.get("ID")
				biotype = attrd.get("gene_biotype")
				if biotype=="protein_coding":
					genebuffer[gene_id] = line
				else:
					gene_line = line.strip() + ";product={}\n".format(biotype)
					sys.stdout.write(gene_line)
			elif feature=="CDS":
				gene_id = attrd.get("Parent")
				if attrd.get("pseudo","x")=="true":
					pseudo_counter[gene_id] = True
					continue
				gene_product = attrd.get("product")
				gene_line = genebuffer.get(gene_id,"").strip()
				gene_line = gene_line + ";product={}\n".format(gene_product)
				sys.stdout.write(gene_line)
				sys.stdout.write(line)
			else: 
				sys.stdout.write(line)
	sys.stderr.write("# Converted {} genes\n".format( len(genebuffer) ) )
	if pseudo_counter:
		sys.stderr.write("# Ignored {} pseudogenes\n".format( len(pseudo_counter) ) )

#

