#!/usr/bin/env python
# convert_embl_to_useful_gff.py v1 created 2020-05-27

'''convert_embl_to_useful_gff.py  last modified 2020-05-27

    convert default NCBI or EMBL eukaryotic genome GFF to a more useful format
    input can be .gz

convert_embl_to_useful_gff.py GCA_900157425.1_version_2_genomic.gff > GCA_900157425.1_version_2_genomic.useful.gff
'''

import sys
import gzip

if len(sys.argv)<2:
	sys.exit(__doc__)
else:
	gene_to_tx_code = {} # key is gene ID, value is modified transcript ID
	#linebuffer = [] # list of lines, reset with each gene
	cds_counter = 0
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
			elif feature=="intron":
				continue
			elif feature=="gene":
			#	if linebuffer:
			#		sys.stdout.write("\n".join(linebuffer) )
			#		linebuffer = [] # reset buffer
				gene_id = attrd.get("ID")
				biotype = attrd.get("gene_biotype")
				if biotype=="protein_coding":
					tx_splits = line.strip().split("\t")
					tx_splits[2] = "mRNA"
					tx_id = gene_id.replace("gene","rna")
					tx_attributes = "ID={1};Parent={0};Name={1};gene_biotype=protein_coding;locus_tag={0}\n".format(gene_id, tx_id)
					tx_splits[8] = tx_attributes
				#	linebuffer.append(line)
					sys.stdout.write(line)
					tx_line = "\t".join(tx_splits)
				#	linebuffer.append(tx_line)
					sys.stdout.write(tx_line)
					gene_to_tx_code[gene_id] = tx_id
				else:
					gene_line = line.strip() + ";product={}\n".format(biotype)
					sys.stdout.write(gene_line)
			elif feature=="exon" or feature=="CDS":
				gene_id = attrd.get("Parent")
				if attrd.get("pseudo","x")=="true":
					pseudo_counter[gene_id] = True
					continue
				gene_to_tx = gene_to_tx_code.get(gene_id,None)
				if gene_to_tx: # otherwise should be tRNA
					cds_counter += 1
					fixed_attrs = lsplits[8].replace(gene_id, gene_to_tx_code[gene_id])
					lsplits[8] = fixed_attrs
				cds_line = "{}\n".format( "\t".join(lsplits) )
				sys.stdout.write(cds_line)
			else: 
				sys.stdout.write(line)
	sys.stderr.write("# Converted {} genes, with {} CDS features\n".format( len(gene_to_tx_code), cds_counter ) )
	if pseudo_counter:
		sys.stderr.write("# Ignored {} pseudogenes\n".format( len(pseudo_counter) ) )

#

