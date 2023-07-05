#!/usr/bin/env python
#

'''append_kegg_to_prodigal_gff.py  last modified 2023-07-05
    append KEGG annotations from blastKOALA to prodigal gff
    KEGG annotations will include gene name, description, and KEGG category
    and will create a gene feature for each CDS

    prodigal predicts proteins with option -a, but these are unannotated
      see:
      https://github.com/hyattpd/Prodigal
    
    KEGG annotations for those proteins can be easily acquired with blastKOALA
      see:
      https://www.kegg.jp/blastkoala/

append_kegg_to_prodigal_gff.py ko_annotation.txt prodigal.gff > prodigal_w_ko.gff
'''

import sys
import re

def clean_description(description_line):
	'''replace symbols and return a string'''
	symbolset = "[],/"
	for symbol in symbolset:
		description_line = description_line.replace(symbol,"_")
	return description_line

##Query	KO	Definition	Score	Second best
#CP025958.1_1 (250)					
#CP025958.1_2 (46)			1	K02346	1
#CP025958.1_3 (268)					
#CP025958.1_4 (45)			2	K01208	1
#CP025958.1_5 (360)			28	K03196	1
#CP025958.1_6 (500)	K19132	csb2; CRISPR-associated protein Csb2	22		
#CP025958.1_7 (376)	K19131	csb1; CRISPR-associated protein Csb1	184		
#CP025958.1_8 (280)	K19135	csx14; CRISPR-associated protein Csx14	13		
#CP025958.1_9 (860)	K07012	cas3; CRISPR-associated endonuclease/helicase Cas3 [EC:3.1.-.- 3.6.4.-]	162

if len(sys.argv)<2:
	sys.exit( __doc__ )
else:
	kegg_desc = {} # key is gene ID for GFF, value is KEGG definition
	kegg_genes = {} # key is gene ID for GFF, value is gene name of match
	kegg_cat = {} # key is gene ID, value is K number
	kegg_EC = {} # key is gene ID, value is optional EC number
	linecounter = 0
	full_description = None	
	sys.stderr.write("# reading annotations from {}\n".format( sys.argv[1] ) )
	for line in open(sys.argv[1],'r'):
		linecounter += 1
		if line and line[0]!="#":
			lsplits = line.split("\t")

			# attempt to extract from table
			# take CP025958.1_1 from CP025958.1_1 (250)
			try:
				protid = lsplits[0].split(" ")[0]
				ko_number = lsplits[1]
				genedesc = lsplits[2]
			except IndexError:
				continue # skip line, as it means duplicated function

			# try to parse description
			if genedesc:
				dsplits = genedesc.split(";")
				if len(dsplits)>1:
					# csb2; CRISPR-associated protein Csb2
					geneid, full_description = dsplits
				else: # for cases like:
					# metallo-beta-lactamase family protein
					geneid = None
					full_description = dsplits[0]
				# remove symbols from description
				full_description = clean_description( full_description.strip() )
				kegg_genes[protid] = geneid
				kegg_desc[protid] = full_description
				kegg_cat[protid] = lsplits[1]
			# if there is no description in the 3rd column, then check for a secondary functional prediction in 5th column
			elif lsplits[4].strip():
				full_description = clean_description( lsplits[4].strip() )
				kegg_desc[protid] = full_description
				kegg_cat[protid] = full_description
	sys.stderr.write("# found {} gene IDs on {} lines\n".format( len(kegg_desc),linecounter ) )

	#CP025958.1	Prodigal_v2.6.3	CDS	6916	7665	2.5	+	0	ID=1_1;partial=00;start_type=GTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.648;conf=64.14;score=2.53;cscore=6.01;sscore=-3.48;rscore=2.21;uscore=-0.38;tscore=-5.31;
	#CP025958.1	Prodigal_v2.6.3	CDS	7662	7799	2.8	+	0	ID=1_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.688;conf=65.42;score=2.77;cscore=-1.89;sscore=4.66;rscore=-9.81;uscore=0.49;tscore=3.38;
	#CP025958.1	Prodigal_v2.6.3	CDS	8912	9715	14.2	+	0	ID=1_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.655;conf=96.29;score=14.17;cscore=6.62;sscore=7.54;rscore=-5.30;uscore=-0.44;tscore=6.25;
	#CP025958.1	Prodigal_v2.6.3	CDS	13477	13611	11.0	-	0	ID=1_4;partial=00;start_type=ATG;rbs_motif=GGxGG;rbs_spacer=5-10bp;gc_cont=0.667;conf=92.65;score=11.02;cscore=11.80;sscore=-0.78;rscore=-3.54;uscore=0.11;tscore=3.30;
	#CP025958.1	Prodigal_v2.6.3	CDS	14057	15136	39.3	-	0	ID=1_5;partial=00;start_type=GTG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;gc_cont=0.623;conf=99.99;score=39.28;cscore=44.09;sscore=-4.80;rscore=0.04;uscore=-0.37;tscore=-5.31;

	sys.stderr.write("# reading GFF from {}\n".format( sys.argv[2] ) )
	for line in open(sys.argv[2],'r'):
		if line and line[0]!="#":
			lsplits = line.split("\t")
			scaffold = lsplits[0]
			feature = lsplits[2]
			attributes = lsplits[8].strip()
			prodigal_geneid = re.search('ID=([\w\d._|-]+);', attributes).group(1)
			geneid = "{}_{}".format( scaffold, prodigal_geneid.split("_")[1] ) # gene ID invented as scaffold_prodigal-number
			parentattrs = "ID={};Name={};gene_biotype=protein_coding;".format(geneid, kegg_genes.get(geneid,"None"))
			gene_name = kegg_genes.get(geneid,None)
			if gene_name is not None:
				gene_name = gene_name.split(',')[0]
				parentattrs += "gene={};".format(gene_name.replace('-','')) # for names like RP-L15, MRPL15, rplO
			genedesc = kegg_desc.get(geneid,None)
			if genedesc is not None:
				parentattrs += "Description={};".format(genedesc)
			geneacc = kegg_cat.get(geneid,None)
			if geneacc is not None:
				parentattrs += "Accession={};".format(geneacc)
			parentline = "{0}\t{1}\tgene\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format( lsplits[0], lsplits[1], lsplits[3], lsplits[4], lsplits[5], lsplits[6], lsplits[7], parentattrs )
			sys.stdout.write(parentline)
			cds_attr_list = attributes.split(";")
			cdsattrs = "Parent={};{}".format( geneid, ";".join(cds_attr_list[1:]) )
			cdsline = "{0}\t{1}\tCDS\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format( lsplits[0], lsplits[1], lsplits[3], lsplits[4], lsplits[5], lsplits[6], lsplits[7], cdsattrs )
			sys.stdout.write(cdsline)
	sys.stderr.write("# done parsing GFF\n")

	#CP025958.1	Prodigal_v2.6.3	gene	6916	7665	2.5	+	0	ID=1_1;Name=None
	#CP025958.1	Prodigal_v2.6.3	CDS	6916	7665	2.5	+	0	Parent=1_1;partial=00;start_type=GTG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.648;conf=64.14;score=2.53;cscore=6.01;sscore=-3.48;rscore=2.21;uscore=-0.38;tscore=-5.31;
	#CP025958.1	Prodigal_v2.6.3	gene	7662	7799	2.8	+	0	ID=1_2;Name=None;Description=K02346;Accession=K02346
	#CP025958.1	Prodigal_v2.6.3	CDS	7662	7799	2.8	+	0	Parent=1_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.688;conf=65.42;score=2.77;cscore=-1.89;sscore=4.66;rscore=-9.81;uscore=0.49;tscore=3.38;
	#CP025958.1	Prodigal_v2.6.3	gene	8912	9715	14.2	+	0	ID=1_3;Name=None
	#CP025958.1	Prodigal_v2.6.3	CDS	8912	9715	14.2	+	0	Parent=1_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.655;conf=96.29;score=14.17;cscore=6.62;sscore=7.54;rscore=-5.30;uscore=-0.44;tscore=6.25;
	#CP025958.1	Prodigal_v2.6.3	gene	13477	13611	11.0	-	0	ID=1_4;Name=None;Description=K01208;Accession=K01208
	#CP025958.1	Prodigal_v2.6.3	CDS	13477	13611	11.0	-	0	Parent=1_4;partial=00;start_type=ATG;rbs_motif=GGxGG;rbs_spacer=5-10bp;gc_cont=0.667;conf=92.65;score=11.02;cscore=11.80;sscore=-0.78;rscore=-3.54;uscore=0.11;tscore=3.30;
	#CP025958.1	Prodigal_v2.6.3	gene	14057	15136	39.3	-	0	ID=1_5;Name=None;Description=K03196;Accession=K03196
	#CP025958.1	Prodigal_v2.6.3	CDS	14057	15136	39.3	-	0	Parent=1_5;partial=00;start_type=GTG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;gc_cont=0.623;conf=99.99;score=39.28;cscore=44.09;sscore=-4.80;rscore=0.04;uscore=-0.37;tscore=-5.31;





