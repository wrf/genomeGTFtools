#!/usr/bin/env python
#
# keg_output_to_table.py

"""keg_output_to_table.py  last modified 2023-07-04

    convert KO links to tabular format

keg_output_to_table.py ko00001.keg > ko00001.tab

    data retrieved from:
https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=

"""

color_description = """Functional category (used in global pathway maps and genome maps)
09101 Carbohydrate metabolism	      	#0000ee
09102 Energy metabolism	      	#9933cc
09103 Lipid metabolism	      	#009999
09104 Nucleotide metabolism	      	#ff0000
09105 Amino acid metabolism	      	#ff9933
09106 Metabolism of other amino acids	      	#ff6600
09107 Glycan biosynthesis and metabolism	      	#3399ff
09108 Metabolism of cofactors and vitamins	      	#ff6699
09109 Metabolism of terpenoids and polyketides	      	#00cc33
09110 Biosynthesis of other secondary metabolites	      	#cc3366
09111 Xenobiotics biodegradation and metabolism	      	#ccaa99
09120 Genetic information processing	      	#ffcccc
09130 Environmental information processing	      	#ffff00
09140 Cellular processes	      	#99cc66
09150 Organismal systems	      	#99cc66
09160 Human diseases	      	#99cc66
09181 Protein families: metabolism	      	#ccffff
09182 Protein families: genetic information processing	      	#ffcccc
09183 Protein families: signaling and cellular processes	      	#99cc66
09191 Unclassified: metabolism	      	#ccffff
09192 Unclassified: genetic information processing	      	#ffcccc
09193 Unclassified: signaling and cellular processes	      	#99cc66
"""

import sys


if len(sys.argv) < 2:
	sys.exit( __doc__ )
else:
	color_index = { "09101": "#0000ee" , "09102": "#9933cc" , "09103": "#009999" , "09104":"#ff0000" , "09105": "#ff9933" , "09106": "#ff6600" , "09107": "#3399ff" , "09108": "#ff6699" , "09109": "#00cc33" , "09110": "#cc3366" , "09111": "#ccaa99" , "09120": "#ffcccc" , "09130": "#ffff00" , "09140": "#99cc66" , "09150": "#99cc66" , "09160": "#99cc66" , "09181": "#ccffff" , "09182": "#ffcccc" , "09183": "#99cc66" , "09191": "#ccffff" , "09192": "#ffcccc" , "09193": "#99cc66" }
	linecounter = 0
	writecounter = 0
	sys.stderr.write("# reading {}\n".format( sys.argv[1], writecounter ) )
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line:
			first_char = line[0]
			if first_char in "#!+A": # skip 
				continue
			if first_char in "BCD":
				lsplits = line.split(None,2)
				if len(lsplits) > 1:
					linecounter += 1
					k_code = lsplits[1]
					description = lsplits[2]
					if first_char=='B':
						top_level_code = k_code
						top_level_desc = description
						color_code = color_index.get(top_level_code)
					elif first_char=='C':
						pathway_code = k_code
						pathway_desc = description
					elif first_char=='D':
						outline = "{}\t{}\t{}\t{}\t{}\n".format( top_level_code, color_code, pathway_code, k_code, description )
						sys.stdout.write( outline )
						writecounter += 1
					else: # should never happen
						pass
	sys.stderr.write("# found {} lines, wrote {} KEGG entries\n".format( linecounter, writecounter ) )

#
