#!/usr/bin/env Rscript
# convert GFF annotation format to polygon blocks for an entire bacterial genome
# for bacteria, genes are usually entire CDS and are better rendered with polygons
# draw style of bacterial genomes, meaning polygons in a row
# meaning like below, all on the same line
# |||> ||> ||> |||||> <||
#
# last modified 2023-07-05

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/aliivibrio_fisheri_PROK/GCF_000011805.1_ASM1180v1_genomic.gff"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz$","",inputfile,perl=TRUE),perl=TRUE)
if (inputfile==outputfile) { stop("cannot parse input file to generate output file name, add a unique 3-letter suffix") }

genome_gff = read.table(inputfile, header=FALSE, sep="\t", quote='', stringsAsFactors=FALSE)

# explicitly use "region" features, if they are given
chr_regions = genome_gff[genome_gff$V3=="region",]
chr_names = chr_regions$V1
# otherwise get from scaffold column
if (nrow(chr_regions)==0){ chr_names = unique( genome_gff$V1 ) }
n_chrs = length(chr_names)

feature_counts = table(genome_gff$V3)
# for GenBank format, this will take CDS mRNA tRNA pseudogene, and others
# for de novo annotations, like with prodigal, only CDS are present
feature_type = as.character( ifelse( !is.na(feature_counts["gene"]) , "gene" , "CDS") )

# default page, mostly for debugging
p=1;i=1;

offset_width = 500 # bp, minimum gene size to draw polygon, otherwise makes triangle
offset_height = 1 # relates to polygon height

use_kegg_index = !is.na(args[2]) # can type anything
# this file generated with another script:
# keg_output_to_table.py ko00001.keg.gz | gzip > ko00001.tab.gz
if (use_kegg_index){
  kegg_index_table = read.table("~/git/genomeGTFtools/test_data/ko00001.tab.gz", header=FALSE, sep="\t", stringsAsFactors = FALSE, quote = "", comment.char = "")
}
# colors from https://www.genome.jp/kegg/kegg1c.html
kegg_cats = c("Carbohydrate", "Energy", "Lipid", "Nucleotide", "Amino acid", "Other AA",
              "Glycan", "Cofactors", "Terpenoid/PKS", "Secondary", "Xenobiotic",
              "DNA/RNA", "Environment", "Cellular", "Other")
kegg_colors = c("#0000ee", "#9933cc", "#009999", "#ff0000", "#ff9933", "#ff6600", 
                "#3399ff", "#ff6699", "#00cc33", "#cc3366", "#ccaa99", "#ffcccc",
                "#ffff00", "#99cc66", "#888888")

##############
# draw the PDF
pdf(file=outputfile, width=8, height=11.5, paper="a4")
par(mar=c(1,4,1,5.5))

# for each chromosome
for (i in 1:n_chrs){
  current_chr = chr_names[i]
  cds_features = genome_gff[(genome_gff$V1==current_chr & genome_gff$V3==feature_type),]
  max_genome_size = chr_regions[i,5]
  # if there are no region features, take the last feature as max length
  if (is.na(max_genome_size)){max_genome_size = max(cds_features$V5)}
  n_pages = ceiling(max_genome_size/1000000)
  
  # for each 1 million bp, make a separate page
  for (p in 1:n_pages){
    # plot 1 million bp per page
    page_offset = (p-1)*1000000
    cds_page_num = ceiling(cds_features$V4/1000000)
    cds_on_page = cds_features[cds_page_num==p,]
    # get gene names assuming as attribute gene=
    gene_names = gsub("^.*gene=(\\w+);.*$","\\1",cds_on_page$V9)
    no_gene_name_index = grep("ID=",gene_names)
    gene_names[no_gene_name_index] = NA
    
    if (use_kegg_index){
      ko_accessions = gsub("^.*Accession=(\\w+);.*$","\\1",cds_on_page$V9)
      acc_to_ko_index = match(ko_accessions, kegg_index_table$V4)
      has_accession = which(!is.na(acc_to_ko_index))
      kegg_color = kegg_index_table$V2[acc_to_ko_index]
    }

    # adjust numbers of start and end positions depending on page
    cds_starts = cds_on_page$V4 - page_offset
    cds_ends = cds_on_page$V5 - page_offset
    cds_strands = cds_on_page$V7
    is_forward = which(cds_strands=="+")
    is_reverse = which(cds_strands=="-")
    is_strandless = which(cds_strands==".")

    # color rRNA, tRNA, and proteins differently
    # tx_type works for NCBI GFFs where gene_biotype is given
    tx_type = gsub("^.*gene_biotype=(\\w+);.*$","\\1",cds_on_page$V9)
    tx_colors = rep("#eefFee", nrow(cds_on_page) )
    tx_colors[tx_type=="rRNA"] = "#025a8d"
    tx_colors[tx_type=="tRNA"] = "#881a25"
    tx_colors[(tx_type=="protein_coding" & cds_strands=="+")] = "#7b7b7b"
    tx_colors[(tx_type=="protein_coding" & cds_strands=="-")] = "#cecece"
    if (use_kegg_index){ tx_colors[has_accession] = kegg_color[has_accession] }
    # otherwise assume that features are CDS
    tx_colors[(cds_on_page$V3=="CDS" & cds_strands=="+")] = "#2b693a"
    tx_colors[(cds_on_page$V3=="CDS" & cds_strands=="-")] = "#6cca87"

    # adjustments for position on page
    cds_x_offset = floor((cds_starts)/50000)*50000
    cds_y_index = 105-(ceiling((cds_starts)/50000)*5)

    plot(0,0, type='n', axes=FALSE, frame.plot=FALSE, 
         ylim=c(0,100), xlim=c(0,50200), 
         xlab="", ylab="",
         main=paste(basename(inputfile),"-",current_chr,"-",p))
    lh_labels = rev((seq(0,19,1)*50000+1)+page_offset)
    rh_labels = as.integer(rev(seq(1,20,1)*50000)+page_offset)
    is_fit_onto_page = lh_labels <= max_genome_size
    axis(2, at=seq(5,100,5)[is_fit_onto_page], labels=lh_labels[is_fit_onto_page], cex.axis=0.9, tick = FALSE, las=1)
    axis(4, at=seq(5,100,5)[is_fit_onto_page], labels=rh_labels[is_fit_onto_page], cex.axis=0.9, tick = FALSE, las=1)
    
    # draw forward polygons
    for (g in is_forward) {
      genelen = cds_ends[g]-cds_starts[g]
      used_offset = ifelse(genelen < offset_width, genelen, offset_width)
      x_offset = cds_x_offset[g]
      forward_x = c( cds_ends[g], cds_ends[g]-used_offset, cds_starts[g], cds_starts[g], cds_ends[g]-used_offset, cds_ends[g])
      forward_y = c( cds_y_index[g], cds_y_index[g]-offset_height, cds_y_index[g]-offset_height, cds_y_index[g]+offset_height, cds_y_index[g]+offset_height, cds_y_index[g])
      polygon( forward_x-x_offset, forward_y, col=tx_colors[g] )
    }
    # draw reverse polygons
    for (g in is_reverse) {
      genelen = cds_ends[g]-cds_starts[g]
      used_offset = ifelse(genelen < offset_width, genelen, offset_width)
      reverse_x = c( cds_starts[g], cds_starts[g]+used_offset, cds_ends[g], cds_ends[g], cds_starts[g]+used_offset, cds_starts[g])
      x_offset = cds_x_offset[g]
      reverse_y = c( cds_y_index[g], cds_y_index[g]-offset_height, cds_y_index[g]-offset_height, cds_y_index[g]+offset_height, cds_y_index[g]+offset_height, cds_y_index[g])
      polygon( reverse_x-x_offset, reverse_y, col=tx_colors[g] )
    }
    # draw rectangles for strandless features, usually this is an error
    for (g in is_strandless) {
      x_offset = cds_x_offset[g]
      rect(cds_starts[g], cds_y_index[g]-offset_height, 
           cds_ends[g], cds_y_index[g]+offset_height, col=tx_colors[g] )
    }
    text(cds_starts-cds_x_offset+offset_width, cds_y_index-offset_height, gene_names, cex=0.5, adj=c(1,1), srt=45)
    if (use_kegg_index){ legend(-5000,min(cds_y_index)-4, legend = kegg_cats,  bty = 'n', ncol = 5, 
           col = kegg_colors, pch = 15, xpd = TRUE ) }
  }
}
dev.off()

#