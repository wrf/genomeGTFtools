# generate quick histogram of synteny block sizes from microsynteny.py
#
# v1.0 2018-07-04 basic barplot with command line input
# v1.1 2023-06-16 make double plot

# called from command line as:
# Rscript ~/git/genomeGTFtools/synteny_block_length_plot.R dpulex_v_dmagna.blastp.microsynteny.tab 0.7
# output PDF should be automatically named
# additional argument is to specify color between 0 and 1

# block expected in format of 12 column tabular:
#Avas-scaffold_001	JAKMXF010000055.1	blk-1	Avas.s001.g76.i1	289852	296181	-	LOD99_14495.mRNA.1	174868	179817	-	313.0
#Avas-scaffold_001	JAKMXF010000055.1	blk-1	Avas.s001.g77.i1	309436	315869	+	LOD99_14496.mRNA.1	179847	182507	-	101.0
#Avas-scaffold_001	JAKMXF010000055.1	blk-1	Avas.s001.g78.i1	316567	319320	+	LOD99_14496.mRNA.1	179847	182507	-	105.0

# columns are:
block_headers = c("q_scaffold", "s_scaffold", "block_ID", 
                  "q_gene_ID", "q_gene_pos_start", "q_gene_pos_end", "q_gene_strand", 
                  "s_gene_ID", "s_gene_pos_start", "s_gene_pos_end", "s_gene_strand", 
                  "bitscore" )

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/genomeGTFtools/test_data/acropora_vs_styllophora_microsynteny.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz$","",inputfile,perl=TRUE),perl=TRUE)
if (inputfile==outputfile) { stop("cannot parse input file to generate output file name, add a unique 3-letter suffix") }

print(paste("# Reading block info file", inputfile, ", writing PDF to", outputfile))
blockdata = read.table(inputfile, col.names = block_headers, sep="\t")
#head(blockdata)

counttable = table(blockdata[,3])
totalblk = length(counttable)
totalgenes = length(unique(blockdata[,4]))

# define function, and count unique genes in each block
# for comparison of block length, fragmentation, etc
unilen = function(x){ return( length(unique(x)) ) }
len_by_block = aggregate(blockdata[,4], list(blockdata[,3]), length )
print("# Calculating unique genes per block")
q_genes_by_block = aggregate(blockdata[,4], list(blockdata[,3]), unilen )
s_genes_by_block = aggregate(blockdata[,8], list(blockdata[,3]), unilen )

# define point color
if ( !is.na(args[2]) ) {
  # use H value specified by user, from 0 to 1
  if ( as.numeric(args[2]) >= 0 & as.numeric(args[2]) <= 1) {
    user_hue = as.numeric(args[2])
  } else { # meaning user entered something outside 0-1
    user_hue = 0.39 # default green
  }
} else {
  # value is not given at all in command line
  user_hue = 0.39 # default green
}

point_color = rainbow(1, s = 0.8, v = 0.3, start = user_hue, alpha = 0.1)
bar_color = rainbow(1, s = 0.8, v = 0.7, start = user_hue, alpha = 0.9)
#barcolor = c("#76ee94") # green
#barcolor = c("#fc8d62") # red/orange
#barcolor = c("#80b1d3") # blue
print( paste("# Using color",point_color,"for points") )

# generate PDF with both graphs
print( paste("# Making PDF of two graphs", outputfile) )
pdf(outputfile, width=8, height=11)
par(mfrow=c(2,1), mar=c(4.2,5,2,2))

# scatterplot of genes per block in different genomes
plot(s_genes_by_block$x, q_genes_by_block$x, cex.axis=1.3, cex.lab=1.4,
     main=basename(inputfile), 
     xlab = "Genes within block for subject species", ylab = "Genes within block for query species",
     pch=16, cex = 3, col=point_color )
abline(a = 0,b = 1, lwd=2, lty=2, col="#00000055")
text(min(q_genes_by_block$x), max(q_genes_by_block$x) , 
     paste(length(unique(blockdata[,4])),"total query genes"), pos = 4 )
text(max(s_genes_by_block$x) , min(s_genes_by_block$x), 
     paste(length(unique(blockdata[,8])),"total subject genes"), pos = 2 )

# make settings for barplot
longestblock = max(counttable)
lencap = min(longestblock, 20) # take whichever is smaller, longest block or 20
counttable[counttable>lencap] = lencap
counthist = hist(counttable, breaks=seq(2.5,longestblock+0.5,1), 
                 axes=FALSE, plot=TRUE, xlim=c(3,lencap), 
                 main=basename(inputfile), 
                 xlab="Number of genes in block", ylab="Number of blocks", 
                 cex.lab=1.5, col=bar_color )
mostcommon = max(counthist$counts)
axis(1,at=seq(3,lencap,3), labels=seq(3,lencap,3), cex.axis=1.2 )
axis(2, cex.axis=1.3)
text(counthist$breaks+0.5, counthist$counts+(0.025*mostcommon), counthist$counts)
totallab = paste(totalblk, "total blocks", sep=" ")
genelab = paste(totalgenes, "total genes in blocks", sep=" ")
text(lencap, mostcommon*0.65, totallab, cex=1.5, pos=2)
text(lencap, mostcommon*0.50, genelab, cex=1.5, pos=2)
dev.off()


#