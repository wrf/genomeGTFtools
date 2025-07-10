#!/usr/bin/env Rscript
# convert tabular annotation format to rectangle blocks
# or polygon blocks for single CDS genes, like bacteria
# last modified 2025-07-10 changed xpd to on by default

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/mnemiopsis_leidyi/ml0011_63k-88k_annot_short.tab"
#inputfile = "~/git/genomeGTFtools/test_data/lux_locus_short_annot.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

annottab = read.table(inputfile,header=FALSE,sep="\t",stringsAsFactors=FALSE)

num_tx = 30
seqcount = args[2]
if (!is.na(seqcount)) {
	num_tx = as.integer(seqcount)
}

# mode to draw style of bacterial genomes, meaning polygons in a row
# meaning like below, all on the same line
# |||> ||> ||> |||||> <||
bmode = FALSE

# default is to draw rectangles as exons with an arrow
# for bacteria, genes are usually entire CDS and are better rendered with polygons
draw_polygons = FALSE

# separate features by category
axistype = annottab[which(annottab[,2]=="axis"),]
# for some reason this has to be forced to numeric sometimes
window_start = as.numeric(axistype[1,3])
window_end   = as.numeric(axistype[1,4])
if (axistype[1,5]=="True") {
  bmode = TRUE
}
axis_width = window_end - window_start
offset_width = axis_width * 0.05
print(paste("# plotting axis of distance",axis_width,"for up to",num_tx,"genes"))

# mRNA or transcripts
mrnatypes = annottab[which(annottab[,2]=="mRNA"),]
# for cases where no mRNA is given, assume meant to be entire genes as blocks
if (dim(mrnatypes)[1]==0) {
	mrnatypes = annottab[which(annottab[,2]=="gene"),]
	draw_polygons = TRUE
	offset_width = axis_width * 0.01
}

is_forward = which(mrnatypes[,5]=="+")
forward_tx = mrnatypes[is_forward,]
is_reverse = which(mrnatypes[,5]=="-")
reverse_tx = mrnatypes[is_reverse,]
is_strandless = which(mrnatypes[,5]==".")
nostrand_tx = mrnatypes[is_strandless,]

exontypes = annottab[which(annottab[,2]=="exon"),]
# for cases where no mRNA is given, assume meant to be entire genes as blocks
# so rects are drawn for genes, not exons
if (dim(exontypes)[1]==0) {
	exontypes = annottab[which(annottab[,2]=="gene"),]
}

tx_names = unique(mrnatypes[,1])
print(paste("# counted", length(tx_names) ,"transcripts"))

tx_index = match(mrnatypes[,1],tx_names)
exon_index = match(exontypes[,1],tx_names)

# draw the PDF
pdf(file=outputfile, width=8, height=10)
plot(0,0, type='n', axes=FALSE, frame.plot=FALSE, ylim=c(0,num_tx), xlim=c(window_start,window_end), xlab="", ylab="")
par(mar=c(4.5,2,1,1), xpd=TRUE)
axis(1, cex.axis=1.4)
# if drawing genes alone, no exons, then draw them as polygons
if (draw_polygons==TRUE) {
	# draw forward polygons
	for (yval in tx_index[is_forward]) {
		genelen = mrnatypes[yval,4]-mrnatypes[yval,3]
		used_offset = ifelse(genelen < offset_width, genelen, offset_width)
		forward_x = c( mrnatypes[yval,4], mrnatypes[yval,4]-used_offset, mrnatypes[yval,3], mrnatypes[yval,3], mrnatypes[yval,4]-used_offset, mrnatypes[yval,4])
		if (bmode){yval = num_tx/2}
		forward_y = c( yval, yval-0.3, yval-0.3, yval+0.3, yval+0.3, yval)
		polygon( forward_x, forward_y, col="#25893a")
	}
	# draw reverse polygons
	for (yval in tx_index[is_reverse]) {
		genelen = mrnatypes[yval,4]-mrnatypes[yval,3]
		used_offset = ifelse(genelen < offset_width, genelen, offset_width)
		reverse_x = c( mrnatypes[yval,3], mrnatypes[yval,3]+used_offset, mrnatypes[yval,4], mrnatypes[yval,4], mrnatypes[yval,3]+used_offset, mrnatypes[yval,3])
		if (bmode){yval = num_tx/2}
		reverse_y = c( yval, yval-0.3, yval-0.3, yval+0.3, yval+0.3, yval)
		polygon( reverse_x, reverse_y, col="#25893a")
	}
	# draw rectangles for strandless features, usually this is an error
	for (yval in tx_index[is_strandless]) {
		rect(mrnatypes[yval,3], yval-0.3, mrnatypes[yval,4], yval+0.3, col="#25893a")
	}
# otherwise draw arrows marking direction, and each exon becomes a box
# strandless exons should be drawn but have no arrow
} else {
	# draw as arrows and boxes
	arrows(forward_tx[,3], tx_index[is_forward], forward_tx[,4]+offset_width, tx_index[is_forward], lwd=3, angle=15, length=0.1, col="#bbbbbb")
	arrows(reverse_tx[,4], tx_index[is_reverse], reverse_tx[,3]-offset_width, tx_index[is_reverse], lwd=3, angle=15, length=0.1, col="#bbbbbb")
	rect(exontypes[,3], exon_index-0.3, exontypes[,4], exon_index+0.3, col="#25893a")
}
# write the names of each transcript
if (bmode) { text(mrnatypes[,3], rep(num_tx/2, length(tx_index)), tx_names, adj=c(1,1), srt=45)
  } else {text(mrnatypes[,3], tx_index, tx_names, pos=2)}
# write the scaffold name on the margin
mtext(axistype[1,1], side=1, at=window_start-(1.5*offset_width), cex=1.8, line=-0.4)
#
dev.off()





#
