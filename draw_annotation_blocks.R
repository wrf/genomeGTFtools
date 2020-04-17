#!/usr/bin/env Rscript
# convert tabular annotation format to rectangle blocks
# or polygon blocks for single CDS genes, like bacteria
# last modified 2020-04-17

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/mnemiopsis_leidyi/ml0011_63k-88k_annot_short.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

annottab = read.table(inputfile,header=FALSE,sep="\t",stringsAsFactors=FALSE)

num_tx = 30
seqcount = args[2]
if (!is.na(seqcount)) {
	num_tx = seqcount
}

draw_polygons = FALSE

# separate features by category
axistype = annottab[which(annottab[,2]=="axis"),]
axis_width = axistype[1,4] - axistype[1,3]
offset_width = axis_width * 0.05
print(paste("# plotting axis of distance",axis_width,"for up to",num_tx,"genes"))

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


pdf(file=outputfile, width=8, height=10)
plot(0,0, type='n', axes=FALSE, frame.plot=FALSE, ylim=c(0,num_tx), xlim=c(axistype[1,3],axistype[1,4]), xlab="", ylab="")
par(mar=c(4.5,2,1,1))
axis(1, cex.axis=1.4)
if (draw_polygons==TRUE) {
	for (yval in tx_index[is_forward]) {
		genelen = mrnatypes[yval,4]-mrnatypes[yval,3]
		if (genelen < offset_width) {
			used_offset = genelen
		} else {
			used_offset = offset_width
		}
		forward_x = c( mrnatypes[yval,4], mrnatypes[yval,4]-used_offset, mrnatypes[yval,3], mrnatypes[yval,3], mrnatypes[yval,4]-used_offset, mrnatypes[yval,4])
		forward_y = c( yval, yval-0.3, yval-0.3, yval+0.3, yval+0.3, yval)
		polygon( forward_x, forward_y, col="#25893a")
	}
	for (yval in tx_index[is_reverse]) {
		genelen = mrnatypes[yval,4]-mrnatypes[yval,3]
		if (genelen < offset_width) {
			used_offset = genelen
		} else {
			used_offset = offset_width
		}
		reverse_x = c( mrnatypes[yval,3], mrnatypes[yval,3]+used_offset, mrnatypes[yval,4], mrnatypes[yval,4], mrnatypes[yval,3]+used_offset, mrnatypes[yval,3])
		reverse_y = c( yval, yval-0.3, yval-0.3, yval+0.3, yval+0.3, yval)
		polygon( reverse_x, reverse_y, col="#25893a")
	}
} else {
	# draw as arrows and boxes
	arrows(forward_tx[,3], tx_index[is_forward], forward_tx[,4]+offset_width, tx_index[is_forward], lwd=3, angle=15, length=0.1)
	arrows(reverse_tx[,4], tx_index[is_reverse], reverse_tx[,3]-offset_width, tx_index[is_reverse], lwd=3, angle=15, length=0.1)
	rect(exontypes[,3], exon_index-0.3, exontypes[,4], exon_index+0.3, col="#25893a")
}
text(mrnatypes[,3], tx_index, tx_names, pos=2)
mtext(axistype[1,1], side=1, at=axistype[1,3]-(1.5*offset_width), cex=1.8, line=-0.4)
dev.off()





#