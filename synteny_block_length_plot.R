# generate quick histogram of synteny block sizes from microsynteny.py
#
# block expected in format of 11 column tabular:
# tethya01_7684	Contig13491	blk-1	twi_ss.23503	143941	146353	-	Aqu2.34641_001	403531	405934	+
#
# last modified 2018-07-04

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "acropora_vs_styllophora_microsynteny.tab"

blockdata = read.table(inputfile, sep="\t")
counttable = table(blockdata[,3])

longestblock = max(counttable)
lencap = min(longestblock, 20) # take whichever is smaller, longest block or 20

counttable[counttable>lencap] = lencap

outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)
pdf(outputfile, width=6, height=5)
par(mar=c(4.2,5,2,2))

barcolor = c("#76ee94") # green
#barcolor = c("#fc8d62") # red/orange
#barcolor = c("#80b1d3") # blue

counthist = hist(counttable, breaks=seq(2.5,longestblock+0.5,1), axes=FALSE, plot=TRUE, xlim=c(3,lencap), main=inputfile, xlab="Number of genes in block", ylab="Number of blocks", cex.lab=1.5, col=barcolor )
mostcommon = max(counthist$counts)

axis(1,at=seq(3,lencap,3), labels=seq(3,lencap,3), cex.axis=1.2 )
axis(2, cex.axis=1.3)
text(counthist$breaks+0.5, counthist$counts+(0.025*mostcommon), counthist$counts)

totalblk = length(counttable)
totallab = paste(totalblk, "total blocks", sep=" ")

totalgenes = length(unique(blockdata[,4]))
genelab = paste(totalgenes, "total genes in blocks", sep=" ")
text(lencap, mostcommon*0.65, totallab, cex=1.5, pos=2)
text(lencap, mostcommon*0.50, genelab, cex=1.5, pos=2)

dev.off()

#