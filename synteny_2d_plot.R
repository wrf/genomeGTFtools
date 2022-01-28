# synteny_2d_plot.R
# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# created by WRF 2019-04-01
# last modified 2019-04-01

args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py
all2Dfile = args[1]

# read optional species names
genome1_arg = args[2]
genome2_arg = args[3]

if (!is.na(genome1_arg)) {
xlab = paste( gsub("-", " ", genome1_arg),"(total Mb)")
} else {
xlab = "Genome 1 (total Mb)"
}

if (!is.na(genome2_arg)) {
ylab = paste( gsub("-", " ", genome2_arg),"(total Mb)")
} else {
ylab = "Genome 2 (total Mb)"
}

# read all data in a single file
all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
#head(all2Ddata)

# breakdown categories of either two genomes s1 and s2, or gene hits g
categories = all2Ddata[,1]

is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
#scafdata1

longestscaf1 = max(scafdata1[,6])
#longestscaf1
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0009)
#is_longscafs1

longscafs1 = c(0, scafdata1[,6][is_longscafs1] )

is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0009)
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )

is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
#head(pointsdata)

# determine which genome is longer, for correct orientation on paper
if (longestscaf1 > longestscaf2) {
	genome_x = pointsdata[,7]
	genome_y = pointsdata[,6]
	xmax = tail( pretty(longscafs2), n=1)
	ymax = tail( pretty(longscafs1), n=1)
	longscafs_x = longscafs2
	longscafs_y = longscafs1
	nscafs_x = length(longscafs2)
	nscafs_y = length(longscafs1)
} else {
	genome_x = pointsdata[,6]
	genome_y = pointsdata[,7]
	xmax = tail( pretty(longscafs1), n=1)
	ymax = tail( pretty(longscafs2), n=1)
	longscafs_x = longscafs1
	longscafs_y = longscafs2
	nscafs_x = length(longscafs1)
	nscafs_y = length(longscafs2)
}

xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)

# larger bitscores make larger points
pointsize = log10(as.numeric(pointsdata[,8])) / 4

# make PDF
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",all2Dfile,perl=TRUE)
pdf(file=outputfile, width=8, height=11)
par( mar=c(4.5,4.5,1,1) )

# dotcolor = "#18935188" # green
# dotcolor = "#c51b8a88" # pink
# dotcolor = "#1c909988" # teal
# dotcolor = "#2071d388" # blue
# dotcolor = "#88419d88" # purple

plot(genome_x, genome_y, pch=16, col="#18935188", cex=pointsize, main=all2Dfile, xlab=xlab, ylab=ylab, axes=FALSE, cex.lab=1.4)

tickpoints = pretty(c(0,xmax_mb))
axis(1, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3)
barpos_x = rep(c( ymax*-0.01, ymax*-0.02),round(nscafs_x)/2)
segments( longscafs_x[1:(nscafs_x-1)], barpos_x[1:(nscafs_x-1)], longscafs_x[2:nscafs_x], barpos_x[0:(nscafs_x-1)], lwd=3)
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.1, col="#777777")


tickpoints = pretty(c(0,ymax_mb))
axis(2, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3, line=0.5)
barpos_y = rep(c( xmax*-0.01, xmax*-0.02),round(nscafs_y)/2)
segments( barpos_y[1:(nscafs_y-1)], longscafs_y[1:(nscafs_y-1)], barpos_y[0:(nscafs_y-1)], longscafs_y[2:nscafs_y], lwd=3)
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.1, col="#777777")

# display numbers beside y-axis scaffold segments
#textpos_y = rep(c( xmax*-0.02, xmax*-0.01 ),round(24)/2)
#textmidbar = as.numeric(scafdata1[1:24,6]) - as.numeric(scafdata1[1:24,4])/2
#text(textpos_y, textmidbar, 1:24, cex=0.5)

dev.off()




#