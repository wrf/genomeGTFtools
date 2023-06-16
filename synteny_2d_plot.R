# synteny_2d_plot.R
# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# created by WRF 2019-04-01
# last modified 2022-12-13

# command line arguments:
# Rscript synteny_2d_plot.R tabular_file.tab Label-1 Label-2 HUE
# example:
# Rscript synteny_2d_plot.R synteny.tab Genus-species1 Genus-species2 128
# Species-1 and Species-2 are text strings to be used as labels for the plot
# HUE is a value between 1-255 to change the hue
args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py
all2Dfile = args[1]
# in case of .gz input, remove the .gz and rename otherwise
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz$","",all2Dfile,perl=TRUE),perl=TRUE)
if (all2Dfile==outputfile) { stop("cannot parse input file to generate output file name, add a unique 3-letter suffix") }

# read optional species names
# should match between query and subject in scaffold_synteny.py, as -f and -F
genome1_arg = args[2]
genome2_arg = args[3]

#genome1_arg = NA
#genome2_arg = NA

if (!is.na(genome1_arg)) {
genome1_lab = paste( gsub("-", " ", genome1_arg),"(total Mb)")
} else {
genome1_lab = "Genome 1 (total Mb)"
}

if (!is.na(genome2_arg)) {
genome2_lab = paste( gsub("-", " ", genome2_arg),"(total Mb)")
} else {
genome2_lab = "Genome 2 (total Mb)"
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
longscafs1_names = scafdata1[is_longscafs1,2]

is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0009)
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
longscafs2_names = scafdata2[is_longscafs2,2]

is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]

# determine which genome is longer, for correct orientation on paper
# if genome 1 is longer, then swap the positions, axes, and labels
if (longestscaf1 > longestscaf2) {
	genome_x = pointsdata[,7]
	genome_y = pointsdata[,6]
	xmax = tail( pretty(longscafs2), n=1)
	ymax = tail( pretty(longscafs1), n=1)
	longscafs_x = longscafs2
	longscafs_y = longscafs1
	nscafs_x = length(longscafs2)
	nscafs_y = length(longscafs1)
	xlab = genome2_lab
	ylab = genome1_lab
	short_names1 = gsub("scaffolds_", "", longscafs2_names)
	short_names2 = gsub("scaffolds_", "", longscafs1_names)
} else {
	genome_x = pointsdata[,6]
	genome_y = pointsdata[,7]
	xmax = tail( pretty(longscafs1), n=1)
	ymax = tail( pretty(longscafs2), n=1)
	longscafs_x = longscafs1
	longscafs_y = longscafs2
	nscafs_x = length(longscafs1)
	nscafs_y = length(longscafs2)
	xlab = genome1_lab
	ylab = genome2_lab
	short_names1 = gsub("Scaffolds_", "", longscafs1_names)
	short_names2 = gsub("scaffolds_", "", longscafs2_names)
}

xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)

# larger bitscores make larger points
pointsize = log10(as.numeric(pointsdata[,8])) / 4

# dotcolor = "#18935188" # green
# dotcolor = "#c51b8a88" # pink
# dotcolor = "#1c909988" # teal
# dotcolor = "#2071d388" # blue
# dotcolor = "#88419d88" # purple

# approximate s/v, with varied h
# this uses an optional 4th argument, of value from 1 to 256
# or use random color with either "r" or "R"
# otherwise uses default green
dot_color_set = rainbow(256, s=0.8, v=0.7, alpha = 0.53)
if ( !is.na(args[4]) ) {
  # use H value specified by user, from 1 to 256
  # can be converted to integer, and is between 1 and 256
  if ( !is.na(as.integer(args[4])) & as.integer(args[4]) > 0 & as.integer(args[4]) <= 256) {
    dotcolor = dot_color_set[ as.integer(args[4]) ]
  } else if ( args[4]=="r" | args[4]=="R") {
    # make random color, upon user request
    dotcolor = dot_color_set[sample(1:length(dot_color_set), 1)]
  } else {
    # meaning user entered something outside 1-256, or not letter r
    dotcolor = "#18935188" # default green
  }
} else {
  # value is not given at all in command line
  dotcolor = "#18935188" # default green
}


# make PDF
pdf(file=outputfile, width=8, height=11) # a4 size
#pdf(file=outputfile, width=16, height=23) # a2 size
#pdf(file=outputfile, width=32, height=45) # a0 size

par( mar=c(4.5,4.5,1,1) )

plot(genome_x, genome_y, pch=16, col=dotcolor, cex=pointsize, 
     main=all2Dfile, xlab=xlab, ylab=ylab, 
     axes=FALSE, cex.lab=1.4)

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

# display numbers beside axis scaffold segments
# this controls bars along the left side, so actually Y axis
# max_x_n = 25 # max numbers to display
# textpos_x = rep(c( ymax*-0.02, ymax*-0.01 ),round(max_x_n)/2)
# textmidbar = as.numeric(scafdata1[1:max_x_n,6]) - as.numeric(scafdata1[1:max_x_n,4])/2
# #text(textpos_x, textmidbar, 1:25, cex=0.5)
# text(textpos_x, textmidbar, short_names1[1:max_x_n], cex=0.5)
# 
# # this controls bars along the bottom side, so actually X axis
# max_y_n = 30
# textpos_y = rep(c( xmax*-0.02, xmax*-0.01 ),round(max_y_n)/2)
# textmidbar = as.numeric(scafdata2[1:max_y_n,6]) - as.numeric(scafdata2[1:max_y_n,4])/2
# #text(textmidbar, textpos_y, 1:25, cex=0.5)
# text(textmidbar, textpos_y, short_names2[1:max_y_n], cex=0.5)

dev.off()




#