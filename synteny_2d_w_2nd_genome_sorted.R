# synteny_2d_plot.R
# make dot plot of synteny between two genomes, based on unidirectional blast hits (i.e. not reciprocal)
# created by WRF 2019-04-01
# last modified 2022-12-14

args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py AND MUST HAVE USED OPTION --local-positions
all2Dfile = args[1]
# in case of .gz input, remove the .gz and rename otherwise
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz","",all2Dfile),perl=TRUE)
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
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0009)
longscafs1 = c(0, scafdata1[,6][is_longscafs1] )
longscafs1_names = scafdata1[is_longscafs1,2]

is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0009)
scafdata2_long = scafdata2[is_longscafs2,]
longscafs2 = c(0, scafdata2[,6][is_longscafs2] )
longscafs2_names = scafdata2_long[,2]

is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
is_both_longscaf = !is.na(match(pointsdata[,3], longscafs1_names)) & !is.na(match(pointsdata[,5], longscafs2_names))
pointsdata_long = pointsdata[is_both_longscaf,]

scaffold_match_data = data.frame( sc1 = pointsdata_long[,3], sc2 = pointsdata_long[,5] )
match_freq_table = table(scaffold_match_data)

scaffold_best_match = apply(match_freq_table, 2, which.max)
#scaffold_best_match
# moves alphabetic to length order
sc2_alpha_to_by_length = match( names(scaffold_best_match) , scafdata2[is_longscafs2,2] )

# moves length to alphabetic
sc2_by_length_to_alpha = match( scafdata2[is_longscafs2,2] , names(scaffold_best_match) )

sc2_sorted_index = sort(scaffold_best_match[sc2_by_length_to_alpha], index.return = TRUE)
scafdata2_reorder = scafdata2_long[sc2_sorted_index$ix,]
sc2_reorder_index = match( scafdata2_reorder[,2], names(sc2_sorted_index$x) )

# recount total length
longscafs1 = cumsum(c(0,as.numeric(scafdata1[is_longscafs1,4])))
longscafs2 = cumsum(c(0,as.numeric(scafdata2_reorder[,4])))
longestscaf1 = max(longscafs1)
longestscaf2 = max(longscafs2)
# reassign dots by scaffold
adj_points1_positions = pointsdata_long[,6] + longscafs1[match( pointsdata_long[,3], longscafs1_names ) ]
adj_points2_positions = pointsdata_long[,7] + longscafs2[match( pointsdata_long[,5], scafdata2_reorder$V2 ) ]


# determine which genome is longer, for correct orientation on paper
# if genome 1 is longer, then swap the positions, axes, and labels
if (longestscaf1 > longestscaf2) {
  genome_x = adj_points2_positions
  genome_y = adj_points1_positions
  xmax = tail( pretty(longscafs2), n=1)
  ymax = tail( pretty(longscafs1), n=1)
  longscafs_x = longscafs2
  longscafs_y = longscafs1
  nscafs_x = length(longscafs2)
  nscafs_y = length(longscafs1)
  xlab = genome2_lab
  ylab = genome1_lab
} else {
  genome_x = adj_points1_positions
  genome_y = adj_points2_positions
  xmax = tail( pretty(longscafs1), n=1)
  ymax = tail( pretty(longscafs2), n=1)
  longscafs_x = longscafs1
  longscafs_y = longscafs2
  nscafs_x = length(longscafs1)
  nscafs_y = length(longscafs2)
  xlab = genome1_lab
  ylab = genome2_lab
}

xmax_mb = round(xmax / 1000000)
ymax_mb = round(ymax / 1000000)

# larger bitscores make larger points
pointsize = log10(as.numeric(pointsdata[,8])) / 4

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
    dotcolor = "#93183688" # default red
  }
} else {
  # value is not given at all in command line
  dotcolor = "#93183688" # default red
}

# make PDF
pdf(file=outputfile, width=8, height=11) # a4 size
#pdf(file=outputfile, width=16, height=23) # a2 size
#pdf(file=outputfile, width=32, height=45) # a0 size
par( mar=c(4.5,4.5,1,1) )

plot(genome_x, genome_y, pch=16, 
     xlim = c(0,60000000), ylim = c(0,80000000),
     col=dotcolor, cex=0.5, cex.lab=1.4, 
     main="", xlab=xlab, ylab=ylab, 
     axes=FALSE )

tickpoints = pretty(c(0,xmax_mb))
axis(1, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3)
barpos_x = rep(c( ymax*-0.01, ymax*-0.02),round(nscafs_x)/2)
segments( longscafs_x[1:(nscafs_x-1)], barpos_x[1:(nscafs_x-1)], longscafs_x[2:nscafs_x], barpos_x[0:(nscafs_x-1)], lwd=3, col = "#00000088")
segments(longscafs_x, 0, longscafs_x, longscafs_y[nscafs_y], lwd=0.1, col="#88888888")

tickpoints = pretty(c(0,ymax_mb))
axis(2, at=tickpoints*1000000, labels=tickpoints, cex.axis=1.3, line=0.5)
barpos_y = rep(c( xmax*-0.01, xmax*-0.02),round(nscafs_y)/2)
segments( barpos_y[1:(nscafs_y-1)], longscafs_y[1:(nscafs_y-1)], barpos_y[0:(nscafs_y-1)], longscafs_y[2:nscafs_y], lwd=3, col = "#00000088")
segments( 0, longscafs_y, longscafs_x[nscafs_x], longscafs_y, lwd=0.8, col="#00000088")

# make shortened names
# short_names1 = gsub("scaffold_", "", longscafs1_names)
# short_names2 = gsub("scaffold_", "", longscafs2_names)

# display numbers beside axis scaffold segments
# this controls bars along the left side, so actually Y axis
max_x_n = 25 # max numbers to display
textpos_x = rep(c( ymax*-0.02, ymax*-0.01 ),round(max_x_n)/2)
textmidbar = as.numeric(scafdata1[1:max_x_n,6]) - as.numeric(scafdata1[1:max_x_n,4])/2
#text(textpos_x, textmidbar, short_names1[1:max_x_n], cex=0.5)

# this controls bars along the bottom side, so actually X axis
max_y_n = length(longscafs_x)-1
scaf2_over_1M = as.numeric(scafdata2_reorder[,4]) > 1000000
textpos_y = rep(c( xmax*-0.02, xmax*-0.01 ),round(max_y_n)/2)
textmidbar = ( as.numeric(longscafs_x[1:(max_y_n)]) + as.numeric(longscafs_x[2:(max_y_n+1)]) ) / 2
#text(textmidbar[scaf2_over_1M], textpos_y[scaf2_over_1M], short_names2[scaf2_over_1M], cex=0.5)

dev.off()

#
