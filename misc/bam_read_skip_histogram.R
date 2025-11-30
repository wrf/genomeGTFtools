#!/usr/bin/env Rscript
# bam_read_skip_histogram.R
# make histogram of read skip S from CIGAR string, to detect trans spliced leader sequence
# created by WRF 2025-11-28
# last modified 2025-11-28

# map long reads with mapper, such as minimap2
# use with get_read_skip_from_bam.py, first or last S in CIGAR string is extracted for each mapped read
# ~/samtools-1.9/samtools view UCSC_Hcal_v1_B1_LR.sorted.bam | get_read_skip_from_bam.py - > UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab
# output from get_read_skip_from_bam.py is used as the input here

# Rscript bam_read_skip_histogram.R UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab

args = commandArgs(trailingOnly=TRUE)

skipdata_file = args[1]
#skipdata_file = "~/genomes/hormiphora_californensis/UCSC_Hcal_v1_B1_LR.sorted.read_skip.tab"
# in case of .gz input, remove the .gz and rename otherwise
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz$","",skipdata_file,perl=TRUE),perl=TRUE)
if (skipdata_file==outputfile) { stop("cannot parse input file to generate output file name, add a unique 3-letter suffix") }

skipdata = read.table(skipdata_file, sep="\t")

forward_skips = table(skipdata[,2])
#forward_skips[1:100]
reverse_skips = table(skipdata[,3])

local_max = forward_skips[which(forward_skips[]==max(forward_skips[(-1:-5)]))] * 1.2
above_1pct = which(forward_skips>(0.02*local_max))

pdf(file=outputfile, height=6, width=7, title="Read skip from mapping")
plot(0,0,type='n', xlim=c(0,100), ylim=c(0,local_max), 
     ylab="Count", xlab="S length (bp)", main=basename(skipdata_file) )
lines( as.integer(unlist(dimnames(forward_skips))), forward_skips , lwd=3 , col="#12348b77" )
lines( as.integer(unlist(dimnames(reverse_skips))), reverse_skips , lwd=3 , col="#7c123477" )
text( as.integer(unlist(dimnames(forward_skips[above_1pct]))), forward_skips[above_1pct], 
      as.integer(unlist(dimnames(forward_skips[above_1pct]))), pos=3)
text( 100, local_max*0.95, pos=2,
      paste( "30-50bp\n", 
             round(100*( sum(forward_skips[30:50]) + sum(reverse_skips[30:50]) ) / 
                     sum(c(forward_skips,reverse_skips) /2 ),digits=2 ), "%" ) )
dev.off()

# must divide by 2 since both forward and reverse are possible
# Hcal
# ( sum(forward_skips[35:48]) + sum(reverse_skips[35:48]) ) / sum(c(forward_skips,reverse_skips)/2)
#[1] 0.5636931



#