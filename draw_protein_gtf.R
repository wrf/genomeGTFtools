#!/usr/bin/env Rscript
#
# convert protein GTF format to rectangle blocks for domains
# last modified 2022-04-28

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/git/genomeGTFtools/test_data/nidogen_full_prots.clan.gff"
outputfile = gsub("([\\w/]+).g..$","\\1.pdf",inputfile,perl=TRUE)

print(paste("Reading domains from", inputfile, ", writing to", outputfile))
gtftab = read.table(inputfile,header=FALSE,sep="\t")


# reverse order to match fasta file
protnames = unique(gtftab[,1])
numprots = length(protnames)

if (numprots < 13) {
  # for 12 or less, squeeze onto one page
  yscale_max = numprots*11+10
  countgroups = list( seq(1,numprots,1) )
} else {
  # implicitly this will split to 2 or more pages
  yscale_max = 110
  countgroups = split(1:numprots, ceiling(seq_along(1:numprots)/10))
  print(paste("PDF should contain", length(countgroups), "pages"))
}


# separate features by category
proteintypes = gtftab[which(gtftab[,3]=="protein"),]
signaltypes = gtftab[which(gtftab[,3]=="signal_peptide"),]
PFAMtypes = gtftab[which(gtftab[,3]=="PFAM"),]

# residue types must be generated separately with alignmentpos2gff.py
# this feature was used to indicate triangles of Prolines in HIF and related factors
# for the paper Mills, Francis, et al 2018 eLife
residuetypes = gtftab[which(gtftab[,3]=="modified residue"),]
motiftypes = gtftab[which(gtftab[,3]=="motif"),]

# domain codes for legend and color coding
domaincodes = gsub("ID=([\\w\\d]+)\\..*","\\1",PFAMtypes[,9],perl=TRUE)
domnames = unique(domaincodes)
numdomains = length(domnames)
domainindex = match(domaincodes,domnames)
colvec = rainbow(numdomains, s=0.82, alpha=0.8)
domainids = gsub("ID=([\\w\\d]+)\\.([\\w\\d-]+)\\..+","\\2",PFAMtypes[,9],perl=TRUE)

# use same scalebar and maximum for all pages
maxend = max(gtftab[,5])
protlenrange = pretty(c(0,maxend))

# set offsets for block width
offs_u = 8
offs_l = 12
offs_m = 10

# OUTPUT
pdf(outputfile, width=8, height=11, title = basename(inputfile) )
par(mar=c(2.5,2.5,1,2) )
for (i in 1:length(countgroups)) {
  # makes list of numbers for each group, should return something like:
  # 1 2 3 4 5 6 7 8 9 10  or  11 12 13 14 15 16
  prots_per_page = as.integer(unlist(countgroups[i]))
  # the names for that group, in reverse order, so 1st prot in the FASTA file appears at the top
  # restrict matches so that only ones on this page are found by match()
  # just a list of protein IDs, like:
  # "NID2_MOUSE" "NID2_HUMAN" "NID1_HUMAN" "NID1_MOUSE"
  protnames_per_page = rev(as.character(protnames[prots_per_page]))
  
  # get domain position information
  protindex = match(PFAMtypes[,1],protnames_per_page)
  yupper = protindex*10 - offs_u
  ylower = protindex*10 - offs_l
  domstart = PFAMtypes[,4]
  domend = PFAMtypes[,5]
  
  plot(0,0,type='n', xlim=c(0,max(protlenrange)), ylim=c(0,yscale_max),
       frame.plot=FALSE, xlab="", ylab="", axes=FALSE)
  axis(1)
  
  # make black lines for length of the protein
  # this works even if no protein types are given, instead just does not draw lines
  ymiddles = match(proteintypes[,1],protnames_per_page)*10 - offs_m
  linexcoords = c(rbind(proteintypes[,4],proteintypes[,5],NA,NA))
  lineycoords = c(rbind(ymiddles,ymiddles,NA,NA))
  lines(linexcoords,lineycoords, lwd=2)
  
  # make signal peptides as black boxes
  signalindex = match(signaltypes[,1],protnames_per_page)
  sigupper = signalindex*10 - offs_u
  siglower = signalindex*10 - offs_l
  sigstart = signaltypes[,4]
  sigend = signaltypes[,5]
  rect(sigstart, siglower, sigend, sigupper, col="#000000")
  
  # draw domains as rectangles
  #axis(1,at=protlenrange)
  rect(domstart, ylower, domend, yupper, col=colvec[domainindex])
  #text(domstart, protindex*10-1, domainids, pos=4, offset=0)
  text(c(0),match(protnames_per_page,protnames_per_page)*10 - (offs_u-2), 
       protnames_per_page, pos=4, offset=0)
  legend(0, yscale_max, legend=domainids[match(domnames,domaincodes)],
         col=colvec, pch=15, x.intersp=0.5, ncol=4, bty='o')
  
  # draw modified residues as downward pointing triangles
  residueindex = match(residuetypes[,1],protnames_per_page)
  points(residuetypes[,4],residueindex*10-4, pch=6, cex=1.33)
  
  # draw selected motifs as letters themselves
  motifindex = match(motiftypes[,1],protnames_per_page)
  motifnames = gsub("Note=(\\w+)","\\2",motiftypes[,9], perl=TRUE)
  points(motiftypes[,4],motifindex*10-7, pch=motifnames)
}
dev.off()
#



#