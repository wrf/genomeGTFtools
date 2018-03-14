#!/usr/bin/env Rscript

# convert protein GTF format to rectangle blocks for domains

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/sycon_ciliatum/28757.cdd.gtf"
outputfile = gsub("([\\w/]+).g.f","\\1.pdf",inputfile,perl=TRUE)

gtftab = read.table(inputfile,header=FALSE,sep="\t")

# separate features by category
proteintypes = gtftab[which(gtftab[,3]=="protein"),]
signaltypes = gtftab[which(gtftab[,3]=="signal_peptide"),]
PFAMtypes = gtftab[which(gtftab[,3]=="PFAM"),]

# residue types must be generated separately with alignmentpos2gff.py
residuetypes = gtftab[which(gtftab[,3]=="modified residue"),]
motiftypes = gtftab[which(gtftab[,3]=="motif"),]


protnames = unique(gtftab[,1])
numprots = length(protnames)

# get domain position information
protindex = match(PFAMtypes[,1],protnames)
yupper = protindex*10-3
ylower = protindex*10-7
domstart = PFAMtypes[,4]
domend = PFAMtypes[,5]

# domain codes for legend and color coding
domaincodes = gsub("ID=([\\w\\d]+)\\..*","\\1",PFAMtypes[,9],perl=TRUE)
domnames = unique(domaincodes)
numdomains = length(domnames)
domainindex = match(domaincodes,domnames)
colvec = rainbow(numdomains, alpha=0.8)
domainids = gsub("ID=([\\w\\d]+)\\.([\\w\\d-]+)\\..+","\\2",PFAMtypes[,9],perl=TRUE)

# OUTPUT
pdf(outputfile, width=10, height=numprots+2)
par(mar=c(2.5,2.5,1,2) )

maxend = max(gtftab[,5])
protlenrange = pretty(c(0,maxend))

plot(0,0,type='n',xlim=c(0,max(protlenrange)),ylim=c(0,numprots*11+10),frame.plot=FALSE,xlab="",ylab="", axes=FALSE)
axis(1)

# make black lines for length of the protein
# this works even if no protein types are given, instead just does not draw lines
ymiddles = match(proteintypes[,1],protnames)*10-5
linexcoords = c(rbind(proteintypes[,4],proteintypes[,5],NA,NA))
lineycoords = c(rbind(ymiddles,ymiddles,NA,NA))
lines(linexcoords,lineycoords, lwd=2)

# make signal peptides as black boxes
signalindex = match(signaltypes[,1],protnames)
sigupper = signalindex*10-3
siglower = signalindex*10-7
sigstart = signaltypes[,4]
sigend = signaltypes[,5]
rect(sigstart, siglower, sigend, sigupper, col="#000000")

# draw domains as rectangles
#axis(1,at=protlenrange)
rect(domstart, ylower, domend, yupper, col=colvec[domainindex])
#text(domstart, protindex*10-1, domainids, pos=4, offset=0)
text(c(0),match(protnames,protnames)*10-1, protnames, pos=4, offset=0)
legend(0,numprots*11+10,legend=domainids[match(domnames,domaincodes)],col=colvec,pch=15,x.intersp=0.5,ncol=6)

# draw modified residues as triangles
residueindex = match(residuetypes[,1],protnames)
points(residuetypes[,4],residueindex*10-4, pch=6, cex=1.33)

# draw selected motifs as letters themselves
motifindex = match(motiftypes[,1],protnames)
motifnames = gsub("Note=(\\w+)","\\2",motiftypes[,9], perl=TRUE)
points(motiftypes[,4],motifindex*10-7, pch=motifnames)

dev.off()