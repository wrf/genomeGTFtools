# calculation of probabilities of matching genes between scaffolds, as having derived from a random distribution
# by WRF 2020-06-15
# last modified 2022-12-27

# analysis is modeled from this explanation from Srivastava 2008
# "Scaffold_8 has genes belonging to 174 of the 5,487 ancestral gene clusters, while human chromosome segment 10.4
# has genes belonging to 165 of them. 54 ancestral gene clusters have genes from both Scaffold_8 and human 
# chromosome segment 10.4. Fisher's exact test rejects the null model that the gene complements of Trichoplax 
# scaffold_8 and human chromosome segment 10.4 are drawn independently from the set of ancestral genes, and 
# that all ancestral genes have an equal probability of having a descendant on a given segment"

#s2008_matrix = matrix(data=c(54, 174, 165, 5487), nrow=2)
#s2008_ft = fisher.test(test_matrix)
#s2008_ft$p.value
# [1] 4.75784e-30

# define function to calculate p values for each scaffold pair across both genomes
calculate_ftest_pvalues <- function(freq_matrix) {
	# start from empty vector
	p_values = c()
	nscaf_g1 = dim(freq_matrix)[1]
	nscaf_g2 = dim(freq_matrix)[2]
	match_total = sum(freq_matrix)
	# loop over each matrix position, calculate Fishers test
	for ( i in 1:nscaf_g1) {
		for (j in 1:nscaf_g2) {
			h1h1 = freq_matrix[i,j]
			h1h0 = sum(freq_matrix[i])
			h0h1 = sum(freq_matrix[j])
			h0h0 = match_total - h1h0 - h0h1 + h1h1
			test_matrix = matrix( data=c(h1h1,h1h0,h0h1,h0h0), nrow=2)
			ft = fisher.test(test_matrix)
			p_values = c(p_values, ft$p.value)
		}
	}
	return(p_values)
}


args = commandArgs(trailingOnly=TRUE)

# read data file from scaffold_synteny.py
all2Dfile = args[1]

# read optional species names
genome1_arg = args[2]
genome2_arg = args[3]

all2Ddata = read.table(all2Dfile, sep="\t", stringsAsFactors=FALSE)
categories = all2Ddata[,1]
is_scaf1 = which(categories=="s1")
scafdata1 = all2Ddata[is_scaf1,]
longestscaf1 = max(scafdata1[,6])
is_longscafs1 = which(as.numeric(scafdata1[,5]) > 0.0001)
is_scaf2 = which(categories=="s2")
scafdata2 = all2Ddata[is_scaf2,]
longestscaf2 = max(scafdata2[,6])
is_longscafs2 =  which(as.numeric(scafdata2[,5]) > 0.0001)
is_points = which(categories=="g")
pointsdata = all2Ddata[is_points,]
if (sum(is.na(pointsdata[,3])) > 0 | sum(is.na(pointsdata[,5])) > 0){
  print("WARNING: NAs DETECTED IN POINTS, OUTPUT MAY NOT WORK")}
scaffold_data = data.frame( sc1 = pointsdata[,3], sc2 = pointsdata[,5])

match_freq_table = table(scaffold_data)
# must na.omit in case some scaffold has no matches
scaf1_names = scafdata1[,2]
sc1_reorder_index = na.omit( match(scaf1_names,row.names(match_freq_table))[is_longscafs1] )
scaf2_names = scafdata2[,2]
sc2_reorder_index = na.omit( match(scaf2_names,colnames(match_freq_table))[is_longscafs2] )
# maybe need to switch indexing
#match(row.names(match_freq_table), scaf1_names[is_longscafs1])

reordered_table = match_freq_table[sc1_reorder_index,sc2_reorder_index]

match_freq_matrix = as.matrix(reordered_table)

nscaf_g1 = dim(match_freq_matrix)[1]
nscaf_g2 = dim(match_freq_matrix)[2]
match_total = sum(match_freq_matrix)

ntests = nscaf_g1 * nscaf_g2

ftest_p_values = calculate_ftest_pvalues(match_freq_matrix)
ftest_p_value_adj = ftest_p_values * ntests
#range(ftest_p_values)
ftest_p_value_log = floor(-log10(ftest_p_value_adj))
ftest_hist = hist(ftest_p_value_log, breaks=20, plot=FALSE) # ylim=c(0,500),

ftest_p_value_log[ftest_p_value_log<1] = 1
ftest_p_value_log[ftest_p_value_log>100] = 100


# final matrix for graph
ftest_p_value_mat = matrix( data=ftest_p_value_log, nrow=nscaf_g1, byrow=TRUE )

# raw p value matrix for csv
raw_p_value_mat = matrix( data=ftest_p_values, nrow=nscaf_g1, byrow=TRUE )
colnames(raw_p_value_mat) = scaf2_names[sc2_reorder_index]
# write csv
output_csv = gsub("([\\w/]+)\\....$","\\1.pvalue.csv",gsub(".gz$","",all2Dfile,perl=TRUE),perl=TRUE)
print( paste("writing csv to",output_csv) )
write.csv(raw_p_value_mat, file=output_csv, row.names=scaf1_names[sc1_reorder_index] )
#


# count blocks that are significant against the randomized dataset
has_randomized = FALSE
if (has_randomized == TRUE) {
	# Nb with lower value than the global min of the randomized set
	length( ftest_p_values[ftest_p_values<(min(ftest_p_values)/ntests)] )
	# NBcS, Bonferroni corrected scaffolds, each relative to the randomized
	length( ftest_p_values[(ftest_p_values / ftest_p_values) < (0.01/ntests)] )
}



# make color slices for the scalebar
scalebar_range = c(1,rev(pretty(ftest_p_value_log))[1])
scale_slice_width = 0.4/scalebar_range[2]
scale_x_starts = seq(0.5, 0.9-scale_slice_width, scale_slice_width)
scale_x_endss = seq(0.5+scale_slice_width, 0.9, scale_slice_width)

colorscheme = colorRampPalette(c("#f1eef6","#66a4c2","#4d004b"))(scalebar_range[2])


outputfile = gsub("([\\w/]+)\\....$","\\1.pvalue.pdf",gsub(".gz$","",all2Dfile,perl=TRUE),perl=TRUE)
print( paste("creating pdf",outputfile) )

main_label = paste("Fisher's exact test P values\n", genome1_arg, "vs", genome2_arg)
# create pdf as full page on a4
pdf(file=outputfile, width=8, height=11, paper="a4")
par( mar=c(5.5,4.5,4,1) )
# draw heatmap
image(z=ftest_p_value_mat, col=colorscheme, axes=FALSE, zlim=scalebar_range, main=main_label)
# label scaffolds
mtext( scaf1_names, side=1, at=seq(0,1,1/(length(scaf1_names)-1)), las=2, cex=0.95)
mtext( scaf2_names,side=2,at=seq(0,1,1/(length(scaf2_names)-1)), las=1, cex=0.75)
# make scalebar
rect( scale_x_starts, rep(0.9,scalebar_range[2]+1), scale_x_endss , rep(0.99,scalebar_range[2]+1), col=colorscheme, border=FALSE)
text(0.5,0.85, "p=0.1", cex=2)
text(0.9,0.85, paste0("p<1e-",scalebar_range[2]),  cex=2)
# notes at bottom of page
range_text = paste("Min corrected P-value:", formatC(min(ftest_p_value_adj),format="e", digits=4 )  )
mtext(range_text,side=1, line=3, at=0.1)
dev.off()





#
