# display KEGG BlastKOALA pie chart with and without the unannotated proteins
# to give an accurate representation of the uncertainty when reporting common categories
# created by WRF 2025-02-07

category_data_file = "~/git/genomeGTFtools/test_data/avas_bin1_kegg_category_counts_w_unknowns.txt"
category_data = read.table(category_data_file, header=TRUE, sep="\t")

pdf(file="~/git/genomeGTFtools/test_data/avas_bin1_kegg_category_counts_w_unknowns.pdf", width=8, height=6, title="KEGG annotation with and without unknowns")
par(mfrow=c(1,2), mar=c(7,1,3,1), xpd=TRUE)
pie(category_data$count[-1] , radius = 1.0, clockwise = TRUE, col=paste0("#",category_data$color[-1]) , labels=NA , main="KEGG annotation\nBlastKOALA output" )
legend(-1,-1, legend=category_data$category[2:9], bty = 'n' , 
       pch=22, pt.bg = paste0("#",category_data$color[2:9]) , pt.cex = 2 )
pie(category_data$count , radius = 1.0, clockwise = TRUE, col=paste0("#",category_data$color) , labels=NA, main="Actual fractions\nincluding 43% unknowns" )
legend(-1,-1, legend=category_data$category[c(10:15,19,20)], bty = 'n', 
       pch=22, pt.bg = paste0("#",category_data$color[c(10:15,19,20)]) , pt.cex = 2 )
dev.off()



#