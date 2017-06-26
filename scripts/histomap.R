##This script allows you diplay results of an enrichment analysis of a specific contrast. The idea is to selct a few pathway of relevance and display the foldchanges of their most significant genes along with a clear belonging to specific pathways.
##The required dataframes are standard output from limma and PIANO enrichment analysis

library(gridExtra)
library(grid)
library(RColorBrewer)
library(readr)
library(ggplot2)
library(pheatmap)
library(reshape)

##import consensus table
consensus_mm_no_reac <- read_csv("~/Documents/Thorsten_liver_model/results/consensus_mm_no_reac.csv")

##import the results of differential expression
tableTop <- read.csv("~/Documents/Thorsten/tableTop.csv")
tableTop <- tableTop[complete.cases(tableTop),]

##import the gene to term table used
gene_to_term <- read_csv("~/Documents/Thorsten_liver_model/src/gene_to_term.csv")

##if a mapping table was used to match identifiers when the enrichment analysis was performed, you will need it here too
mapping_uniprot_to_name <- read_delim("~/Documents/Thorsten_liver_model/src/mapping_uniprot_to_name.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
mapping_uniprot_to_name$name <- toupper(mapping_uniprot_to_name$name)
names(mapping_uniprot_to_name)[2] <- "gene"

gene_to_term <- merge(gene_to_term, mapping_uniprot_to_name)
gene_to_term <- gene_to_term[,c(3,2)]

names(gene_to_term)[1] <- "X" #to match the identifier column of the tabletop

tableTop <- merge(tableTop, gene_to_term) #now each omic feature will be associated with the corresponding terms

##Here you should provide the relevant terms for your analysis
sub_terms <- c("KEGG_MM_METABOLIC_PATHWAYS","KEGG_MM_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION","KEGG_MM_FATTY_ACID_METABOLISM","KEGG_MM_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM","KEGG_MM_CITRATE_CYCLE","KEGG_MM_OXIDATIVE_PHOSPHORYLATION")

##Subset the tabletop for the relevant terms
sub_df <- tableTop[tableTop$term %in% sub_terms,]

##This is optional. The p value threshold depends on how many genes you want to show. More than 200 genes will become really messy.
sub_df <- sub_df[sub_df$adj.P.Val <= 0.1,]
names(sub_df)[1] <- "X"

##now we create the heatmap that will tell which gene belong to which term
heatmap <- cast(sub_df, term~X, value = "logFC")
row.names(heatmap) <- heatmap$term
heatmap <- heatmap[,-1]
heatmap[length(heatmap[,1])+1,] <- colMeans(heatmap,na.rm = TRUE)
heatmap[is.na(heatmap)] <- 1
heatmap <- heatmap[,order(heatmap[length(heatmap[,1]),])]

##Then we create the barplot to show the foldchange of the selected genes
to_bar <- heatmap[length(heatmap[,1]),]
to_bar <- t(to_bar)
to_bar <- cbind(to_bar, rep(1,length(heatmap[1,])))
to_bar[to_bar[,1] < 0 ,2] <- 0
to_bar <- as.data.frame(to_bar)
to_bar$X <- names(heatmap)
to_bar$X <- factor(to_bar$X, levels = to_bar$X)
to_bar$V1 <- abs(to_bar$V1)
to_bar$couleur <- rep("#31a354",length(heatmap[1,]))
to_bar[to_bar$V2 == 0,4] <- "#de2d26"

##Clean the heatmap a bit
heatmap <- heatmap[-length(heatmap[,1]),]
heatmap[heatmap != 1] <- 0
row.names(heatmap) <- gsub("KEGG_MM_","",row.names(heatmap))
row.names(heatmap) <- gsub("_"," ", row.names(heatmap))

##Create the plot object for the heatmap
a <- pheatmap(heatmap, cluster_rows = FALSE, cluster_cols = FALSE, colorRampPalette(rev(brewer.pal(n = 5, name ="Blues")))(2), fontsize_col = 5)

##Create the plot object for the barplot
b <- ggplot(to_bar, aes(x = X, y = V1)) + geom_bar(stat = "identity", width = 0.84, fill = to_bar$couleur) + theme_void() +
  theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = NULL) +
  scale_y_continuous(expand = c(0,0))

##Plot the barplot on top of the heatmap
lay <- rbind(c(1,NA),c(3,3))
grid.arrange(b, a$gtable)

##Then at this point you will need to export the resulting plot. I recommend using Inkscape to adjust the plots togethere