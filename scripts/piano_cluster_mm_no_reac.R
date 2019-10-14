#Copyright (C) 2019  Aurelien Dugourd

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(piano)
library(reshape)
library(readr)

tableTop <- read.csv("/home/ad234505/Thorsten_liver_model/results/tableTop.csv")
tableTop <- tableTop[complete.cases(tableTop),]

gene_to_term <- as.data.frame(read_csv("/home/ad234505/Thorsten_liver_model/src/gene_to_term.csv"))

mapping_unip_to_geneName <- read.delim("/home/ad234505/Thorsten_liver_model/src/mapping_uniprot_to_name.txt")
names(mapping_unip_to_geneName) <- c("unip","gene")
mapping_unip_to_geneName$gene <- toupper(mapping_unip_to_geneName$gene)

gene_to_term <- merge(gene_to_term, mapping_unip_to_geneName, all = FALSE)

gene_to_term_no_reac <- gene_to_term[grep("REACTOME",gene_to_term[,"term"], invert = TRUE),]
geneSet <- loadGSC(gene_to_term_no_reac[,c(3,2)])

myFC <- tableTop$logFC
names(myFC) <- tableTop$X

myPval <- tableTop$adj.P.Val
names(myPval) <- tableTop$X

myTval <- tableTop$t
names(myTval) <- tableTop$X

table(gene_to_term_no_reac$term)
hist(table(gene_to_term_no_reac$term), breaks= 1000)

###Run the GSA
gsaRes1 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "mean")
gsaRes2 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "median")
gsaRes3 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "sum")
gsaRes4 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "maxmean")
gsaRes5 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "reporter")
gsaRes6 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "tailStrength")
gsaRes7 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "wilcoxon")
gsaRes8 <- runGSA(myTval, gsc=geneSet, adjMethod = "fdr", geneSetStat = "page")
gsaRes9 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "fisher")
gsaRes10 <- runGSA(geneLevelStats = myPval, directions = myFC, gsc=geneSet, adjMethod = "fdr", geneSetStat = "stouffer")

resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes8, gsaRes9, gsaRes10)
names(resList) <- c("mean","median","sum","maxmean","reporter","tailStrength","wilcoxon","page", "fisher", "stouffer")

ch <- consensusHeatmap(resList,cutoff=50,method="median", ncharLabel = 50, cellnote = "medianPvalue", cex = 0.2, plot = FALSE) ##The results are strange

consensus_mm_no_reac <- ch$pMat

write.csv(consensus_mm_no_reac,"/home/ad234505/Thorsten_liver_model/results/consensus_mm_no_reac.csv")
