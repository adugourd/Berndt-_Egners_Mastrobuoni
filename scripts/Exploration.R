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

library(readxl)
library(limma)
library(ggplot2)
library(ggfortify)

#Get the dataset
Proteom_Daten_Bl6_vs_16w_AG_Kempa <- as.data.frame(read_excel("~/Documents/Thorsten_liver_model/data/Proteom-Daten Bl6 vs 16w AG Kempa.xlsx"))
names(Proteom_Daten_Bl6_vs_16w_AG_Kempa) <- Proteom_Daten_Bl6_vs_16w_AG_Kempa[1,]

#Extract the measurement matrix, which will be used for the differential expression analysis 
batches <- Proteom_Daten_Bl6_vs_16w_AG_Kempa[-1,c(7:34)]
batches <- as.data.frame(apply(batches,2,function(x) as.numeric(x)) )
row.names(batches) <- gsub("[;].*","",Proteom_Daten_Bl6_vs_16w_AG_Kempa[-1,1])
batches[batches == 0] <- NA
batches <- log2(batches)

#This function allows you to visualise the samples, and how they compare to each others
magicPlotMaker(batches,"~/Documents/Thorsten_liver_model/results/data_visualisation/")

#Just renaming the column so it look nicer
control_names <- gsub("lfq.intensity.unknown-black6-ctrl_ctrl","Control_", names(batches))
cancer_names <- gsub("lfq.intensity.m-16w-wt","Cancer", names(batches))

#This part averages the technical replicates
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27))
{
  batches[[as.character(i)]] <- rowMeans(batches[,c(i,i+1)], na.rm = TRUE) 
}

batches <- batches[,c(29:42)]

#this part filters out the protein that are not relevant fo differential analysis (measured in only one of the two conditions)
batches <- batches[rowSums(is.na(batches)) < 13,]
names(batches)[c(1,2,3,4,5,6,7)] <- control_names[c(1,3,5,7,9,11,13)]
names(batches)[c(8,9,10,11,12,13,14)] <- cancer_names[c(15,17,19,21,23,25,27)]
batches <- batches[rowSums(is.na(batches[,c(1,2,3,4,5,6,7)])) < 7 & rowSums(is.na(batches[,c(8,9,10,11,12,13,14)])) < 7,]

#This part is to make a nice PCA plot
targets$color <- c("#DE2D26","#DE2D26","#DE2D26","#DE2D26","#DE2D26","#DE2D26","#DE2D26","#31A354","#31A354","#31A354","#31A354","#31A354","#31A354","#31A354")
pca <- prcomp(t(batches[complete.cases(batches),]))
explained <- (pca$sdev)^2 / sum(pca$sdev^2)

autoplot(pca, data = targets, colour = 'color', 
         xlab = paste('PC1 (',round(100*explained[1],digits=2),'%)',sep=''),
         ylab = paste('PC2 (',round(100*explained[2],digits=2),'%)',sep=''), size = 8, alpha = 0.5) + 
  #theme(plot.title = element_text(hjust = 0.5)) +
  xlim(c(-1,1)) +
  ylim(c(-1,1)) +
  theme_minimal() 

barplot(explained, col = "lightblue") #Plot the explained variance of the components

magicPlotMaker(batches,"~/Documents/Thorsten_liver_model/results/data_visualisation/")

###LIMMA, This is the differential analysis

#You need to make a target file that repsent the information that you have of your samples
targets <- readTargets("~/Documents/Thorsten_liver_model/src/target.txt", row.names = "Sample")
f <- factor(targets$condition, levels = c("control","cancer")) 

#create your design matrix. You can add more parameters to be taken into account by the linear model (here we only consider the two experimental condition : f)
design <- model.matrix(~0+f) 

corfit <- duplicateCorrelation(batches,design)
fit <- lmFit(batches, design,block=targets$Sample,correlation=corfit$consensus)


###Creation of the contrast matrix.
cont.matrix <- makeContrasts(tumVsCtrl = fcancer-fcontrol,
                             levels = design)

###Analysis
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#Extract the result table
tableTop <- topTable(fit2, coef = 1, adjust = "fdr", n = 1886)
results <- decideTests(fit2, p.value = 0.1, adjust.method = "fdr")

#Make nice volcano plots
volcanoplot(fit2, names = fit2$genes$ID)
volcano_nice(tableTop, 0.1, 0.5, 1, 5)

write.csv(tableTop, "~/Documents/Thorsten_liver_model/results/tableTop.csv")

tableTop <- tableTop[complete.cases(tableTop),]
###in what time point are protein changing between T and W
vennDiagram(results)

####Fold changes distributions
mat <- as.data.frame(tableTop[,1])
temp <- melt(mat)
ggplot(temp, aes(x = value, fill = variable)) + geom_density(alpha = 0.4) + ylab("Density") + xlab("Log2 Fold Change") + ggtitle("Fold changes distributions") + scale_fill_discrete(name="Conditions")

###Bubleplot
library(venneuler)
library(readr)

mouse_prot <- as.data.frame(read_delim("~/Documents/Thorsten/src/proteome_mouse", "\t", escape_double = FALSE, trim_ws = TRUE))
metabolic_prot <- as.data.frame(read_delim("~/Documents/Thorsten/src/metabolic_prot_reviewed", "\t", escape_double = FALSE, trim_ws = TRUE))
mouse_prot_reviewed <- read_delim("~/Documents/Thorsten/src/mouse_prot_reviewed", "\t", escape_double = FALSE, trim_ws = TRUE)

mouse_prot <- as.data.frame(mouse_prot$Entry)
mouse_prot$set <- rep("SwissProt+TREMBL: 81545",length(mouse_prot))
names(mouse_prot) <- c("prot","set")

mouse_prot_reviewed <- as.data.frame(mouse_prot_reviewed$Entry)
mouse_prot_reviewed$set <- rep("SwissProt: 16853",length(mouse_prot_reviewed))
names(mouse_prot_reviewed) <- c("prot","set")

metabolic_prot <- as.data.frame(metabolic_prot$Entry)
metabolic_prot$set <- rep("Metabolism: 8786", length(metabolic_prot))
names(metabolic_prot)[1] <- "prot"

detected <- as.data.frame(gsub("[;].*","",Proteom_Daten_Bl6_vs_16w_AG_Kempa$protein.ids))
detected$set <- rep("Detected: 4415", length(detected))
names(detected)[1] <- "prot"

tested <- as.data.frame(row.names(tableTop))
tested$set <- rep("Tested: 1886", length(tested))
names(tested)[1] <- "prot"

significant <- as.data.frame(row.names(tableTop[tableTop$adj.P.Val < 0.0005,]))
significant$set <- rep("FDR < 0.05: 1018", length(significant))
names(significant)[1] <- "prot"

sets <- rbind(mouse_prot,mouse_prot_reviewed)
sets <- rbind(sets,metabolic_prot)
sets <- rbind(sets,detected) 
sets <- rbind(sets,tested)
sets <- rbind(sets,significant)

plot(venneuler(sets))

sum(detected$prot %in% metabolic_prot$prot)
sum(tested$prot %in% metabolic_prot$prot)
sum(significant$prot %in% metabolic_prot$prot)
