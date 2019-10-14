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

library(BioNet)
library(igraph)
library(XML)
library(readr)

tableTop <- as.data.frame(read_csv("~/Documents/Thorsten_liver_model/results/tableTop.csv"))
tableTop <- as.data.frame(tableTop[complete.cases(tableTop),])
omniPKN <- read_delim("~/Documents/Thorsten/omniPKN_16_03_17.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
STRING_PPI_800 <- read_delim("~/Documents/Thorsten/src/STRING_PPI_800.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
STRING_PPI_800 <- as.data.frame(STRING_PPI_800)

STRING_uniprot_to_Gene <- read_delim("~/Documents/Thorsten/src/STRING_uniprot_to_Gene.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ttop_uniprot_to_Gene <- read_delim("~/Documents/Thorsten_liver_model/src/ttop_unip_to_name.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(ttop_uniprot_to_Gene) <- c("unip","gene")

names(STRING_uniprot_to_Gene)[1] <- "Accession"
STRING_PPI_800 <- merge(STRING_PPI_800,STRING_uniprot_to_Gene)
STRING_PPI_800 <- STRING_PPI_800[,c(7,2,3,4,5,6)]
names(STRING_uniprot_to_Gene) <- c("Accession2","Gene2")
STRING_PPI_800 <- merge(STRING_PPI_800,STRING_uniprot_to_Gene)
STRING_PPI_800 <- STRING_PPI_800[,c(2,7,3,4,5,6)]

names(ttop_uniprot_to_Gene)[1] <- "X1"
tableTop <- merge(tableTop,ttop_uniprot_to_Gene)
tableTop <- tableTop[,c(8,2,3,4,5,6,7)]


row.names(tableTop) <- tableTop$gene
pval <- as.vector(tableTop$P.Val)
names(pval) <- row.names(tableTop)
fdr <- as.vector(tableTop$adj.P.Val)
names(fdr) <- row.names(tableTop)
logFC <- tableTop$logFC
names(logFC) <- row.names(tableTop)

#interactome <- graph_from_data_frame(omniPKN, directed = TRUE, vertices = NULL)
interactome <- graph_from_data_frame(STRING_PPI_800, directed = TRUE, vertices = NULL)

subnet <- subNetwork(row.names(tableTop), interactome,neighbors="first")
subnet <- rmSelfLoops(subnet)
subnet<-largestComp(subnet)
subnet_nodes <- (V(subnet)$name)
stats_genes <- intersect(subnet_nodes,row.names(tableTop))

fb <- fitBumModel(pval)
fb
scores <- scoreNodes(subnet, fb, fdr = 0.0005)
scores <- round(scores, 4)

writeHeinzEdges(network=subnet, file="~/Documents/Thorsten_liver_model/src/Network_edges.txt", use.score=FALSE)
writeHeinzNodes(network=subnet, file="~/Documents/Thorsten_liver_model/src/Network_nodes.txt", node.scores=scores)

run_heinz_cmd <- paste('heinz-mc -e ', "~/Documents/Thorsten_liver_model/src/Network_edges.txt", ' -n ', "~/Documents/Thorsten_liver_model/src/Network_nodes.txt", ' -o ', "~/Documents/Thorsten_liver_model/results/Heinz_solution.txt",sep='')
system(run_heinz_cmd)

module <- readHeinzGraph(node.file="~/Documents/Thorsten_liver_model/results/Heinz_solution.txt", network=subnet)
nodeDataDefaults(module, attr='score') <- ''
nodeDataDefaults(module, attr='FC') <- ''
nodeDataDefaults(module, attr='FDR') <- ''
nodeData(module, n=nodes(module), attr='score') <- scores[nodes(module)]
nodeData(module, n=nodes(module), attr='FC') <- logFC[nodes(module)]
nodeData(module, n=nodes(module), attr='FDR') <- fdr[nodes(module)]
edgeDataDefaults(module, attr = 'sign') <- ''
temp <- STRING_PPI_800[STRING_PPI_800[,1] %in% nodes(module) & STRING_PPI_800[,2] %in% nodes(module),]
temp <- temp[temp[,1] != temp[,2],]
#temp[temp[,4] == "activation",4] <- 1
#[temp[,4] == "inhibition",4] <- -1
edgeData(module, temp[,1], temp[,2], attr = 'sign') <- as.character(temp[,4])
#par(mfrow = c(1,1))
#plotModule(module, scores=scores, diff.expr = logFC)

module_igraph <- igraph.from.graphNEL(module)

write_graph(module_igraph, "~/Documents/Thorsten_liver_model/results/bionet_tumor_over_control.graphml", format = "graphml")
saveNetwork(module,file="~/Documents/Thorsten_liver_model/results/network_tumor_over_control", type='XGMML')

write(as.character(nodes(module)),"~/Documents/Thorsten_liver_model/src/module_nodes.txt")
write(as.character(tableTop$gene),"~/Documents/Thorsten_liver_model/src/background.txt")
