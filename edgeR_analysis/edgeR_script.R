# First part of the edgeR analysis; lines with # signs are just control commands to check the anlysis

library(edgeR)

pval = 0.01 # set pvalue threshold (only relevant for plotting)
pval_adjust = "BH"  # set method for p-value adjustment for multiple testing (= FDR calculation)
cutoff = c(-1,1)  # set log2 fold change line in plots (e.g. "c(-1,1)" means lines will be drawn at 2-fold up and down)
res_list <- list()
res_list2 <- list()
DGElists <- list()


ID_Names <- read.table('ID_Names_Asp_Kle.txt', row.names = 1) # Import table for addition of gene names
names(ID_Names) <- c('Name')

targets_Asp <- readTargets(file = "Targets_Asp.txt") # this imports the data in the Targets_Asp.txt file = experiment setup
targets_Kle <- readTargets(file = "Targets_Kle.txt") # this imports the data in the Targets_Kle.txt file = experiment setup

targets = list(targets_Asp, targets_Kle)

for (i in targets) {

  d <- readDGE(i, header = F) # this reads in the count data
  #d$samples
  #head(d$counts)
  #summary(d$counts)
  #dim(d)

  # Filter low expression tags (cpm<1)
  keep <- rowSums(cpm(d)> 1) >= 3 # this creates an R object with the genes/transcripts that have more than 1 counts per million - change the 1 if you want to modify the filter; the 3 at the end is the group size - change this if you have more or less replicates than 3
  d <- d[keep,] # this filters out all genes/transcripts with less than 1 counts per million
  #dim(d)
  d$samples$lib.size <- colSums(d$counts) # this calculates the library size again after the filtering step
  #d$samples

  # Normalization (TMM)
  d <- calcNormFactors(d) # this calculates the normalization factors used during the diffrential expression analysis
  #d$samples

  # Estimating the dispersions and plot them
  d <- estimateCommonDisp(d, verbose=F) # this calculates the common dispersion over all genes and samples - is used during the diff. expr. anaylsis
  d <- estimateTagwiseDisp(d) # this calculates the dispersion for each gene/transcript - is used during the diff. expr. anaylsis
  DGElists <- c(DGElists, list(d))
  
  pdf(paste(i[1,2], "_", i[length(i[,2]),2], "_disp.pdf", sep = ""))  # change the name of the .pdf file; this opens a pdf file to store the next plot
  plotBCV(d)
  dev.off()  # closes and saves the pdf file with the plot
  
  # Differential expression, lines with # signs are just control commands to check the analysis
  
  et <- exactTest(d, pair= i[c(1,length(i[,2])),2]) # change the names in quotes to your group names; this command performs the differential expression analysis
  res_list2 <- c(res_list2, list(et))
  
  #topTags(et, n=20)
  #detags <- rownames(topTags(et, n=20))
  #cpm(d)[detags,]
  result <- summary(de <- decideTestsDGE(et, p=pval, adjust=pval_adjust)) # this command gives the number of differentially expressed genes; add "lfc = 1" for logFC -1/1 as additional cutoff
  #result
  detags <- rownames(d)[as.logical(de)] # this command saves the row names of differentially expressed genes as "detags"
  pdf(paste(i[1,2], "_", i[length(i[,2]),2], ".pdf", sep = ""))  # change the name of the .pdf file; this opens a pdf file to store the next plot
  plotSmear(et, de.tags=detags, xlab = "Average log2 counts per million", ylab = "log2 fold change") # this plots the analysis results; use "panel.first = NULL" to remove grid
  abline(h = cutoff, col = "blue")
  dev.off() # closes and saves the pdf file with the plot
  
  
  # Export results
  
  
  
  res <-as.data.frame(topTags(et, n=Inf)) # this saves the fold change table in "tt"
  res$FC <- 2^res$logFC   # include a fold change column (non-log)
  
  resNames <- merge(res, ID_Names, by = 0, all.x = T)
  names(resNames) <- c('ID', 'log2_FoldChange', 'log2_counts_per_million', 'pvalue', 'adj_pvalue', 'FoldChange', 'Name')
  
  write.table(resNames[,c(1,7,2,6,5,3)], file=paste(i[1,2], "_", i[length(i[,2]),2], ".txt", sep = ""), sep="\t", quote=FALSE, row.names = F) # this exports the fold change table; change the file name of the .txt file as you want
  
  cpmill = as.data.frame(d$pseudo.counts) # this saves the counts per million table in a new R object called "pcounts_1h" - change this name as you want
  names(cpmill) = i$description # this changes the column names in the previously created "pcounts_1" table;
  write.table(cpmill, file=paste(i[1,2], "_", i[length(i[,2]),2], "_cpm.txt", sep = ""), sep="\t", quote=FALSE) # this exports the counts per million table; change the file name as you want
  
  res_list <- c(res_list, list(res, cpmill))
  
}

### PCA analysis
### Aspergillus

Asp_transp <- t(res_list[[2]])
mds <- cmdscale(dist(Asp_transp), k=3, eig=T)
eig_pc <- mds$eig *100 /sum(mds$eig)
png(file="PCA_PropExplainedVariance_Asp.png")
barplot(eig_pc, las=1, xlab="Dimensions", ylab="Proportion of explained variance (%)", y.axis=NULL, col="darkgrey")
dev.off()

mds <- cmdscale(dist(Asp_transp))
rownames(mds) <- c("Asp1", "Asp2", "Asp3", "Asp+Kle1", "Asp+Kle2", "Asp+Kle3")
plot(mds[,1], -mds[,2], type="p", xlab="Dimension 1", ylab="Dimension 2", main="", pch = 19, cex = 1.5)
text(mds[,1], -mds[,2], rownames(mds), cex=1, adj = c(0.6,1), pch=2)

### Klabsiella

Kle_transp <- t(res_list[[4]])
mds <- cmdscale(dist(Asp_transp), k=3, eig=T)
eig_pc <- mds$eig *100 /sum(mds$eig)
png(file="PCA_PropExplainedVariance_Kle.png")
barplot(eig_pc, las=1, xlab="Dimensions",  ylab="Proportion of explained variance (%)", y.axis=NULL, col="darkgrey")
dev.off()

mds <- cmdscale(dist(Kle_transp))
rownames(mds) <- c("Kle1", "Kle2", "Kle3", "Kle+Asp1", "Kle+Asp2", "Kle+Asp3")
plot(mds[,1], -mds[,2], type="p", xlab="Dimension 1", ylab="Dimension 2", main="", pch = 19, cex = 1.5)
text(mds[,1], -mds[,2], rownames(mds), cex=1, adj = c(0.6,1), pch=2)
