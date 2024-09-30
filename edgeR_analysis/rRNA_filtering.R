## Filter count files for rRNA genes

rRNAs <- read.table(file.choose())  ## select a file containing one column with all the rRNA gene IDs without a header
rRNAs$ID <- rep("rRNA", length(rRNAs[,1]))

files <- list.files(pattern = "crop.txt")
ext1 <- "filtAf.txt"
ext2 <- "filtKp.txt"
## files_out <- lapply(files, function(x) paste(x, ext, sep = "."))  ## not necessary, can be done in loop

for (i in files) {
  s <- read.table(i)
  sm <- merge(s, rRNAs, all.x = T)
  smf <- sm[is.na(sm$ID),]
  smf <- smf[,c(1,2)]
  smfA1 <- smf[grep("^Afu", smf$V1),]   ### Filter out all fungal genes
  smfA2 <- smf[grep("^U1", smf$V1),]
  smfA <- rbind(smfA1, smfA2)
  smfK <- smf[grep("^gene", smf$V1),]   ### Filter out all bacterial genes
  write.table(smfA, file = paste(i, ext1, sep = "."), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(smfK, file = paste(i, ext2, sep = "."), col.names = F, row.names = F, quote = F, sep = "\t")
}
