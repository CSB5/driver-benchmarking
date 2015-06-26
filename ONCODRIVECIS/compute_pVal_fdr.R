#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: directory

argv <- commandArgs(TRUE)

# Change to OncodriveCIS working directory
setwd(argv[1])

# Read in results file
amp <- read.table("OncoCNA.AMP", header=T, check.names=F)
del <- read.table("OncoCNA.DEL", header=T, check.names=F)

# Compute pVal
## amp
amp.pval <- data.frame(matrix(, nrow=dim(amp)[1], ncol=7))
colnames(amp.pval) <- c(colnames(amp), "pVal")
amp.pval[,1:6] <- amp
amp.pval[,7] <- pnorm(-(amp.pval[,6]))
## del
del.pval <- data.frame(matrix(, nrow=dim(del)[1], ncol=7))
colnames(del.pval) <- c(colnames(del), "pVal")
del.pval[,1:6] <- del
del.pval[,7] <- pnorm(del.pval[,6])

# Combine matrix
combined <- data.frame(matrix(, nrow=(dim(amp)[1]+dim(del)[1]), ncol=8))
colnames(combined) <- c(colnames(amp.pval), "FDR")
combined[,1:7] <- rbind(amp.pval, del.pval)

# Compute FDR
combined[,8] <- p.adjust(combined[,7], "fdr")

# Sort results by FDR
combined <- combined[order(combined$FDR),]

# Saving results to file
write.table(combined, file="OncoCNA.combined.txt", quote=F, sep="\t", row.names=F, col.names=T)

