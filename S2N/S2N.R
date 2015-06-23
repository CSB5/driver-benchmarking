#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: expression matrix
#argv[2]: cnv matrix
#argv[3]: outDir

library(CNAmet)
argv <- commandArgs(TRUE)

expr <- read.table(argv[1], header = TRUE)
cna <- read.table(argv[2], header = TRUE)
exprMatrix   <- as.matrix(expr, rownames.force=TRUE)  
cnaMatrix   <- as.matrix(cna, rownames.force=TRUE)  
results <- CNAmet(exprMatrix = exprMatrix, cghMatrix = cnaMatrix,
                  methylMatrix = NULL, perms = 1000, 
                  na.limit = 0.1, gainData = TRUE, favorSynergetic = TRUE,
                  strictChecks = FALSE, strictLim = 0.05)
                  
write.table(results, file=paste(argv[3], "S2N.result", sep="/"), sep="\t", row.names=T, col.names=T, quote=FALSE)
