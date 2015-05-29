#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: Root directory to DATA: eg. /mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/DATA
#argv[2]: cancer census file
#argv[3]: output directory
argv <- commandArgs(TRUE)
library(vioplot)

blca.all <- read.table(paste(argv[1], "BLCA", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
coad.all <- read.table(paste(argv[1], "COAD", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
gbm.all <- read.table(paste(argv[1], "GBM", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
lihc.all <- read.table(paste(argv[1], "LIHC", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
luad.all <- read.table(paste(argv[1], "LUAD", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
ov.all <- read.table(paste(argv[1], "OV", "mutations_per_sample-all_genes.dat", sep="/"), header=T)
prad.all <- read.table(paste(argv[1], "PRAD", "mutations_per_sample-all_genes.dat", sep="/"), header=T)

# Mutations per sample
labels = c(paste("BLCA",dim(blca.all)[1], sep="\n"), paste("COAD",dim(coad.all)[1], sep="\n"), paste("GBM",dim(gbm.all)[1], sep="\n"), paste("LIHC",dim(lihc.all)[1], sep="\n"), paste("LUAD",dim(luad.all)[1], sep="\n"), paste("OV",dim(ov.all)[1], sep="\n"), paste("PRAD",dim(prad.all)[1], sep="\n"))
#boxplot(blca.all[,2], coad.all[,2], gbm.all[,2], lihc.all[,2], luad.all[,2], ov.all[,2], prad.all[,2], main="Mutations per sample", xlab="Cancer Type", ylab="No. of mutations", labels=F, col="cornflowerblue", log="y") 
#text(x =  seq_along(labels), y = par("usr")[3] - 1, srt = 60, adj = 1.5, labels = labels, xpd = TRUE)

pdf(file=paste(argv[3], "mutations_per_sample-all_genes.pdf", sep="/"), width=8, height=8)
par(las=1,bty="l")  ## my preferred setting
par(mgp=c(3.2,1.5,0))
## set up empty plot
plot(0:1,0:1,type="n",xlim=c(0.5,7.5),ylim=range(log10(c(blca.all[,2], coad.all[,2], gbm.all[,2], lihc.all[,2], luad.all[,2], ov.all[,2], prad.all[,2], 10000))),axes=FALSE,ann=FALSE)
vioplot(log10(blca.all[,2]), log10(coad.all[,2]), log10(gbm.all[,2]), log10(lihc.all[,2]), log10(luad.all[,2]), log10(ov.all[,2]), log10(prad.all[,2]), col="cornflowerblue",add=TRUE)
axis(side=1,at=1:7,labels=labels)
axis(side=2,at=0:4,labels=10^(0:4))
title(main="Mutations per sample", xlab="Cancer Types", ylab="Mutations in log10 scale")
abline(h=log10(400), col="indianred", lwd=5, lty="dashed")
abline(h=log10(1000), col="black", lwd=5, lty="dashed")
#abline(h=log10(1000), col="indianred", lwd=5, lty="dashed")
dev.off()


# % of snv per sample
pdf(file=paste(argv[3], "percentage_snvs_per_sample-all_genes.pdf", sep="/"), width=8, height=8)
par(las=1,bty="l")  ## my preferred setting
par(mgp=c(3.2,1.5,0))
boxplot(blca.all[,3], coad.all[,3], gbm.all[,3], lihc.all[,3], luad.all[,3], ov.all[,3], prad.all[,3], main="% of SNV Mutations per sample", xlab="Cancer Type", ylab="% of SNV mutations", names=labels, col="cornflowerblue")
dev.off();

pdf(file=paste(argv[3], "percentage_snvs_per_sample-all_genes-vioplot.pdf", sep="/"), width=8, height=8)
par(las=1,bty="l")  ## my preferred setting
par(mgp=c(3.2,1.5,0))
plot(0:1,0:1,type="n",xlim=c(0.5,7.5),ylim=c(0,1),axes=FALSE,ann=FALSE)
vioplot(blca.all[,3], coad.all[,3], gbm.all[,3], lihc.all[,3], luad.all[,3], ov.all[,3], prad.all[,3], col="cornflowerblue", add=TRUE)
axis(side=1,at=1:7,labels=labels)
axis(side=2)
title(main="% of SNV Mutations per sample", xlab="Cancer Types", ylab="% of SNV mutations")		
dev.off();


# read cancer census file
CC <- read.table(argv[2], header=T, colClasses=c("character", rep("NULL",15)), sep="\t")

geneFreq <- function(inFile){
	data <- read.table(inFile, header=T, colClasses=c("character", "double", "NULL"), sep="\t")
	rownames(data) <- data$GeneName
	common <- intersect(CC[,1], rownames(data))
	plot(density(data$Frequency), col="forestgreen")
	lines(density(data[common,2]), col="indianred")
}