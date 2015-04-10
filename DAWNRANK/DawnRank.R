#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: expression matrix
#argv[2]: mutation matrix
#argv[3]: adjacency network matrix

library(DawnRank);

argv <- commandArgs(TRUE)

expressionMatrix <- read.table(argv[1], header=T);
mutationMatrix <- read.table(argv[2], header=T);
adjMatrix <- read.table(argv[3], header=T);


#To run
#Recommandd parameters
dawnRankScore<-DawnRank(adjMatrix=adjMatrix, expressionMatrix=expressionMatrix, mutationMatrix=mutationMatrix, mu = 20, maxit = 100, epsilon = 1e-04,goldStandard = NULL);
#Test data set command
#dawnRankScore<-DawnRank((adjMatrix=adjMatrix, expressionMatrix=expressionMatrix, mutationMatrix=mutationMatrix, mu=3,goldStandard=NULL)

#To get the output
print("## Moving on to aggregate step.##\n")
aggregateDawnRankScore<-condorcetRanking(scoreMatrix=dawnRankScore[[2]], mutationMatrix=mutationMatrix);
write.table(aggregateDawnRankScore[[2]], file="driver_list.dat", sep="\t", row.names=T, col.names=T, quote=FALSE);



