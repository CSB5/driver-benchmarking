#!/mnt/software/bin/Rscript-3.1.0 --slave

# how to run this script:
# ./ActiveDriver.R -a '/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/DATA/GBM/annovar.hg19_multianno.txt' -f '/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/SOFTWARE_TESTBED/ActiveDriver/ActiveDriver.results'

# Check if ActiveDriver is installed and install it otherwise.
if(!require("ActiveDriver")){
  install.packages("http://individual.utoronto.ca/reimand/ActiveDriver/ActiveDriver_0.0.10.tar.gz", repo=NULL, type="source")
}

suppressPackageStartupMessages(library(ActiveDriver))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-a", "--annovar_file"), action="store", default=NA, type='character',
              help="directory containing .results files"),
  make_option(c("-f", "--fout"), action="store", default=NA, type='character',
              help=".results output file path")
)
opt = parse_args(OptionParser(option_list=option_list))

# directory with phosphosite/sequence and sequence disorder files
dname = "/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/SOFTWARE_TESTBED/ActiveDriver/ActiveDriverData/Reimand_NatSciRep2013/"

# Map Gene name to Longest Isoform
finlongestTranscript <-"/mnt/projects/bertrandd/oncoimpact/MUTATION_BENCHMARK/SOFTWARE_TESTBED/ActiveDriver/ActiveDriverData/Reimand_NatSciRep2013/gene_symbol_to_refseq.tab"

longestTranscript <- read.table(finlongestTranscript, stringsAsFactors = FALSE)
colnames(longestTranscript)[1] <- 'Gene.refGene'

dAnnovar <- read.delim(opt$a, stringsAsFactors = FALSE, header = FALSE)
dAnnovar <- dAnnovar[ ,c(39,41,42,ncol(dAnnovar))]
colnames(dAnnovar) <- c(dAnnovar[1,1],dAnnovar[1,2],dAnnovar[1,3],'sample_id')
dAnnovar <- dAnnovar[-1,]

dAnnovar <- merge(dAnnovar, longestTranscript)
dAnnovar <- dAnnovar[grepl('nonsynonymous SNV', dAnnovar$ExonicFunc.refGene), ]
rownames(dAnnovar) <- NULL


# the lines which contain more than one "synonymous SNV" entry in the ExonicFunc.knownGene field
# separated with ; contain repeated info in the  ExonicFunc.knownGene and AAChange.knownGene fields
# I could only find 1 corresponding entry for the sample_id in the .maf file, so I am not sure why the
# repeated info is in the file. The script below just takes the first instance of the repeated data.

muts <- mdply(.data = 1:nrow(dAnnovar), .fun = function(x){
  AAchange <- unlist(strsplit(dAnnovar$AAChange.refGene[x],';'))[1]
  temp <- regmatches(AAchange, regexpr(paste0(dAnnovar$rseq[x],':.*?p.([A-Z])([0-9]+)([A-Z])'),AAchange))
  if(length(temp) == 0){#case that the longest Transcript is not found
    return(NULL)
  }
  else{
    temp <- gsub('.*p.([A-Z])([0-9]+)([A-Z])',"\\1\\2\\3",temp,perl=TRUE)
    temp.cha <- unlist(str_extract_all(temp, "[A-Z]"))
    position <- as.numeric(str_extract(temp, "[0-9]+"))
    return(data.frame(gene = dAnnovar$Gene.refGene[x], cancer_type = 'GBM', sample_id = dAnnovar$sample_id[x], position, wt_residue = temp.cha[1], mut_residue = temp.cha[2]))
  }
})

# load required datasets
sites = read.table(paste0(dname,"all_phosphosites_for_TCGA_pancancer.tab"))
seqs = read_fasta(paste0(dname,"sequences_for_TCGA_pancancer.fa"))
disorder = read_fasta(paste0(dname,"sequence_disorder_for_TCGA_pancancer.fa"))

# run ActiveDriver
psnv_info = ActiveDriver(seqs, disorder, muts, sites)

# save genes ranked by p-values in a .result file
dout <- psnv_info$all_gene_based_fdr
dout <- dout[order(dout$fdr),]
dout <-dout[(dout$fdr <= 0.01),] # threshold used in Reimand 2013, doi:10.1038/srep02651
dout <- data.frame(Gene_name = dout$gene, Sample = rep('ALL',nrow(dout)), Rank = 1:nrow(dout), pvalue = dout$p , FDR = dout$fdr)
write.table(dout, opt$f, sep='\t', quote=F, row.names=F)
