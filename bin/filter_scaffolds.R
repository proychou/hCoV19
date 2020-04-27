#!/usr/bin/env Rscript

library(Biostrings)

## keep command line argument parsing simple for now - argparse would be better
args <- commandArgs(TRUE)
sample <- args[1]
scaffname <- args[2]
scaffname_filtered <- args[3]
filtered_report <- args[4]
min_len <- as.numeric(args[5])
min_cov <- as.numeric(args[6])

## min_width <- 200
## ## min_cov <- 10
## ## temporary fix for spades
## min_cov <- 5

## import scaffolds and filter by length (>200) and coverage (>10x)
contigs<-readDNAStringSet(scaffname, format='fasta')
contigs<-contigs[width(contigs) > min_len]
cov<-unlist(lapply(names(contigs), function(x){
  as.numeric(strsplit(x,'_cov_')[[1]][2])
}))
report<-data.frame(
  'sample'=sample,
  'mean_coverage'=mean(cov),
  'mean_length'=mean(width(contigs)))
contigs<-contigs[cov > min_cov]  ## fix for spades
report$'scaffolds'<-length(contigs)
write.csv(report, filtered_report, row.names=FALSE)
writeXStringSet(contigs, scaffname_filtered)
