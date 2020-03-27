#!/usr/bin/env Rscript

library(Biostrings)

## keep command line argument parsing simple for now - argparse would be better
args <- commandArgs(TRUE)
scaffname <- args[1]
scaffname_filtered <- args[2]
min_len <- as.numeric(args[3])
min_cov <- as.numeric(args[4])

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
contigs<-contigs[cov > min_cov]  ## fix for spades
writeXStringSet(contigs, scaffname_filtered)
