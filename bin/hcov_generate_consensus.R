#!/usr/bin/env Rscript
## RSV : This script imports bam files and makes a consensus sequence
## Pavitra Roychoudhury
## Adapted from hsv_generate_consensus.R on 6-Mar-19

library(Rsamtools);
library(GenomicAlignments);
library(ShortRead);
library(Biostrings);

##Takes in a bam file, produces consensus sequence
generate_consensus<-function(bamfname){
  require(Rsamtools)
  require(GenomicAlignments)
  require(Biostrings)
  require(parallel)
  ncores<-detectCores()

  ##for testing this function--comment out or remove
  ## bamfname<-'./testing/ABI-HHV6A_S385_L001_A.sorted.bam'

  if(!is.na(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){

    ##Index bam if required
    if(!file.exists(paste(bamfname,'.bai',sep=''))){
      baifname<-indexBam(bamfname);
    }else{
      baifname<-paste(bamfname,'.bai',sep='');
    }

    ##Import bam file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);
    ## summary(gal);

    ##Remove any contigs with mapq <2 -- this leads to a loss of a lot of the DR seqs even though there are reads there
    ## gal<-gal[mcols(gal)$mapq>2];

    ##First lay reads on reference space--this doesn't include insertions
    qseq_on_ref<-sequenceLayer(mcols(gal)$seq,cigar(gal),from="query",to="reference");

    ##Make a consensus matrix and get a consensus sequence from the aligned scaffolds
    ## cm<-consensusMatrix(qseq_on_ref,as.prob=T,shift=start(gal)-1,width=seqlengths(gal))[c('A','C','G','T','N','-'),];
    ## cm['N',colSums(cm)==0]<-1;

    ##Edit to include a coverage threshold
    cm<-consensusMatrix(qseq_on_ref,as.prob=F,shift=start(gal)-1,width=seqlengths(gal))[c('A','C','G','T','N','-'),];
    poor_cov<-which(colSums(cm)<10);
    cm<-apply(cm,2,function(x)x/sum(x));
    cm[,poor_cov]<-0;
    cm['N',poor_cov]<-1;

    tmp_str<-strsplit(consensusString(cm,ambiguityMap='?',threshold=0.5),'')[[1]];
    ambig_sites<-which(tmp_str=='?');
    ambig_bases<-unlist(lapply(ambig_sites,function(i){mixedbase<-paste(names(cm[,i])[cm[,i]>0],collapse='');
      if(mixedbase%in%IUPAC_CODE_MAP) return(names(IUPAC_CODE_MAP)[IUPAC_CODE_MAP==mixedbase])
      else return('N')}));
    tmp_str[ambig_sites]<-ambig_bases
    con_seq<-DNAStringSet(paste0(tmp_str,collapse=''));
    names(con_seq)<-sub('.bam','_consensus',basename(bamfname));
    rm(tmp_str);

    ##Remove gaps and leading and trailing Ns to get final sequence
    con_seq_trimmed<-DNAStringSet(gsub("N*N$",'',gsub("^N*",'',as.character(con_seq))));
    con_seq_final<-DNAStringSet(gsub('-','',as.character(con_seq_trimmed)));
    names(con_seq_final)<-sub('.bam','_consensus',basename(bamfname));

    ##Delete bai file
    file.remove(baifname);

    return(con_seq_final);
  }else{
    return(NA)
  }
}

#Return the number of mapped reads in a bam file
n_mapped_reads<-function(bamfname){
  require(Rsamtools)
  indexBam(bamfname)
  if(file.exists(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    return(idxstatsBam(bamfname)$mapped)
  }else{
    return(NA)
  }
}

args <- commandArgs(TRUE)

sample <- args[1]

## consensus is generated from remapped_bam
remapped_bam <- args[2]

## stats are calculated from mappedtoref_bam
mappedtoref_bam <- args[3]
refname <- args[4]

## outputs
final_cons <- args[5]
stats_csv <- args[6]

con_seq <- generate_consensus(remapped_bam)
writeXStringSet(con_seq, file=final_cons, format='fasta')

num_Ns <- sum(letterFrequency(con_seq,c('N','+')))
width <- width(con_seq)

mapping_stats <- data.frame(
    sample=sample,
    ref=refname,
    remapped_bam=NA,
    mappedtoref_bam=NA,
    # mapped_reads_ref=unlist(lapply(mappedtoref_bam, n_mapped_reads)),
    mapped_reads_assemblyref=unlist(lapply(remapped_bam, n_mapped_reads)),
    perc_Ns=100 * num_Ns / width,
    num_Ns=num_Ns,
    width=width,
    stringsAsFactors=FALSE
)

write.csv(mapping_stats, file=stats_csv, row.names=FALSE)


