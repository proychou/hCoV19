#!/usr/bin/env Rscript
## RSV : This script imports bam files and makes a consensus sequence
## Pavitra Roychoudhury
## Adapted from hsv_generate_consensus.R on 6-Mar-19

## Built to be called from hhv6_wgs_pipeline.sh with input arguments specifying input filename
## Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

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


clean_consensus_hcov<-function(sampname,remapped_bamfname,mappedtoref_bamfname,ref){
  require(Rsamtools);
  require(GenomicAlignments);
  require(Biostrings);
  mapping_stats<-data.frame(ref=ref,
                            remapped_bam=remapped_bamfname,
                            mappedtoref_bam=mappedtoref_bamfname,
                            mapped_reads_ref=0,mapped_reads_assemblyref=0,perc_Ns=0,num_Ns=0,width=0,
                            stringsAsFactors=F);

  #Import mapped reads + assembly and generate consensus
  con_seq<-generate_consensus(mapping_stats$remapped_bam);
  if(!dir.exists('./consensus_seqs')) dir.create('./consensus_seqs');
  writeXStringSet(con_seq,file=paste('./consensus_seqs/',sampname,'.fasta',sep=''),format='fasta');

  #Compute #mapped reads and %Ns
  mapping_stats$mapped_reads_ref<-unlist(lapply(mapping_stats$mappedtoref_bam,n_mapped_reads));
  mapping_stats$mapped_reads_assemblyref<-unlist(lapply(mapping_stats$remapped_bam,n_mapped_reads));
  mapping_stats$num_Ns<-sum(letterFrequency(con_seq,c('N','+')));
  mapping_stats$width<-width(con_seq);
  mapping_stats$perc_Ns<-100*mapping_stats$num_Ns/mapping_stats$width;
  if(!dir.exists('./stats/')) dir.create('./stats/');
  write.csv(mapping_stats,file=paste('./stats/',sampname,'_mappingstats.csv',sep=''),row.names=F);

  return(TRUE)
}

##Get args from command line
args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
    print(args[[i]])
  }
}

##For testing (these args should come from command line)
## sampname='2016-01040_S451_L001'
## ref='NC_016842'
## remapped_bamfname

##Files, directories, target site
mapped_reads_folder<-'./mapped_reads/';

##Make consensus sequence--returns TRUE if this worked
conseq<-clean_consensus_hcov(sampname,remapped_bamfname,mappedtoref_bamfname,ref);

##Prepare seqs for annotation -- will make separate folders for A and B
if(conseq==TRUE){
  if(!dir.exists('./annotations_prokka')) dir.create('./annotations_prokka');

  ##Write consensus seq to folder for prokka
  fname<-paste('./consensus_seqs/',sampname,'.fasta',sep='')
  con_seq<-readDNAStringSet(fname);
  names(con_seq)<-substring(names(con_seq),1,20); ##prokka needs contig name to be <=20 chars long
  sampdir<-paste('./annotations_prokka/',sampname,sep='');
  if(!dir.exists(sampdir)) dir.create(sampdir); ##create folder for the sample
  writeXStringSet(con_seq,file=paste(sampdir,'/',sampname,'.fa',sep=''),format='fasta');

}else{
  print('Failed to generate consensus sequences.')
}
