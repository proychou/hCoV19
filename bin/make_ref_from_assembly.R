#!/usr/bin/env Rscript

library(Biostrings)
require(Rsamtools)
require(GenomicAlignments)
require(parallel)

make_ref_from_assembly<-function(bamfname, reffname, outfname){
  ## Make a new reference from scaffolds
  ## Pavitra Roychoudhury
  ## from https://github.com/proychou/ViralWGS/blob/master/wgs_functions.R

  ## bamfname - filename for bam of filtered scaffolds
  ## reffname - filename for reference fasta
  ## outfname - filename for output fasta

  ncores<-detectCores();

  ## Read reference sequence
  ref_seq<-readDNAStringSet(reffname);

  if(!is.na(bamfname)&class(try(scanBamHeader(bamfname),silent=T))!='try-error'){
    ## Index bam if required
    if(!file.exists(paste(bamfname,'.bai',sep=''))){
      baifname<-indexBam(bamfname);
    }else{
      baifname<-paste(bamfname,'.bai',sep='');
    }

    ## Import bam file
    params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                         what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq'));
    gal<-readGAlignments(bamfname,index=baifname,param=params);

    ## Remove any contigs with width <200 bp
    gal<-gal[width(gal)>200];

    ## First lay contigs on reference space--this removes insertions and produces a seq of the same length as ref
    qseq_on_ref<-sequenceLayer(mcols(gal)$seq,cigar(gal),from="query",to="reference");
    qseq_on_ref_aligned<-stackStrings(
        qseq_on_ref,1,max(mcols(gal)$pos+qwidth(gal)-1,width(ref_seq)),
        shift=mcols(gal)$pos-1,Lpadding.letter='N',Rpadding.letter='N');

    ## Make a consensus matrix and get a consensus sequence from the aligned scaffolds
    cm<-consensusMatrix(qseq_on_ref_aligned,as.prob=T,shift=0)[c('A','C','G','T','N','-'),];
                                        # cm[c('N','-'),]<-0;
    cm['N',]<-0;
    cm<-apply(cm,2,function(x)if(all(x==0))return(x) else return(x/sum(x)));
    cm['N',colSums(cm)==0]<-1;
    con_seq<-DNAStringSet(gsub('\\?','N',consensusString(cm,threshold=0.25)));
    con_seq<-DNAStringSet(gsub('\\+','N',con_seq));


    ## Now fill in the Ns with the reference
    temp<-as.matrix(con_seq);
    temp[temp=='N']<-as.matrix(ref_seq)[temp=='N'];
    con_seq<-DNAStringSet(paste0(temp,collapse=''));
    names(con_seq)<-sub('.bam','_consensus',basename(bamfname));

    ## Look for insertions in bam cigar string
    cigs_ref<-cigarRangesAlongReferenceSpace(cigar(gal),with.ops=F,ops='I',
                                             reduce.ranges=T,drop.empty.ranges=F,
                                             pos=mcols(gal)$pos);
    cigs_query<-cigarRangesAlongQuerySpace(cigar(gal),ops='I',with.ops=F,
                                           reduce.ranges=T,drop.empty.ranges=F);
    all_ins<-mclapply(c(1:length(cigs_query)),function(i)
      extractAt(mcols(gal)$seq[i],cigs_query[[i]])[[1]]);

    ## Merge all insertions
    all_ins_merged<-do.call('rbind',mclapply(c(1:length(cigs_ref)),function(i)
      return(data.frame(
          start_ref=start(cigs_ref[[i]]),end_ref=end(cigs_ref[[i]]),
          start_query=start(cigs_query[[i]]),end_query=end(cigs_query[[i]]),
          ins_seq=all_ins[[i]],width_ins=width(all_ins[[i]]))),
      mc.cores=ncores));
    all_ins_merged<-all_ins_merged[order(all_ins_merged$end_ref),];
    ## write.csv(all_ins_merged,'./testing/all_ins.csv',row.names=F);

    ## TO DO: Check for overlaps--should be minimal since scaffolds don't usually overlap that much
    if(any(table(all_ins_merged$start_ref)>1)){
      browser()

      ## Now the beauty part of inserting the strings back in
      ## Split ref seq by the insert positions
    }else if(nrow(all_ins_merged)!=0){
      new_strs<-DNAStringSet(rep('',nrow(all_ins_merged)+1))
      for(i in 1:nrow(all_ins_merged)){
        if(i==1){
          new_strs[i]<-paste0(
              extractAt(con_seq,IRanges(start=1,end=all_ins_merged$end_ref[i]))[[1]],
              all_ins_merged$ins_seq[i]);
        }else{
          new_strs[i]<-paste0(
              extractAt(con_seq,IRanges(start=all_ins_merged$start_ref[i-1],
                                        end=all_ins_merged$end_ref[i]))[[1]],
              all_ins_merged$ins_seq[i]);
        }
      }

      ## Last bit
      new_strs[i+1]<-paste0(
          extractAt(con_seq,IRanges(start=all_ins_merged$start_ref[i],
                                    end=width(con_seq)))[[1]])
      temp_str<-paste0(as.character(new_strs),collapse='');

      ## Remove gaps to get final sequence
      con_seq_final<-DNAStringSet(gsub('-','',temp_str));

      ## No insertions
    }else{
      con_seq_final<-con_seq;
    }

    writeXStringSet(con_seq_final, outfname);

    ## Delete bai file
    file.remove(baifname);

  }else{
    print('Bam file could not be opened.')
    return(NA)
  }
}

## keep command line argument parsing simple for now - argparse would be better
args <- commandArgs(TRUE)
make_ref_from_assembly(bamfname=args[1], reffname=args[2], outfname=args[3])
