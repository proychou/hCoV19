#!/bin/bash
#Pipeline for whole genome sequence assembly and annotation for COVID 
#Feb 2020
#Pavitra Roychoudhury

#### Load required modules ####
#do this before submitting the sbatch command
# cd /fh/fast/jerome_k/COVID19_WGS/ 
# SLURM_CPUS_PER_TASK=14
# module load BBMap/38.44-foss-2016b
# module load FastQC/0.11.8-Java-1.8
# module load BWA/0.7.17-foss-2016b
# module load bowtie2/2.2.5
# module load SAMtools/1.8-foss-2016b 
# module load R/3.6.1-foss-2016b-fh1
# module load prokka/1.13-foss-2016b-BioPerl-1.7.0
# module load parallel/20170222-foss-2016b
# module load BLAST+/2.7.1-foss-2016b
# wget https://github.com/tseemann/prokka/raw/master/binaries/linux/tbl2asn -O $HOME/.local/bin/tbl2asn

#### One time: First build reference for bowtie and make a copy of the ref seqs ####
# module load bowtie2/2.2.5
# module load BWA/0.7.17-foss-2016b
# module load prokka/1.13-foss-2016b-BioPerl-1.7.0
# bowtie2-build './refs/NC_045512.2.fasta' ./refs/NC_045512.2
# bwa index ./refs/NC_045512.2.fasta
# prokka-genbank_to_fasta_db ./refs/NC_045512.2.gb > ./refs/hcov_proteins.faa


#### Usage ####
#For paired-end library
#		covid_wgs_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz \
#								-s samplename -aqpf
#For single-end library
#		covid_wgs_pipeline.sh -u yourreads.fastq.gz -s samplename -aqpf
#This is meant to be run on the cluster (typically through sbatch) so if run locally,
#first set the environment variable manually, e.g.
#		SLURM_CPUS_PER_TASK=8
#or whatever is the number of available processors 
 
#Load required tools
#Note that spades is locally installed and need to be updated manually as required


PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.14.0-Linux/bin:
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK
# echo "Path: "$PATH

while getopts ":1:2:u:s:fapq" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		u) in_fastq="$OPTARG"
			paired="false"
		;;
		f) filter="true"
		;;
		a) adapter_trim="true"
		;;
		q) qual_trim="true"
		;;
		p) primer_trim="true"
		;;
		s) sampname="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@

ref_fasta='./refs/NC_045512.2.fasta'
ref_bowtie='NC_045512.2'

##  PAIRED-END  ##
## paired end not tested as everything right now is single-end 
if [[ $paired == "true" ]]
then


if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi


#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK  


#Adapter trimming with bbduk
if [[ $adapter_trim == "true" ]]
then
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq1='./preprocessed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz'
tmp_fastq2='./preprocessed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed_r2.fastq.gz'

bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2 out1=$tmp_fastq1 out2=$tmp_fastq2 ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 

bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2 ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm $tmp_fastq1 $tmp_fastq2

else
processed_fastq1=$in_fastq_r1 
processed_fastq2=$in_fastq_r2
fi


# #Quality trimming
if [[ $qual_trim == "true" ]]
then
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
processed_fastq_old1=$processed_fastq1
processed_fastq_old2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_preprocessed_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_preprocessed_r2.fastq.gz'

bbduk.sh in1=$processed_fastq_old1 in2=$processed_fastq_old2 out1=$processed_fastq1 out2=$processed_fastq2 t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

fi


#Primer trimming: settings based on discussions in SPHERES consortium-- this assumes longer reads, so setting minimum read length to 75 if using primer trimming
if [[ $primer_trim == "true" ]]
then
printf "\n\nPrimer trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq1=$processed_fastq1
tmp_fastq2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed2_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed2_r2.fastq.gz'

bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2  ref=./refs/swift_primers.fasta k=18 ktrim=l hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictleft=30 t=$SLURM_CPUS_PER_TASK minlen=75

tmp_fastq1=$processed_fastq1
tmp_fastq2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed3_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed3_r2.fastq.gz'
bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2 ref=./refs/swift_primers.fasta k=18 ktrim=r hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictright=30 t=$SLURM_CPUS_PER_TASK minlen=75
rm $tmp_fastq1 $tmp_fastq2

fi




#Map reads to reference
printf "\n\nMapping reads to reference ... \n\n\n"
mkdir -p ./mapped_reads
mappedtoref_bam='./mapped_reads/'$sampname'.bam'
bowtie2 -x ./refs/$ref_bowtie -1 $processed_fastq1 -2 $processed_fastq2 -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > $mappedtoref_bam
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 
rm $mappedtoref_bam 
mv './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 


#Use bbduk to filter viral reads 
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hcov_refs.fasta ... \n\n\n"
processed_fastq_old1=$processed_fastq1
processed_fastq_old2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_matched_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_matched_r2.fastq.gz'
unmatched_fastq1='./filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' 
unmatched_fastq2='./filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' 
filter_stats='./preprocessed_fastq/'$sampname'_stats_filtering.txt'
bbduk.sh -Xmx80g in1=$processed_fastq_old1 in2=$processed_fastq_old2 out1=$unmatched_fastq1 out2=$unmatched_fastq2 outm1=$processed_fastq1 outm2=$processed_fastq2 ref=$ref_fasta k=31 hdist=2 stats=$filter_stats overwrite=TRUE t=$SLURM_CPUS_PER_TASK
fi
 

#FastQC report on processed reads
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
mkdir -p ./fastqc_reports_preprocessed
fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK $processed_fastq1 $processed_fastq2


#Assemble with SPAdes
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -1 $processed_fastq1 -2 $processed_fastq2 -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}




##  SINGLE-END  ## 
else 
if [[ $paired == "false" ]]
then

if [ -z $in_fastq ]
then
echo "Missing input argument."
fi

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq 

#Adapter trimming with bbduk
if [[ $adapter_trim == "true" ]]
then
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq='./preprocessed_fastq/'$sampname'_trimmed_tmp.fastq.gz'
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed.fastq.gz'

bbduk.sh in=$in_fastq out=$tmp_fastq ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in=$tmp_fastq out=$processed_fastq ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm $tmp_fastq

else
processed_fastq=$in_fastq 
fi

  
#Quality trimming
if [[ $qual_trim == "true" ]]
then
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
processed_fastq_old=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz'

bbduk.sh in=$processed_fastq_old out=$processed_fastq t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

fi


#Primer trimming -- this assumes longer reads, so setting minimum read length to 75 if using primer trimming
if [[ $primer_trim == "true" ]]
then
printf "\n\nPrimer trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed2.fastq.gz'

bbduk.sh in=$tmp_fastq out=$processed_fastq ref=/fh/fast/jerome_k/COVID19_WGS/refs/swift_primers.fasta k=18 ktrim=l hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictleft=30 t=$SLURM_CPUS_PER_TASK minlen=75

tmp_fastq=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed3.fastq.gz'
bbduk.sh in=$tmp_fastq out=$processed_fastq ref=/fh/fast/jerome_k/COVID19_WGS/refs/swift_primers.fasta k=18 ktrim=r hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictright=30 t=$SLURM_CPUS_PER_TASK minlen=75
rm $tmp_fastq

fi


#Use bbduk to filter viral reads 
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hcov_refs.fasta ... \n\n\n"
processed_fastq_old=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_matched.fastq.gz'
unmatched_fastq='./filtered_fastq/'$sampname'_unmatched.fastq.gz' 
filter_stats='./preprocessed_fastq/'$sampname'_stats_filtering.txt'
bbduk.sh -Xmx80g in=$processed_fastq_old out=$unmatched_fastq outm=$processed_fastq ref=$ref_fasta k=31 hdist=2 stats=$filter_stats overwrite=TRUE t=$SLURM_CPUS_PER_TASK
fi

#Map reads to reference
printf "\n\nMapping reads to reference ... \n\n\n"
mkdir -p ./mapped_reads
mappedtoref_bam='./mapped_reads/'$sampname'.bam'
bowtie2 -x ./refs/$ref_bowtie -U $processed_fastq -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > $mappedtoref_bam
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 
rm $mappedtoref_bam 
mv './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 


#FastQC report on processed reads
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
mkdir -p ./fastqc_reports_preprocessed
fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK $processed_fastq


#Assemble with SPAdes
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -s $processed_fastq -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

fi
fi


# Now call an R script that merges assembly and mapping and ultimately makes the consensus sequence 
scaffname='./contigs/'$sampname'/scaffolds.fasta'
Rscript --vanilla hcov_make_seq.R sampname=\"$sampname\" reffname=\"$ref_fasta\" scaffname=\"$scaffname\" ncores=\"$SLURM_CPUS_PER_TASK\"


#Remap reads to "new" reference
printf "\n\nRe-mapping reads to assembled sequence ... \n\n\n"
mkdir -p ./remapped_reads
remapping_btref='./ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref_bowtie

bowtie2-build -q './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref_bowtie'_consensus.fasta' $remapping_btref
remapped_bamfname='./remapped_reads/'$sampname'.bam'
if [[ $paired == "true" ]]
then
echo 'paired end not tested'
bowtie2 -x $remapping_btref -1 $processed_fastq1 -2 $processed_fastq2 -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > $remapped_bamfname
elif [[ $paired == "false" ]]
then
bowtie2 -x $remapping_btref -U $processed_fastq -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > $remapped_bamfname
fi


samtools sort -@ ${SLURM_CPUS_PER_TASK} -o './remapped_reads/'$sampname'.sorted.bam' $remapped_bamfname
rm $remapped_bamfname
mv './remapped_reads/'$sampname'.sorted.bam' $remapped_bamfname 


#Call R script to generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs
mkdir -p ./stats
Rscript --vanilla hcov_generate_consensus.R sampname=\"$sampname\" ref=\"$ref_bowtie\" remapped_bamfname=\"$remapped_bamfname\" mappedtoref_bamfname=\"$mappedtoref_bam\"

#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka
prokka --outdir './annotations_prokka/'$sampname'/' --force --kingdom 'Viruses' --genus 'Betacoronavirus' --usegenus --prefix $sampname --proteins ./refs/hcov_proteins.faa  './annotations_prokka/'$sampname'/'$sampname'.fa'


