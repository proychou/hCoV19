seqs = file('test/UWVL20030402-B01-VIRv1_S2_R1_001.fastq.gz')
reference_fa = file('refs/NC_045512.fasta')
genbank = file('refs/NC_045512.gb')

process bowtie_index {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('reference.fasta') from reference_fa

    output:
        file('*.bt2') into bowtie_ref

    """
    bowtie2-build reference.fasta reference
    """
}

process bwa_index {
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    label 'med_cpu_mem'

    input:
        file('reference.fasta') from reference_fa

    output:
        file('reference.*') into bwa_ref

    """
    bwa index reference.fasta
    """
}

process prokka_index {
    container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'

    label 'med_cpu_mem'

    input:
        file('reference.gb') from genbank

    output:
        file('proteins.faa') into prokka_ref

    """
    prokka-genbank_to_fasta_db reference.gb > proteins.faa
    """
}

// FastQC report on raw reads
process fastqc_raw_report  {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from seqs

    output:
       file('fastqc_reports_raw/*')

    publishDir params.output, overwrite: true

    """
    mkdir fastqc_reports_raw
    fastqc --outdir fastqc_reports_raw --threads 8 seqs.fastq.gz
    """
}

// Adapter trimming with bbduk
process right_trim {
    container 'quay.io/biocontainers/bbmap:38.79--h516909a_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from seqs

    output:
        file('right_trimmed.fastq.gz') into right_trimmed

    """
    bbduk.sh in=seqs.fastq.gz out=right_trimmed.fastq.gz hdist=2 k=21 ktrim=r mink=4 ref=adapters,artifacts
    """
}

process left_trim {
    container 'quay.io/biocontainers/bbmap:38.79--h516909a_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from right_trimmed

    output:
        file('left_trimmed.fastq.gz') into trimmed

    """
    bbduk.sh in=seqs.fastq.gz out=left_trimmed.fastq.gz hdist=2 k=21 ktrim=l mink=4 ref=adapters,artifacts
    """
}

// Quality trimming
process quality_trim {
    container 'quay.io/biocontainers/bbmap:38.79--h516909a_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from trimmed

    output:
        file('preprocessed.fastq.gz') into preprocessed

    """
    bbduk.sh in=seqs.fastq.gz out=preprocessed.fastq.gz maq=10 minlen=20 qtrim=rl trimq=20
    """
}

// Map reads to reference
// TODO: find out of bowtie2 can write directly to sorted bam file
process map_to_ref {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from preprocessed
        file('') from bowtie_ref

    output:
        file('alignment.sam') into sam

    """
    bowtie2 --threads 8 -x reference -U seqs.fastq.gz -S alignment.sam
    """
}

process sam_to_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.sam') from sam

    output:
        file('alignment.bam') into bam

    """
    samtools view -o alignment.bam -b alignment.sam
    """
}

process sort_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.bam') from bam

    output:
        file('sorted.bam') into sorted

    """
    samtools sort -o sorted.bam --threads 8 alignment.bam
    """
}

// Use bbduk to filter viral reads
process filter {
    container 'quay.io/biocontainers/bbmap:38.79--h516909a_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from preprocessed
        file('reference.fasta') from reference_fa

    output:
        file('matched.fastq.gz') into filtered
        file('unmatched.fastq.gz')
        file('stats_filtering.txt')

    """
    bbduk.sh in=seqs.fastq.gz out=unmatched.fastq.gz outm=matched.fastq.gz ref=reference.fasta hdist=2 k=31 stats=stats_filtering.txt
    """
}

// FastQC report on processed reads
process fastqc_processed_report  {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from filtered

    output:
       file('fastqc_reports_preprocessed/*')

    publishDir params.output, overwrite: true

    """
    mkdir fastqc_reports_preprocessed
    fastqc --outdir fastqc_reports_preprocessed --threads 8 seqs.fastq.gz
    """
}

// Assemble with SPAdesa
// FIXME: fix error here with some samples..
process assemble {
    container 'quay.io/biocontainers/spades:3.14.0--h2d02072_0'

    label 'med_cpu_mem'

    memory '8.GB' // 8.GB = 64 Gb

    input:
        file('seqs.fastq.gz') from filtered

    output:
        file('contigs/*') into assembled
        file('contigs/scaffolds.fasta') into scaffolds

    """
    mkdir contigs
    spades.py --careful --memory 64 --threads 8 -s seqs.fastq.gz -o contigs
    """
}

// First import scaffolds and filter by length (>200) and coverage (>10x)
process merge {
    // TODO: This file is a mini pipeline. Break out pieces into nf processes
    container 'hcov19-r-deps:latest'

    label 'med_cpu_mem'

    input:
        file('scaffolds.fasta') from scaffolds
        file('reference.fasta') from reference_fa
        file('') from bwa_ref

    output:
        file('ref_for_remapping/*') into consensus

    """
    hcov_make_seq.R sampname=\\"sample\\"  reffname=\\"reference.fasta\\" scaffname=\\"scaffolds.fasta\\" ncores=\"8\"
    """
}

process bowtie_consensus_build {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('consensus.fasta') from consensus

    output:
        file('*.bt2') into consensus_build

    """
    bowtie2-build -q consensus.fasta reference
    """
}

process map_to_consensus {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from preprocessed
        file('') from consensus_build

    output:
        file('alignment.sam') into remap_sam

    """
    bowtie2 --threads 8 -x reference -U seqs.fastq.gz -S alignment.sam
    """
}

process consensus_sam_to_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.sam') from remap_sam

    output:
        file('alignment.bam') into remap_bam

    """
    samtools view -o alignment.bam -b alignment.sam
    """
}

process sort_consensus_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.bam') from remap_bam

    output:
        file('remapped_sorted.bam') into remap_sorted

    """
    samtools sort -o remapped_sorted.bam --threads 8 alignment.bam
    """
}

process final_seq {
    container 'hcov19-r-deps:latest'

    label 'med_cpu_mem'

    input:
        file('remapped_sorted.bam') from remap_sorted
        file('sorted.bam') from sorted

    output:
        file('annotations_prokka/sample/sample.fa') into final_cons

    """
    hcov_generate_consensus.R sampname=\\"sample\\" ref=\\"reference\\" remapped_bamfname=\\"remapped_sorted.bam\\" mappedtoref_bamfname=\\"sorted.bam\\"
    """
}

process prokka_annnotations {
    container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'

    label 'med_cpu_mem'

    input:
        file('sample.fa') from final_cons
        file('proteins.faa') from prokka_ref

    publishDir params.output, overwrite: true

    output:
        file('*')

    """
    prokka --kingdom 'Viruses' --genus 'Betacoronavirus' --usegenus --prefix annotations --proteins proteins.faa sample.fa
    """
}
