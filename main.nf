// TODO: Channel [sample_name, fastq]
Channel.fromPath('test/*.fastq.gz').into{seqs; for_reporting}
reference_fa = file('refs/NC_045512.fasta')
reference_gb = file('refs/NC_045512.gb')

// TODO: make separate pipeline for reference creation, store in s3
process bowtie_index {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('reference.fasta') from reference_fa

    output:
        file('*.bt2') into bowtie_ref

    """
    bowtie2-build --threads ${task.cpus} reference.fasta reference
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
        file('reference.gb') from reference_gb

    output:
        file('proteins.faa') into prokka_ref

    """
    prokka-genbank_to_fasta_db reference.gb > proteins.faa
    """
}

// FastQC report on raw reads
process fastqc_raw_report  {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    label 'fastqc_mem'

    input:
        file('seqs_*.fastq.gz') from for_reporting.collect()

    output:
       file('fastqc_reports_raw/*')

    publishDir params.output, overwrite: true

    """
    mkdir fastqc_reports_raw
    fastqc --outdir fastqc_reports_raw --threads ${task.cpus} seqs_*.fastq.gz
    """
}

// Adapter trimming with bbduk
process preprocess {
    container 'quay.io/biocontainers/bbmap:38.79--h516909a_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from seqs

    output:
        file('preprocessed.fastq.gz') into (preprocessed, for_reference_mapping, for_consensus_mapping)

    """
    bbduk.sh in=seqs.fastq.gz out=right_trimmed.fastq.gz hdist=2 k=21 ktrim=r mink=4 ref=adapters,artifacts threads=${task.cpus}
    bbduk.sh in=right_trimmed.fastq.gz out=trimmed.fastq.gz hdist=2 k=21 ktrim=l mink=4 ref=adapters,artifacts threads=${task.cpus}
    bbduk.sh in=trimmed.fastq.gz out=preprocessed.fastq.gz maq=10 minlen=20 qtrim=rl trimq=20 threads=${task.cpus}
    """
}

// Map reads to reference
process map_to_ref {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from for_reference_mapping
        file('') from bowtie_ref

    output:
        file('alignment.sam') into sam

    """
    bowtie2 --threads ${task.cpus} -x reference -U seqs.fastq.gz -S alignment.sam
    """
}

process sam_to_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.sam') from sam

    output:
        file('sorted.bam') into sorted

    """
    samtools view --threads ${task.cpus} -b alignment.sam -o alignment.bam
    samtools sort --threads ${task.cpus} -o sorted.bam alignment.bam
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
        file('matched.fastq.gz') into (filtered, for_reporting2)
        file('unmatched.fastq.gz')
        file('stats_filtering.txt')

    """
    bbduk.sh in=seqs.fastq.gz out=unmatched.fastq.gz outm=matched.fastq.gz ref=reference.fasta hdist=2 k=31 stats=stats_filtering.txt --threads=${task.cpus}
    """
}

// FastQC report on processed reads
process fastqc_processed_report  {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    label 'fastqc_mem'

    input:
        file('seqs_*.fastq.gz') from for_reporting2.collect()

    output:
       file('fastqc_reports_preprocessed/*')

    publishDir params.output, overwrite: true

    """
    mkdir fastqc_reports_preprocessed
    fastqc --outdir fastqc_reports_preprocessed --threads ${task.cpus} seqs_*.fastq.gz
    """
}

// Assemble with SPAdesa
// FIXME: fix error here with some samples..
process assemble {
    container 'quay.io/biocontainers/spades:3.14.0--h2d02072_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from filtered

    output:
        // file('contigs/*') into assembled
        file('contigs/scaffolds.fasta') into scaffolds

    """
    mkdir contigs
    spades.py --careful --memory ${task.memory.toGiga() * 8} --threads ${task.cpus} -s seqs.fastq.gz -o contigs
    """
}

// First import scaffolds and filter by length (>200) and coverage (>10x)
// TODO: This file is a mini pipeline. Break out pieces into nf processes, bwa and samtools
// TODO: Make coverage a parameter
process merge {
    container 'hcov19-r-deps:latest'

    label 'med_cpu_mem'

    input:
        file('scaffolds.fasta') from scaffolds
        file('reference.fasta') from reference_fa
        file('') from bwa_ref

    output:
        file('ref_for_remapping/*') into consensus

    """
    hcov_make_seq.R sampname=\\"sample\\"  reffname=\\"reference.fasta\\" scaffname=\\"scaffolds.fasta\\" ncores=\"${task.cpus}\"
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
    bowtie2-build --threads ${task.cpus} consensus.fasta reference
    """
}

process map_to_consensus {
    container 'quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0'

    label 'med_cpu_mem'

    input:
        file('seqs.fastq.gz') from for_consensus_mapping
        file('') from consensus_build

    output:
        file('alignment.sam') into remap_sam

    """
    bowtie2 --threads ${task.cpus} -x reference -U seqs.fastq.gz -S alignment.sam
    """
}

process consensus_sam_to_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    label 'med_cpu_mem'

    input:
        file('alignment.sam') from remap_sam

    output:
        file('sorted.bam') into remap_sorted

    """
    samtools view --threads ${task.cpus} -b alignment.sam -o alignment.bam
    samtools sort --threads ${task.cpus} -o sorted.bam alignment.bam
    """
}

// TODO: paste imported R function directly into this file
process final_seq {
    container 'hcov19-r-deps:latest'

    label 'med_cpu_mem'

    input:
        file('remapped.bam') from remap_sorted
        file('mapped.bam') from sorted

    output:
        file('annotations_prokka/sample/sample.fa') into final_cons

    """
    hcov_generate_consensus.R sampname=\\"sample\\" ref=\\"reference\\" remapped_bamfname=\\"remapped.bam\\" mappedtoref_bamfname=\\"mapped.bam\\"
    """
}

// TODO: check if prokka needs contigs to be >20 chars long
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
    prokka --cpus ${task.cpus} --kingdom 'Viruses' --genus 'Betacoronavirus' --usegenus --prefix annotations --proteins proteins.faa sample.fa
    """
}

// TODO: add ncbi genbank submissions manifest
