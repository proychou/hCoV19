// Author: Pavitra Roychoudhury
Channel.fromPath(params.manifest)
       .splitCsv(header: true)
       .take(params.take)
       .map{ [it['sample'], file(params.datadir + it['fastq'])] }
       .into{samples; for_raw_stats; for_prereport}
reference_fa = file("refs/NC_045512.fasta")
reference_gb = file("refs/NC_045512.gb")

// TODO: make separate pipeline for reference creation, store in s3
process bowtie_build {
    container "quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0"

    input:
        file(ref) from reference_fa

    output:
        file('*.bt2') into bowtie_ref

    """
    bowtie2-build --threads ${task.cpus} ${ref} reference
    """
}

process bwa_index {
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"

    input:
        file(ref) from reference_fa
    output:
        file('*') into bwa_ref

    """
    bwa index ${ref}
    """
}

process prokka_db {
    container "quay.io/biocontainers/prokka:1.14.6--pl526_0"

    input:
        file("reference.gb") from reference_gb
    output:
        file("proteins.faa") into prokka_ref

    """
    prokka-genbank_to_fasta_db reference.gb > proteins.faa
    """
}

// FastQC report on raw reads
// process fastqc_prereport  {
//     container "quay.io/biocontainers/fastqc:0.11.9--0"
//     label "fastqc_mem"
//     publishDir "${params.output}/fastqc/prereport/", mode: "copy", overwrite: true
//
//     input:
//         file("*.fastq.gz") from for_prereport.collect{ it[1] }
//     output:
//         file("*") into prereports
//
//     """
//     zcat *.fastq.gz | fastqc --threads ${task.cpus} stdin:report
//     """
// }

process generate_raw_stats {
    container "python:3.8.2-buster"

    input:
        tuple(val(sample), file(fastq)) from for_raw_stats

    output:
        file("raw_stats.csv") into raw_stats

    """
    raw_stats.py --out raw_stats.csv ${sample} ${fastq}
    """
}

// Adapter trimming with bbduk
process trimming {
    container "quay.io/biocontainers/bbmap:38.79--h516909a_0"

    input:
        tuple(val(sample), file(fastq)) from samples
    output:
        tuple(val(sample), file("trimmed.fastq.gz")) into (trimmed, for_reference_mapping, for_consensus_mapping)

    // bbduk --help: When piping interleaving must be explicitly stated: int=f unpaired, int=t for paired
    """
    df -h
    bbduk.sh in=${fastq} out=stdout.fq hdist=2 interleaved=f k=21 ktrim=r mink=4 ref=adapters,artifacts threads=${(task.cpus/3).intValue()} |
    bbduk.sh in=stdin.fq out=stdout.fq hdist=2 interleaved=f k=21 ktrim=l mink=4 ref=adapters,artifacts threads=${(task.cpus/3).intValue()} |
    bbduk.sh in=stdin.fq out=trimmed.fastq.gz interleaved=f maq=10 minlen=20 qtrim=rl trimq=20 threads=${(task.cpus/3).intValue()}
    """
}

// Map reads to reference
process map_to_ref {
    container "quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0"

    input:
        tuple(val(sample), file(fastq)) from for_reference_mapping
        file('') from bowtie_ref
    output:
        tuple(val(sample), file("alignment.sam")) into sam

    """
    df -h
    bowtie2 --threads ${task.cpus} -x reference -U ${fastq} -S alignment.sam
    """
}

process sam_to_bam {
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    input:
        tuple(val(sample), file("alignment.sam")) from sam
    output:
        tuple(val(sample), file("sorted.bam")) into (sorted, for_counting)

    """
    samtools view --threads ${task.cpus} -b alignment.sam |
    samtools sort -m ${(task.memory.toMega()/task.cpus).intValue()}M --threads ${task.cpus} -o sorted.bam
    """
}

process mapped_reads_ref {
    container "quay.io/biocontainers/pysam:0.15.4--py37hbcae180_0"

    input:
        tuple(val(sample), file(bam)) from for_counting

    output:
        file("count.csv") into mapped_reads_ref_counts

    """
    mapped_reads_ref.py --out count.csv ${sample} ${bam}
    """
}

// Use bbduk to filter viral reads
process filter_viral {
    container "quay.io/biocontainers/bbmap:38.79--h516909a_0"

    input:
        tuple(val(sample), file(fastq)) from trimmed
        file(ref) from reference_fa
    output:
        tuple(val(sample), file("viral.fastq.gz")) into processed
        file("viral.fastq.gz") into for_processed_report
        file("unmatched.fastq.gz")
        file("stats_filtering.txt")

    """
    bbduk.sh in=${fastq} out=unmatched.fastq.gz outm=viral.fastq.gz ref=${ref} hdist=2 k=31 stats=stats_filtering.txt --threads=${task.cpus}
    """
}

// FastQC report on processed reads
// process fastqc_processed_report  {
//     container "quay.io/biocontainers/fastqc:0.11.9--0"
//     publishDir "${params.output}/fastqc/processed/", mode: "copy", overwrite: true
//
//     input:
//         file("*.fastq.gz") from for_processed_report.collect()
//     output:
//        file("*") into processed_reports
//
//     """
//     zcat *.fastq.gz | fastqc --threads ${task.cpus} stdin:report
//     """
// }

// Assemble with SPAdes
// FIXME: fix error here with some samples..
process assemble_scaffolds {
    container "quay.io/biocontainers/spades:3.14.0--h2d02072_0"

    input:
        tuple(val(sample), file(fastq)) from processed
    output:
        tuple(val(sample), file("scaffolds.fasta")) into scaffolds

    """
    spades.py --careful --memory ${task.memory.toGiga() * 8} --threads ${task.cpus} -s ${fastq} -o .
    """
}

// lifted scaffold filtering out of hcov_make_seq.R
process filter_scaffolds {
    container 'bioconductor/release_core2:R3.6.2_Bioc3.10'

    input:
        tuple(val(sample), file(scaffolds)) from scaffolds
    output:
        tuple(val(sample), file("filtered_scaffolds.fasta")) into filtered_scaffolds
        file("stats.csv") into scaffolds_stats

    //TODO: min_len and min_cov as params?
    //TODO: output mean coverage
    // TODO: raise exception if all scaffolds filtered
    """
    filter_scaffolds.R ${sample} ${scaffolds} filtered_scaffolds.fasta stats.csv 200 10
    """
}

// lifted scaffold alignment out of hcov_make_seq.R
process align_scaffolds {
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    input:
        tuple(val(sample), file(scaffold)) from filtered_scaffolds
        file(ref_fa) from reference_fa
        file('') from bwa_ref
    output:
        tuple(val(sample), file('filtered_scaffolds.sam')) into aligned_contigs

    """
    bwa mem ${ref_fa} ${scaffold} > filtered_scaffolds.sam
    """
}

process scaffolds_sam_to_bam {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    input:
        tuple(val(sample), file(contig_sam)) from aligned_contigs
        file(ref_fa) from reference_fa
    output:
        tuple(val(sample), file('filtered_scaffolds_sorted.bam')) into filtered_scaffold_bam

    // TODO: samtools sort --help: -m INT Set maximum memory per thread; suffix K/M/G recognized [768M]
    """
    samtools view -bh -@ ${task.cpus} -o filtered_scaffolds.bam \
       ${contig_sam} -T ${ref_fa}
    samtools sort -@ ${task.cpus} -o filtered_scaffolds_sorted.bam \
       filtered_scaffolds.bam
    """
}

process scaffolds_bam_to_consensus {
    container 'bioconductor/release_core2:R3.6.2_Bioc3.10'

    input:
        tuple(val(sample), file(bam)) from filtered_scaffold_bam
        file(ref_fa) from reference_fa
    output:
        tuple(val(sample), file('scaffold_consensus.fasta')) into consensus

    """
    make_ref_from_assembly.R ${bam} ${ref_fa} scaffold_consensus.fasta
    """
}

process bowtie_build_consensus {
    container "quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0"

    input:
        tuple(val(sample), file("consensus.fasta")) from consensus
    output:
        tuple(val(sample), file("*.bt2")) into consensus_build

    """
    bowtie2-build --threads ${task.cpus} consensus.fasta reference
    """
}

process map_to_consensus {
    container "quay.io/biocontainers/bowtie2:2.4.1--py38he513fc3_0"

    input:
        tuple(val(sample), file(fastq), file("")) from for_consensus_mapping.join(consensus_build)
    output:
        tuple(val(sample), file("alignment.sam")) into remap_sam

    """
    bowtie2 --threads ${task.cpus} -x reference -U ${fastq} -S alignment.sam
    """
}

process consensus_sam_to_bam {
    container "quay.io/biocontainers/samtools:1.10--h9402c20_2"

    input:
        tuple(val(sample), file("alignment.sam")) from remap_sam
    output:
        tuple(val(sample), file("remap_sorted.bam")) into remap_sorted

    """
    samtools view --threads ${task.cpus} -b alignment.sam |
    samtools sort -m ${(task.memory.toMega()/task.cpus).intValue()}M --threads ${task.cpus} -o remap_sorted.bam
    """
}

process final_consensus {
    container 'bioconductor/release_core2:R3.6.2_Bioc3.10'

    input:
        tuple(val(sample), file(mapped_bam), file(remapped_bam)) from sorted.join(remap_sorted)
    output:
        tuple(val(sample), file('final_consensus.fasta')) into final_cons
        file('mapping_stats.csv') into mapping_stats

    // FIXME: final_cons seqname must be sample name for prokka genbank LOCUS
    // TODO: add join(remainder: true) option and allow final_consensus without remap_sam
    """
    hcov_generate_consensus.R ${sample} ${remapped_bam} ${mapped_bam} ${reference_fa} \
        final_consensus.fasta mapping_stats.csv
    """
}

// TODO: check if prokka needs contigs to be >20 chars long
process prokka_annotations {
    container "quay.io/biocontainers/prokka:1.14.6--pl526_0"
    publishDir "${params.output}/${sample}/", mode:"copy", overwrite: true

    input:
        tuple(val(sample), file(fasta)) from final_cons
        file(ref) from prokka_ref
    output:
        val(sample) into submit
        file("genbank/*")

    """
    prokka --cpus ${task.cpus} --kingdom "Viruses" --genus "Betacoronavirus" --usegenus --outdir genbank --prefix ${sample} --proteins ${ref} ${fasta}
    """
}

process report {
    container "python:3.8.2-buster"
    publishDir params.output, mode: "copy", overwrite: true

    input:
        file('*.csv') from raw_stats.concat(
                             mapped_reads_ref_counts,
                             scaffolds_stats,
                             mapping_stats).collect()
    output:
        file("report.csv")

    """
    report.py --out report.csv sample,run,fastq,length,raw_read_count,mapped_reads_ref,mean_coverage,perc_Ns,analysis_date *.csv
    """
}

// process multiqc_report {
//     container "quay.io/biocontainers/multiqc:1.8--py_2"
//     publishDir params.output, mode: "copy", overwrite: true
//
//     input:
//         file("prereport/*") from prereports
//         file("post_report/*") from processed_reports
//     output:
//         file("multiqc_*")
//
//     """
//     multiqc prereport post_report
//     """
// }
