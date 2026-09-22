#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * ATAC-seq Nextflow Pipeline
 * Converted from the original Bash implementation.
 *
 * Default dataset: ERR15487529
 * Reference genome: hg38
 */

// -----------------------------
// Pipeline parameters
// -----------------------------

params.sample = "ERR15487529"
params.genome_index = "${launchDir}/reference/hg38"
params.genome_size = "2.7e9"
params.threads = 4

// Project root defaults to the directory from which Nextflow is launched
params.project_dir = "${launchDir}"

// Input / output directories
params.raw_data_dir = "${params.project_dir}/raw_data"
params.qc_dir = "${params.project_dir}/results/qc_reports"
params.trim_dir = "${params.project_dir}/results/trimmed_data"
params.align_dir = "${params.project_dir}/results/aligned_data"
params.filter_dir = "${params.project_dir}/results/filtered_data"
params.peaks_dir = "${params.project_dir}/results/peaks"
params.bigwig_dir = "${params.project_dir}/results/bigwig_tracks"
params.multiqc_dir = "${params.project_dir}/results/multiqc_report"


// -----------------------------
// Process 1: FastQC
// -----------------------------

process fastqc {

    publishDir "${params.qc_dir}", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}_*.zip", emit: fastqc_files

    script:
    """
    fastqc -o . ${reads}
    """
}


// -----------------------------
// Process 2: Trim Galore
// -----------------------------

process trimGalore {

    publishDir "${params.trim_dir}", mode: 'copy', pattern: "${sample}_{1,2}_val_{1,2}.fq.gz"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample),
          path("${sample}_1_val_1.fq.gz"),
          path("${sample}_2_val_2.fq.gz"),
          emit: trimmed_reads

    script:
    """
    trim_galore \
        --paired \
        --nextera \
        --cores ${params.threads} \
        ${reads[0]} ${reads[1]}
    """
}


// -----------------------------
// Process 3: Bowtie2 Alignment
// -----------------------------

process bowtie2 {

    publishDir "${params.align_dir}", mode: 'copy'

    input:
    tuple val(sample), path(trimmed_reads)

    output:
    path "${sample}.bam", emit: bam

    script:
    """
    bowtie2 \
        -x ${params.genome_index} \
        -1 ${trimmed_reads[0]} \
        -2 ${trimmed_reads[1]} \
        --very-sensitive-local \
        -p ${params.threads} \
        | samtools view -@ ${params.threads} -bS - \
        > ${sample}.bam
    """
}


// -----------------------------
// Process 4: SAMtools Filtering
// -----------------------------

process samtoolsFilter {

    publishDir "${params.filter_dir}", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.baseName}_clean.bam", emit: clean_bam
    path "${bam.baseName}_clean.bam.bai", emit: clean_bam_index

    script:
    """
    samtools sort \
        -@ ${params.threads} \
        -n \
        -o ${bam.baseName}_sorted_name.bam \
        ${bam}

    samtools view \
        -@ ${params.threads} \
        -b \
        -F 1804 \
        ${bam.baseName}_sorted_name.bam \
        | samtools view -@ ${params.threads} -h - \
        | grep -v 'chrM' \
        | samtools view -@ ${params.threads} -b - \
        | samtools sort \
            -@ ${params.threads} \
            -o ${bam.baseName}_clean.bam -

    samtools index \
        -@ ${params.threads} \
        ${bam.baseName}_clean.bam
    """
}


// -----------------------------
// Process 5: MACS2 Peak Calling
// -----------------------------

process macs2 {

    publishDir "${params.peaks_dir}", mode: 'copy'

    input:
    path clean_bam

    output:
    path "${clean_bam.baseName}_ATAC_peaks.*", emit: peaks

    script:
    """
    macs2 callpeak \
        -t ${clean_bam} \
        -f BAMPE \
        -g ${params.genome_size} \
        -n ${clean_bam.baseName}_ATAC \
        -q 0.05 \
        --nomodel \
        --shift -100 \
        --extsize 200 \
        --keep-dup all
    """
}


// -----------------------------
// Process 6: DeepTools
// -----------------------------

process bamCoverage {

    publishDir "${params.bigwig_dir}", mode: 'copy'

    input:
    path clean_bam
    path clean_bam_index

    output:
    path "${clean_bam.baseName}_coverage.bw", emit: bigwig

    script:
    """
    bamCoverage \
        -b ${clean_bam} \
        -o ${clean_bam.baseName}_coverage.bw \
        --normalizeUsing CPM \
        --binSize 10 \
        -p ${params.threads} \
        --ignoreForNormalization chrX chrY chrM
    """
}


// -----------------------------
// Process 7: MultiQC
// -----------------------------

process multiqc {

    publishDir "${params.multiqc_dir}", mode: 'copy'

    input:
    path fastqc_files
    path trimmed_reads
    path peaks

    output:
    path "multiqc_report/*", emit: multiqc_report

    script:
    """
    multiqc . -o multiqc_report
    """
}


// -----------------------------
// Main Workflow
// -----------------------------

workflow {

    /*
     * Input paired-end FASTQ files.
     * Expected naming convention:
     *   {sample}_1.fastq.gz
     *   {sample}_2.fastq.gz
     */

    fastq_ch = Channel.fromFilePairs(
        "${params.raw_data_dir}/${params.sample}_{1,2}.fastq.gz"
    )

    fastq_ch.subscribe { sample, files ->

        if (!files[0].exists() || !files[1].exists()) {
            error "Input files ${files[0]} or ${files[1]} not found."
        }

    }

    // Quality control
    fastqc(fastq_ch)

    // Adapter trimming
    trimGalore(fastq_ch)

    trimmed_reads_ch = trimGalore.out.trimmed_reads.map {
        sample, read1, read2 ->
        [sample, [read1, read2]]
    }

    // Alignment
    bowtie2(trimmed_reads_ch)

    // Filtering
    samtoolsFilter(bowtie2.out.bam)

    // Peak calling
    macs2(samtoolsFilter.out.clean_bam)

    // Coverage track generation
    bamCoverage(
        samtoolsFilter.out.clean_bam,
        samtoolsFilter.out.clean_bam_index
    )

    // Collect files for MultiQC
    trimmed_reads_files = trimGalore.out.trimmed_reads.map {
        sample, read1, read2 ->
        [read1, read2]
    }.collect()

    multiqc(
        fastqc.out.fastqc_files.collect(),
        trimmed_reads_files,
        macs2.out.peaks.collect()
    )

    println """
    Pipeline complete.

    Peaks:      ${params.peaks_dir}
    BigWig:     ${params.bigwig_dir}
    MultiQC:    ${params.multiqc_dir}
    """
}


// -----------------------------
// Optional Trim Galore test workflow
// -----------------------------

workflow trimGaloreTest {

    fastq_ch = Channel.fromFilePairs(
        "${params.raw_data_dir}/${params.sample}_{1,2}.fastq.gz"
    )

    fastq_ch.subscribe { sample, files ->

        if (!files[0].exists() || !files[1].exists()) {
            error "Input files ${files[0]} or ${files[1]} not found."
        }

    }

    trimGalore(fastq_ch)

    emit:
    trimGalore.out.trimmed_reads
}
