nextflow.enable.dsl=2

// ===== PARAMETERS =====
params.samplesheet = params.samplesheet ?: 'samples.csv'
params.adapters    = '/home/raghu/DGE_work/TruSeq3-PE.fa'
params.index       = '/home/raghu/DGE_work/GRCh38_index'
params.gtf         = '/home/raghu/DGE_work/gencode.v46.annotation.gtf'
params.output      = '/home/raghu/DGE_work/nextflow_output'
params.threads     = 8

// -------- PROCESSES --------
process SetupDirectories {
    output:
    path "setup_complete.txt"

    script:
    """
    mkdir -p ${params.output}/{fastqc,trimmed,hisat2,bam,counts}
    echo "Setup completed at \$(date)" > setup_complete.txt
    """
}

process FastQC {
    tag "$sample_id"
    publishDir "${params.output}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.{html,zip}")

    script:
    """
    fastqc --threads ${params.threads} -o . ${read1} ${read2}
    """
}

process Trimmomatic {
    tag "$sample_id"
    publishDir "${params.output}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fq.gz"), path("${sample_id}_2_trimmed.fq.gz")

    script:
    """
    trimmomatic PE -threads ${params.threads} \\
        ${read1} ${read2} \\
        ${sample_id}_1_trimmed.fq.gz ${sample_id}_1_unpaired.fq.gz \\
        ${sample_id}_2_trimmed.fq.gz ${sample_id}_2_unpaired.fq.gz \\
        ILLUMINACLIP:${params.adapters}:2:30:10 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:36
    """
}

process FastQC_Trimmed {
    tag "$sample_id"  
    publishDir "${params.output}/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*.{html,zip}")

    script:
    """
    fastqc --threads ${params.threads} -o . ${read1} ${read2}
    """
}

process HISAT2 {
    tag "$sample_id"
    publishDir "${params.output}/hisat2", mode: 'copy', pattern: "*.hisat2_mapstats.txt"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}.hisat2_mapstats.txt")

    script:
    """
    hisat2 -p ${params.threads} \\
        --very-sensitive \\
        -x ${params.index} \\
        -1 ${read1} -2 ${read2} \\
        --rna-strandness RF \\
        --new-summary \\
        -S ${sample_id}.sam \\
        2> ${sample_id}.hisat2_mapstats.txt
    """
}

process SamtoolsSort {
    tag "$sample_id"
    publishDir "${params.output}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file), path(mapstats)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools view -@ ${params.threads} -bS ${sam_file} | \\
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    rm -f ${sam_file}
    """
}

process FeatureCounts {
    publishDir "${params.output}/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "featureCounts.txt"
    path "featureCounts.txt.summary"

    script:
    """
    featureCounts -T ${params.threads} \\
        -a ${params.gtf} \\
        -o featureCounts.txt \\
        -p --countReadPairs -B -C \\
        -t exon -g gene_id \\
        ${bam_files}
    """
}

workflow {
    samples = Channel
        .fromPath("samples.csv")
        .splitCsv(header:true)
        .map { row ->
            tuple(row.sample_id, file(row.read1), file(row.read2))
        }

    setup = SetupDirectories()

    qc_raw       = samples | FastQC
    trimmed_reads = samples | Trimmomatic
    qc_trimmed   = trimmed_reads | FastQC_Trimmed
    aligned      = trimmed_reads | HISAT2
    sorted_bams  = aligned | SamtoolsSort

    bam_channel = sorted_bams.map { id, bam, bai -> bam }
    featurecounts_out = FeatureCounts(bam_channel)
}




