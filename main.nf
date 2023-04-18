nextflow.enable.dsl = 2

params.read1 = '/fusion/s3/ecd-rosalind/redsheet-data/output/1042608/pubdir/TID_1919_P_NC-Pool-2503_1_2_1_1_S5_L001_R1.trimmed.fq.gz'
params.read2 = '/fusion/s3/ecd-rosalind/redsheet-data/output/1042608/pubdir/TID_1919_P_NC-Pool-2503_1_2_1_1_S5_L001_R2.trimmed.fq.gz'
params.container_bwa_samtools = '297036008099.dkr.ecr.us-west-2.amazonaws.com/ecd-rosalind:bwa-samtools'
params.reference = '/fusion/s3/ecd-rosalind/redsheet-data/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa'
process bwamem1 {
    container params.container_bwa_samtools
    cpus 24
    memory '36 GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(read1), file(read2), file(reference), file(reference_bwt), file(reference_pac), file(reference_sa), file(reference_amb), file(reference_ann)
    output:
    tuple val(sampleID), file("*bam")
    script:
    """
    echo bwa-mem ${sampleID}
    bwa mem -t 32 \
        $reference $read1 $read2 | samtools view -hb | samtools sort -o ${sampleID}.sorted.bam
    echo ---
    """
}
process bwamem2 {
    container params.container_bwa_samtools
    cpus 24
    memory '36 GB'
    tag {"$sampleID"}
    input:
    tuple val(sampleID), file(read1), file(read2), file(reference), file(reference_bwt), file(reference_pac), file(reference_sa), file(reference_amb), file(reference_ann)
    output:
    tuple val(sampleID), file("*bam")
    script:
    """
    echo bwa-mem ${sampleID}
    bwa mem -t 32 \
        $reference $read1 $read2 | samtools view -@ 12 -hb | samtools sort -@ 12 -o ${sampleID}.sorted.bam
    echo ---
    """
}
workflow {
    bwamem1(fastp.out.trimmed_fqs.map { row -> row +[reference, reference_bwt, reference_pac, reference_sa, reference_amb, reference_ann]})
    bwamem2(fastp.out.trimmed_fqs.map { row -> row +[reference, reference_bwt, reference_pac, reference_sa, reference_amb, reference_ann]})
}