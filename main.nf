nextflow.enable.dsl = 2

params.read1 = '/fusion/s3/ecd-rosalind/redsheet-data/output/1042608/pubdir/TID_1919_P_NC-Pool-2503_1_2_1_1_S5_L001_R1.trimmed.fq.gz'
params.read2 = '/fusion/s3/ecd-rosalind/redsheet-data/output/1042608/pubdir/TID_1919_P_NC-Pool-2503_1_2_1_1_S5_L001_R2.trimmed.fq.gz'
params.container_bwa_samtools = '297036008099.dkr.ecr.us-west-2.amazonaws.com/ecd-rosalind:bwa-samtools'
params.reference = '/fusion/s3/ecd-rosalind/redsheet-data/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa'
process bwamem1 {
    container params.container_bwa_samtools
    cpus 24
    memory '36 GB'
    input:
    tuple path(read1), path(read2), file(reference), file(reference_bwt), file(reference_pac), file(reference_sa), file(reference_amb), file(reference_ann)
    script:
    """
    echo bwa-mem1
    bwa mem -t 32 \
        $reference $read1 $read2 | samtools view -hb | samtools sort -o test.sorted.bam
    echo ---
    """
}
process bwamem2 {
    container params.container_bwa_samtools
    cpus 24
    memory '36 GB'
    input:
    tuple path(read1), path(read2), file(reference), file(reference_bwt), file(reference_pac), file(reference_sa), file(reference_amb), file(reference_ann)
    script:
    """
    echo bwa-mem2
    bwa mem -t 32 \
        $reference $read1 $read2 | samtools view -@ 12 -hb | samtools sort -@ 12 -o test.sorted.bam
    echo ---
    """
}
workflow {
    reference = file(params.reference)
    reference_bwt = file(params.reference+'.bwt')
    reference_pac = file(params.reference+'.pac')
    reference_sa = file(params.reference+'.sa')
    reference_amb = file(params.reference+'.amb')
    reference_ann = file(params.reference+'.ann')
    
    bwamem1([params.read1, params.read2, reference, reference_bwt, reference_pac, reference_sa, reference_amb, reference_ann])
    bwamem2([params.read1, params.read2, reference, reference_bwt, reference_pac, reference_sa, reference_amb, reference_ann])
}