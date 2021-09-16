process bwa {
    if (params.bwa_output) {
        publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode:"copy"
    }
    
    
    errorStrategy 'ignore'

    tag {sample_id}
    
    cpus 16

    input:
    tuple val(sample_id), path(reads), path(ref)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    tuple val(sample_id), path("${ref}") , emit: ref
    tuple val(sample_id), path(reads), emit: reads

    script:
    """
    bwa index ${ref}
    bwa mem -t ${task.cpus} ${ref} $reads | samtools view -Sb - | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam - 
    """
}

process snippy {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode:"copy"
    
    errorStrategy 'ignore'

    tag {sample_id}
    
    cpus 8

    input:
    
    tuple val(sample_id), path(reads), path(ref)
    
    output:
    path "${sample_id}"

    script:
    """
    hostname > hostname
    snippy --R1 ${reads[0]} --R2 ${reads[1]} --outdir ${sample_id} --cpus $task.cpus --ref ${ref}
    """
}

process snippy_core {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode:"copy"
    
    errorStrategy 'ignore'
    
    cpus 16

    input:
    path(ref)
    path(snp)
    
    output:
    path "core*"

    script:
    """
    snippy-core --ref ${ref} --prefix core ${snp}
    """
}