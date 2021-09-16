params.pilon_output = false
params.memory = '32G'

process medaka {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    label 'medaka'
    
    tag {sample_id}
    
    container 'quay.io/biocontainers/medaka:1.4.3--py38h130def0_0'

    cpus 12

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}_medaka.fasta")
    
    script:
    """
    medaka_consensus -i ${reads} -d ${contigs} -o output -t ${task.cpus} -m ${params.medaka_model}
    mv output/consensus.fasta ${sample_id}_medaka.fasta
    """
}

process pilon_raw {
    publishDir "${params.outdir}", mode: "copy",
    saveAs: {filename -> if (params.pilon_output) { 
                "${task.process.replaceAll(":","_")}/${filename}"
                } 
            else {
                null
                }
            }

    tag {sample_id}
    
    cpus 8
    // memory params.memory

    input:
    tuple val(sample_id), path(short_reads), path(contigs) 

    output:
    tuple val(sample_id), path("${sample_id}_pilon.fasta"), emit: contigs
    path("*.changes")

    script:
    """
    export _JAVA_OPTIONS="-Xmx${params.memory} -Xms512m"
    bwa index ${contigs}
    
    bwa mem -M -t ${task.cpus} ${contigs} ${short_reads} | samtools view -bS -| samtools sort > alignments.bam
   
    samtools index alignments.bam

    pilon --genome ${contigs} --frags alignments.bam --changes \
    --output ${sample_id}_pilon --fix all
    """
}


process pilon {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12
    
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    path "${sample_id}_pilon.fasta"
    path "*.changes", optional: true
    
    script:
    """
    pilon.py -t ${task.cpus} -n 10 ${reads[0]} ${reads[1]} ${contigs}
    mv final.polished.fasta  ${sample_id}_pilon.fasta
    """
}
