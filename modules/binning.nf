process metabat {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: {filename -> "${sample_id}/${filename}"}
    
    label 'high'

    errorStrategy 'ignore'
    
    tag {sample_id}
    
    cpus 8

    input:
    tuple val(sample_id), path(reads), path(contigs), path(bam)

    output:
    path "${sample_id}.*.fa" , emit: bin
    
    script:
    """
    OMP_NUM_THREADS=${task.cpus} jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    metabat2 -t "${task.cpus}" -i "${contigs}" -a depth.txt -o "${sample_id}" -m ${params.metabat_min_size} --unbinned --seed ${params.metabat_rng_seed}
    """
}

process maxbin {
    publishDir 
    
    tag {}
    
    cpus 

    input:

    output:

    script:
    """
    """
}

process conccoct {
    publishDir 
    
    tag {}
    
    cpus 

    input:

    output:

    script:
    """
    """
}

process dastool {
    publishDir 
    
    tag {}
    
    cpus 

    input:

    output:

    script:
    """
    """
}