process fastp {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.{gz,json}",mode: "copy" 
    
    label 'high'
    
    container 'quay.io/biocontainers/nanoq:0.2.1--h7d875b9_0'

    tag {sample_id}
    
    cpus 16
    memory '64.GB'

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: reads
    path "${sample_id}.fastp.json", emit: json
    
    script:
    """
    hostname > hostname
    fastp -w ${task.cpus} -z 6 \
    -i ${reads[0]} -I ${reads[1]} \
    -o ${sample_id}_R1.fastq.gz -O ${sample_id}_R2.fastq.gz \
    -j ${sample_id}.fastp.json
    """
}
    // --detect_adapter_for_pe \
process filtlong {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.gz",mode: "copy"
    
    tag {sample_id}
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtlong.fastq.gz")
    
    script:
    """
    filtlong --min_length $params.min_length --keep_percent 95 ${reads} | gzip > ${sample_id}_filtlong.fastq.gz 
    """
}


process nanoq {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.gz",mode: "copy"
    
    container 'quay.io/biocontainers/nanoq:0.2.1--h7d875b9_0'

    tag {sample_id}
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_fil.fastq.gz")
    
    script:
    //Assume inputs are gzip
    """
    zcat ${reads} | nanoq -l ${params.min_length} -q ${params.min_qscore} | gzip - > ${sample_id}_fil.fastq.gz
    """
}

process porechop {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.gz",mode: "copy"
    
    tag {sample_id}
    
    cpus 16

    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_porechop.fastq.gz")
    
    script:
    """
    porechop -t ${task.cpus} --discard_middle -i ${reads} -o ${sample_id}_porechop.fastq.gz
    """
}

process multiqc {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "*.html",mode: "copy"
    
    tag "Multiqc"
    
    errorStrategy 'ignore'

    input:
    path(input)
  
    output:
    path "multiqc_report.html" optional true
  
    script:
    """
    multiqc --interactive .  
    """
}

process checkm {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 16
    memory '480 GB'
    input:
    path(bins)
    
    output:
    path "checkm_output", emit: checkm
    path "*.tsv"

    script:
    """
    checkm lineage_wf -t ${task.cpus} -x fasta --tab_table . checkm_output > checkm_output.tsv
    grep -v -e "\\[" checkm_output.tsv | csvtk -t cut -f 1,2,3,12-14 > checkm_output_mqc.tsv
    """
}

process busco {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: {filename -> if (filename.endsWith(".tsv")) { 
        "./reports/${filename}"
    } else { 
        "./${filename}"
        }
    }
    
    container 'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0'

    errorStrategy 'ignore'

    tag {sample_id}
    
    cpus 8
    
    memory '16.GB'
    
    input:
    tuple val(sample_id), path(contigs)
    
    output:
    path("${sample_id}")
    path("${sample_id}_busco.tsv"), emit: report
    script:
    db = params.busco_db_location ? "--download_path ${params.busco_db_location}" : ""
    """
    export TMPDIR=.;
    busco \
    -i ${contigs} \
    -o ${sample_id} \
    --auto-lineage-prok \
    -m geno \
    -c ${task.cpus} \
    $db \
    grep -e "C:" ${sample_id}/short_summary.specific* > ${sample_id}_busco
    awk '{print \$0="${sample_id}"\$0}' ${sample_id}_busco > ${sample_id}_busco.tsv
    """
}


process quast {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {"running"}
    
    container 'quay.io/biocontainers/quast:5.0.2--py36pl5262h30a8e3e_4'
    
    cpus 8

    input:
    path(contigs)
    output:
    path("quast_output")
    
    script:
    """
    quast.py --min-contig 100 -t ${task.cpus} -o quast_output ${contigs}
    """
}

process sendsketch {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 8

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}.sketch")

    script:
    """
    sendsketch.sh in=${contigs} out=${sample_id}.sketch printrfname=t
    """
}

process gunc {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12
    memory '200.GB'

    input:
    tuple val(sample_id), path("contigs.fa")
    
    output:
    path("${sample_id}")
    script:
    """
    mkdir ${sample_id}
    ln -sf contigs.fa ${sample_id}.fa
    gunc run --threads ${task.cpus} \
    -i ${sample_id}.fa \
    --use_species_level \
    --out_dir ${sample_id} \
    --detailed_output \
    --sensitive \
    -r ${params.gunc_db}
    """
}
