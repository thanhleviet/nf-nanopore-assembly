//Default params
params.genome_size = false
params.flye_options = false

process megahit {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
   saveAs: {filename -> if (filename.endsWith(".contigs.fa.gz") || filename.endsWith(".log")) "$filename"
                          else null}

  label 'high'
  
  errorStrategy 'retry'
  
  maxRetries params.maxRetries
  
  tag {sample_id}
  
  cpus 8
  
  input:
  tuple val(sample_id), path(reads)
  output:
  path "${sample_id}.log"
  tuple val(sample_id), path("${sample_id}.contigs.fa.gz"), emit: contigs
  tuple val(sample_id), path(reads), emit: reads
 
  script:

  """
  megahit -t $task.cpus -o megahit --out-prefix ${sample_id} -1 ${reads[0]} -2 ${reads[1]}
  gzip -c "megahit/${sample_id}.contigs.fa" > "${sample_id}.contigs.fa.gz"
  mv megahit/${sample_id}.log ${sample_id}.log 
  """
}

process flye {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    errorStrategy 'ignore'
    
    container 'quay.io/biocontainers/flye:2.9--py27h6a42192_0'

    tag {sample_id}
    
    cpus 12

    input:
    tuple val(sample_id), path(reads), val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: contigs
    tuple path("${sample_id}_assembly_graph.gfa"), path("${sample_id}_assembly_info.txt"), path("${sample_id}_flye.log"), path("${sample_id}_assembly_stats.txt")
    
    script:
    genome_size = params.genome_size ? "--genome-size ${params.genome_size}" : ""
    """
    flye -t ${task.cpus} -o flye_output ${params.flye_input} ${reads} --plasmids ${genome_size} ${params.flye_options}
    mv flye_output/assembly.fasta ${sample_id}_assembly.fasta
    mv flye_output/assembly_graph.gfa ${sample_id}_assembly_graph.gfa
    mv flye_output/assembly_info.txt ${sample_id}_assembly_info.txt
    mv flye_output/flye.log ${sample_id}_flye.log
    tail -n 9 ${sample_id}_flye.log | head -n 8 >  ${sample_id}_assembly_stats.txt
    """
}

process canu {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 12

    input:
    tuple val(sample_id), path(reads), val(genome_size)
    output:
    tuple val(sample_id), path("${sample_id}"), emit: canu_outdir
    tuple val(sample_id), path("${sample_id}/${sample_id}.contigs.fasta"), emit: contigs
    
    script:
    """
    hostname > hostname
    canu -p ${sample_id} -d ${sample_id} genomeSize=$genome_size maxThreads=${task.cpus} ${params.canu_options} -nanopore ${reads}
    """
}


process raven {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 8

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_raven_contigs.fasta")
    
    script:
    """
    hostname > hostname
    raven -t ${task.cpus} ${reads} > ${sample_id}_raven_contigs.fasta
    """
}

process miniasm {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag{ sample_id }
    
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_assembly_miniasm.fasta"), emit: contigs
    path("${sample_id}_graph_miniasm.gfa")
    
    script:
    """
    minimap2 -x ava-ont -t ${task.cpus} ${reads} ${reads} > ovlp.paf
    miniasm -f ${reads} ovlp.paf > ${sample_id}_graph_miniasm.gfa
    awk '/^S/{print ">"\$2"\\n"\$3}' ${sample_id}_graph_miniasm.gfa | fold > ${sample_id}_assembly_miniasm.fasta
    """
}


process trycycler_cluster {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy" 
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 8

    input:
    tuple val(sample_id), path(reads), path(contigs)
    output:
    tuple val(sample_id), path("${sample_id}")
    
    script:
    """
    trycycler cluster --assemblies *.fasta --reads ${reads} --out_dir ${sample_id}
    """
}


process unicycler {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"  
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 12

    input:
    tuple val(sample_id), path(short_reads), path(long_reads)
    
    output:
    tuple path("${sample_id}_unicycler.gfa"), path("${sample_id}_unicycler.log")
    tuple val(sample_id), path("${sample_id}_unicycler.fasta"), emit: contigs

    script:
    """
    unicycler -l ${long_reads} -1 ${short_reads[0]} -2 ${short_reads[1]} -o output -t ${task.cpus} --mode ${params.unicycler_mode}
    mv output/assembly.fasta "${sample_id}_unicycler.fasta"
    mv output/assembly.gfa "${sample_id}_unicycler.gfa"
    mv output/unicycler.log "${sample_id}_unicycler.log"
    """
}

process shovill {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: {filename -> if (filename.endsWith(".fa")){ 
                 null
            }
            else {"${filename}"}
            }

    tag {sample_id}
    
    errorStrategy 'retry'
    
    maxRetries 3

    cpus 12

    memory {64.GB + task.attempt * 8.GB}

    input:
    tuple val(sample_id), path(fastq), val(genome_size)
    
    output:
    path("${sample_id}")
    tuple val(sample_id), path("${sample_id}/contigs.fa"), emit: contigs

    script:
    ram = task.memory.toString().replaceAll(/ GB|G/,'')
    """
    shovill --outdir ${sample_id} \
    --R1 ${fastq[0]} --R2 ${fastq[1]} \
    --cpus $task.cpus \
    --assembler $params.assembler \
    --gsize $genome_size \
    --tmpdir ./tmpdir \
    --minlen 300 
    """
}


process skesa {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"  
    
    tag {sample_id}
    
    errorStrategy 'ignore'

    cpus 8

    input:
    tuple val(sample_id), path(short_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_skesa.fasta"), emit: contigs

    script:
    """
    skesa --cores $task.cpus --use_paired_ends --contigs_out ${sample_id}_skesa.fasta --fastq ${short_reads[0]} ${short_reads[1]}
    """
}