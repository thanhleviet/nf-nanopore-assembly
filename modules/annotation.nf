params.prokka_options = ""

process prokka {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    // saveAs: {filename -> if (file(filename).isDirectory()){ 
    //             "${filename}" 
    //             }
    //          else null
    //          }
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 8

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}")
    path "${sample_id}/${sample_id}.gff" , emit: roary
    
    script:
    """
    prokka --kingdom Bacteria ${params.prokka_options} --outdir ${sample_id} --prefix ${sample_id} ${contigs}
    """
}


process roary {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 16

    input:
    path(gffs)
    output:
    path("roary_results")

    script:
    """
    roary -p ${task.cpus} -e --mafft -f roary_results ${gffs}
    """
}

process mlst {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 1

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}_mlst.tsv")

    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.tsv
    """
}


process staramr {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 8

    input:
    path(contigs)
    output:
    path("staramr_output")

    script:
    """
    staramr search --output-dir staramr *.fa*
    """
}



process amrfinderplus {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    errorStrategy 'ignore'
    
    cpus 2

    input:
    tuple val(sample_id), path(contigs)
    output:
    path("${sample_id}.*")

    script:
    """
    amrfinder -t ${task.cpus} \
    --organism ${params.amrfinderplus_species} \
    --nucleotide ${contigs} \
    --report_common \
    -o ${sample_id}.armfinder \
    --name ${sample_id} \
    --nucleotide_output ${sample_id}.amrgenes.fna \
    --plus
    """
}

process mob_recon {
  publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy",
    saveAs: { filename -> if (filename.endsWith(".tsv")){ 
                "report/${filename}"  } 
            else if (file(filename).isDirectory())
                {"${filename}"} 
            else { 
                "${filename}"
            }
    }
  
  container 'kbessonov/mob_suite:3.0.3'

  tag {sample_id}

  cpus 16 

  input:
  
  tuple val(sample_id), path(contigs)
  
  output:
  path("${sample_id}"), emit: outdir
  path("${sample_id}/*.fasta"), emit: contigs
  path("*.tsv"), optional: true
  
  script:
  mod_db = params.mod_db ? "-d ${params.mob_db}" : ""
  ref_db = params.ref_db ? "-g ${params.ref_db}" : ""
  """
    hostname > hostname
    mob_recon \
    --num_threads ${task.cpus} \
    --run_typer \
    -s ${sample_id} \
    ${mod_db} \
    ${ref_db} \
    --infile ${contigs} \
    --outdir ${sample_id}
    cp ${sample_id}/contig_report.txt ${sample_id}_contig_report.tsv 2>/dev/null || :
    cp ${sample_id}/mobtyper_results.txt ${sample_id}_mobtyper_results.tsv 2>/dev/null || :
  """
}
