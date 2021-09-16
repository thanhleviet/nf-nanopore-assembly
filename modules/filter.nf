params.filter_reads = false

process filter_campy {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 8

    input:
    tuple val(sample_id), path(reads), path(kraken2_output)
   
    output:
    tuple path("${sample_id}_campy.tsv"), path("*.{fq.gz,fastq.gz,fasta}") optional true

    shell:
    if (params.filter_reads) {
        '''
        cat !{kraken2_output} | grep -i "campy" | cut -f 1 > campy
        if [[ -s campy ]]; then 
            filterbyname.sh in=!{reads[0]} in2=!{reads[1]} names=campy out="!{sample_id}_campy_1.fq.gz" out2="!{sample_id}_campy_2.fq.gz" include=t qin=33
            mv campy !{sample_id}_campy.tsv
        fi
        '''
    } else {
        '''
        cat !{kraken2_output} | grep -i "campy" | cut -f 2 > campy
        if [[ -s campy ]]; then 
            filterbyname.sh in=!{reads} names=campy out="!{sample_id}_campy.fastq.gz" include=t qin=33
            mv campy !{sample_id}_campy.tsv
        fi
        '''
    }
    
}