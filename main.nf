#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

MODULES = projectDir + '/modules/'

include {flye} from MODULES + '/assembler.nf'
include {busco} from MODULES + '/qc.nf'
include {quast} from MODULES + '/qc.nf'
include {medaka} from MODULES + '/polishing.nf'
include {nanoq} from MODULES + '/qc.nf'
include {mob_recon} from MODULES + '/annotation.nf'
include {staramr} from MODULES + '/annotation.nf'

//Reads
ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
                  .map {it -> tuple(it.getSimpleName(), file(it))}

//Genome size for Flye
ch_genome_size = Channel.of(params.genome_size)

log.info 'C A M P Y PIPELINE'
log.info "Input: ${params.reads}"
log.info "Out dir: ${params.outdir}"


workflow {
    // Filtering reads with cutoff length and qscore
    nanoq(ch_reads)
    // De-novo assembly with Flye
    flye(nanoq.out.combine(ch_genome_size))
    // Polishing using medaka
    medaka(nanoq.out.join(flye.out.contigs))
    //Evaluate contigs completeness
    if (params.busco) {
        busco(medaka.out)
    }
    //Evaluate contigs contiguity
    quast(medaka.out.map{it -> it[1]}.collect())
    // Reconstruct chromosome and plasmid from the assembly
    if (params.plasmid) {
        mob_recon(medaka.out)
        staramr_ch = mob_recon.out.contigs.flatten()
    } else {
        staramr_ch = medaka.out.map{it -> it[1]}
    }
    
    // Find genes of interest, AMR, MLST
    staramr(staramr_ch.collect())
}