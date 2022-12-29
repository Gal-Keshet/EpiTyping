#!/usr/bin/env nextflow

include { EPITYPING } from "$params.module_dir/main_workflow.nf"

log.info """\
         E P I T Y P I N G   P I P E L I N E    
         ===================================
         run memory             : ${params.mem}
         threads                : ${params.cpu}
         outdir                 : ${params.outdir}

         RNA-seq parameters:
         RNA fastq dir          : ${params.fastq_folder}
         multiple samples       : ${params.multiple_samples}
         RNA single end         : ${params.single}
         RNA strand             : ${params.strandness}
         filter mouse           : ${params.mouse_feeders}

         DNA-seq parameters (taken into account only if DNA-seq is provided):
         DNA sample             : ${params.dna_reference}
         filter mouse from DNA  : ${params.mouse_feeders_dna}
         DNA single_end         : ${params.dna_single}
         DNA dir                : ${params.dna_fastq_folder}
         ==================================
         """
         .stripIndent()


workflow {
    EPITYPING()
}
