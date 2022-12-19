#!/usr/bin/env nextflow

params.keepInter = false

params.module_dir = "$projectDir/modules"

params.multiple_samples = false

params.single = false

params.strandness = 0

params.mouse_feeders = true

params.fastq_folder = "$PWD/fastq"

params.common_adapters = "$projectDir/genome_files/CommonAdapters.fa"

params.index = "$projectDir/genome_files/star_index"

params.mouse_index = "$projectDir/genome_files/star_mouse_index"

params.gtf = "$projectDir/genome_files/star_ref/gencode.v42.primary_assembly.annotation.gtf"

params.ref = "$projectDir/genome_files/star_ref/GRCh38.primary_assembly.genome.fa"

params.ref_dict = "$projectDir/genome_files/star_ref/GRCh38.primary_assembly.genome.dict"

params.dbsnp = "$projectDir/genome_files/dbSNP155Renamed.vcf.gz"

params.exon_intervals = "$projectDir/genome_files/GRCh38_exome.bed.gz"

params.imprinted_intervals = "$projectDir/genome_files/imprinted_genes_sorted.bed"

params.imprinted_genes = "$projectDir/genome_files/imprinted_genes_location.csv"

params.headers = "$projectDir/genome_files/headers_for_annotation.txt"

params.bwa_index = "$projectDir/genome_files/bwa_ref/*"

params.dna_fastq_folder = "$PWD/dna_fastq"

params.dna_reference = ""

params.dna_single = false

params.outdir = "$PWD/output"

include { EPITYPING } from "$params.module_dir/main_workflow.nf"

log.info """\
         E P I T Y P I N G   P I P E L I N E    
         ===================================
         run memory         : ${params.mem}
         threads            : ${params.cpu}
         RNA fastq dir      : ${params.fastq_folder}
         multiple samples   : ${params.multiple_samples}
         RNA single end     : ${params.single}
         RNA strand         : ${params.strandness}
         filter mouse       : ${params.mouse_feeders}
         DNA sample         : ${params.dna_reference}
         DNA single_end     : ${params.dna_single}
         DNA dir            : ${params.dna_fastq_folder}
         outdir             : ${params.outdir}
         ==================================
         """
         .stripIndent()


workflow {
    EPITYPING()
}
