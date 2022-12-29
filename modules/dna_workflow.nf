#!/usr/bin/env nextflow

include { TRIM_DNA_FASTQ;\
 ALIGN_DNA;\
 ALIGN_DNA_MOUSE;\
 SORT_MOUSE_BAM;\
 CUT_BAM;\
 ADD_DNA_READ_GROUPS;\
 SORT_BAM;\
 DNA_XENOFILTER;\
 MARK_DNA_DUPLICATES;\
 DNA_BASE_RECALIBRATION;\
 DNA_HAPLOTYPE_CALLER;\
 DNA_VARIANT_FILTRATION;\
 CNN_SCORE_VARIANTS;\
 FILTER_VARIANT_TRANCHES;\
 ADD_dbSNP_ID_DNA;\
 ADD_DNA_GENE_NAMES;\
 FINAL_DNA_VCF;\
 OUTPUT_DNA_TABLE;\
 LOI_MATRIX_DNA_INTEGRATION } from "$params.module_dir/dna_integration.nf"

workflow DNA {

    main:

    if( params.dna_single == false ) {
      Channel
      .fromFilePairs( "$params.dna_fastq_folder/${params.dna_reference}_{1,2}.fastq.gz", checkIfExists:true )
      .set { dna_read_files }
    }
    else {
      Channel
      .fromFilePairs( "$params.dna_fastq_folder/${params.dna_reference}.fastq.gz", size: 1, checkIfExists:true )
      .set { dna_read_files }
    }

    trimmed_dna_fastq_ch = TRIM_DNA_FASTQ(dna_read_files, params.common_adapters)

    Channel
     .fromPath( params.bwa_index )
     .set { bwa_index_ch }

    Channel
     .fromPath( params.bwa_mouse_index )
     .set { bwa_mouse_index_ch }

    if( params.mouse_feeders_dna == true ) {

      dna_bam_ch = ALIGN_DNA(trimmed_dna_fastq_ch, params.ref, bwa_index_ch.collect())

      dna_mouse_bam_ch = ALIGN_DNA_MOUSE(trimmed_dna_fastq_ch, params.mouse_ref, bwa_mouse_index_ch.collect())
      
      dna_mouse_sorted_ch = SORT_MOUSE_BAM(dna_mouse_bam_ch)

      dna_bam_cut_ch = CUT_BAM(dna_bam_ch, params.exon_intervals)

      dna_bam_rg_ch = ADD_DNA_READ_GROUPS(dna_bam_cut_ch)

      dna_sorted_bam_ch = SORT_BAM(dna_bam_rg_ch)

      dna_no_mouse = DNA_XENOFILTER(dna_sorted_bam_ch, dna_mouse_sorted_ch)

      dna_marked_dup_ch = MARK_DNA_DUPLICATES(dna_no_mouse)
    }
    else {

      dna_bam_ch = ALIGN_DNA(trimmed_dna_fastq_ch, params.ref, bwa_index_ch.collect())

      dna_bam_cut_ch = CUT_BAM(dna_bam_ch, params.exon_intervals)

      dna_bam_rg_ch = ADD_DNA_READ_GROUPS(dna_bam_cut_ch)

      dna_sorted_bam_ch = SORT_BAM(dna_bam_rg_ch)

      dna_marked_dup_ch = MARK_DNA_DUPLICATES(dna_sorted_bam_ch)
    }

    dna_recalibrated_ch = DNA_BASE_RECALIBRATION(dna_marked_dup_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals,\
     "${params.exon_intervals}.tbi", params.dbsnp, "${params.dbsnp}.tbi")

    dna_haplotype_ch = DNA_HAPLOTYPE_CALLER(dna_recalibrated_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals, "${params.exon_intervals}.tbi")

    dna_filtered_ch = DNA_VARIANT_FILTRATION(dna_haplotype_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals, "${params.exon_intervals}.tbi")

    dna_cnn_scores_ch = CNN_SCORE_VARIANTS(dna_filtered_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals, "${params.exon_intervals}.tbi")

    dna_tranches_filtered_ch = FILTER_VARIANT_TRANCHES(dna_cnn_scores_ch, params.exon_intervals, "${params.exon_intervals}.tbi",\
     params.dbsnp, "${params.dbsnp}.tbi")

    dna_dbsnp_ch = ADD_dbSNP_ID_DNA(dna_tranches_filtered_ch, params.dbsnp, "${params.dbsnp}.tbi")

    dna_named_ch = ADD_DNA_GENE_NAMES(dna_dbsnp_ch, params.exon_intervals, "${params.exon_intervals}.tbi", params.headers)

    dna_vcf_ch = FINAL_DNA_VCF(dna_named_ch, params.imprinted_intervals)

    dna_table_ch = OUTPUT_DNA_TABLE(dna_vcf_ch)

    emit:
    dna_table_ch
}
