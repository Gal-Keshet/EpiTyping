#!/usr/bin/env nextflow

include { TRIM_FASTQ;\
 ALIGNMENT_AND_MOUSE_REMOVAL;\
 ALIGNMENT;\
 XENOFILTER;\
 FEATURE_COUNT;\
 ADD_READ_GROUPS;\
 MARK_DUPLICATES;\
 SPLIT_N_CIGAR;\
 BASE_RECALIBRATION;\
 HAPLOTYPE_CALLER;\
 VARIANT_FILTRATION;\
 ADD_GENE_NAMES;\
 ADD_dbSNP_ID;\
 KEEP_ANNOTATED_SNPS }\
 from "$params.module_dir/call_vars.nf"

include { COMBINE_LOI_VCF; OUTPUT_LOI_TABLE; LOI_MATRIX } from "$params.module_dir/loi.nf"

include { COMBINE_X_VCF; OUTPUT_XCI_TABLE; XCI_MATRIX } from "$params.module_dir/xci.nf"

include { DNA } from "$params.module_dir/dna_workflow.nf"

include { LOI_MATRIX_DNA_INTEGRATION } from "$params.module_dir/dna_integration.nf"


workflow EPITYPING {
    
    if( params.single == false ) {
      Channel
      .fromFilePairs( "$params.fastq_folder/*_{1,2}.fastq.gz", checkIfExists:true )
      .set { read_files }
    }
    else {
      Channel
      .fromFilePairs( "$params.fastq_folder/*.fastq.gz", size: 1, checkIfExists:true )
      .set { read_files }
    }

    read_files.view()

    trimmed_fastq_ch = TRIM_FASTQ(read_files, params.common_adapters)

    if( params.mouse_feeders == true ) {
      bam_ch = ALIGNMENT_AND_MOUSE_REMOVAL(trimmed_fastq_ch, params.index, params.mouse_index)
      no_mouse_reads_ch = XENOFILTER(bam_ch)
    }
    else {
      no_mouse_reads_ch = ALIGNMENT(trimmed_fastq_ch, params.index)
    }

    counts_ch = FEATURE_COUNT(no_mouse_reads_ch, params.gtf)

    bam_with_rg_ch = ADD_READ_GROUPS(no_mouse_reads_ch)

    bam_with_marked_dup_ch = MARK_DUPLICATES(bam_with_rg_ch)

    bam_cigar_split_ch = SPLIT_N_CIGAR(bam_with_marked_dup_ch, params.ref, "${params.ref}.fai", params.ref_dict)

    recalibrated_ch = BASE_RECALIBRATION(bam_cigar_split_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals,\
     "${params.exon_intervals}.tbi", params.dbsnp, "${params.dbsnp}.tbi")

    called_haplotypes_ch = HAPLOTYPE_CALLER(recalibrated_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals, "${params.exon_intervals}.tbi")

    vcf_filterd_ch = VARIANT_FILTRATION(called_haplotypes_ch, params.ref, "${params.ref}.fai", params.ref_dict, params.exon_intervals, "${params.exon_intervals}.tbi")

    vcf_with_names_ch = ADD_GENE_NAMES(vcf_filterd_ch, params.exon_intervals, "${params.exon_intervals}.tbi", params.headers)

    vcf_with_dbsnp_ch = ADD_dbSNP_ID(vcf_with_names_ch, params.dbsnp, "${params.dbsnp}.tbi")

    vcf_with_dbsnp_filtered_ch = KEEP_ANNOTATED_SNPS(vcf_with_dbsnp_ch)

    merged_vcf_paths_ch = vcf_with_dbsnp_filtered_ch
    .map { it[0] }
    .collect()

    merged_vcf_index_paths_ch = vcf_with_dbsnp_filtered_ch
    .map { it[1] }
    .collect()

    counts_files_ch = counts_ch.collect()

    if( merged_vcf_paths_ch.map { it -> it.size() > 1} == true ) {
        
      merged_loi_vcf_ch = COMBINE_LOI_VCF(merged_vcf_paths_ch, merged_vcf_index_paths_ch, params.imprinted_intervals)

      var_loi_table_ch = OUTPUT_LOI_TABLE(merged_loi_vcf_ch, params.imprinted_intervals)

      merged_vcf_x_ch = COMBINE_X_VCF(merged_vcf_paths_ch, merged_vcf_index_paths_ch)

      var_x_table_ch = OUTPUT_XCI_TABLE(merged_vcf_x_ch)
    }

    else {
      var_loi_table_ch = OUTPUT_LOI_TABLE(vcf_with_dbsnp_filtered_ch, params.imprinted_intervals)
      var_x_table_ch = OUTPUT_XCI_TABLE(merged_vcf_paths_ch)
    }

    xci_out_ch = XCI_MATRIX(var_x_table_ch, counts_ch.collect())

    if(params.dna_reference != "" ) {
      dna_table_ch = DNA()
      loi_matrix_dna_integrated_ch = LOI_MATRIX_DNA_INTEGRATION(dna_table_ch, var_loi_table_ch, counts_files_ch, params.imprinted_genes)
    }
    else {
      loi_matrix_ch = LOI_MATRIX(var_loi_table_ch, counts_files_ch, params.imprinted_genes)
    }
}

