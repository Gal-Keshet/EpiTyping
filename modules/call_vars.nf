process TRIM_FASTQ {
    debug true

    publishDir "$params.outdir/trimmed_fastq", mode:'copy', enabled: "$params.keepInter" == true 

    input:
    tuple val(sampleId), path(fastq)
    path adapters

    output:
    tuple val("${sampleId}"), path("${sampleId}*")

    script:
    def single = fastq instanceof Path
    if( !single ) {
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE\
        -threads $task.cpus\
        ${fastq[0]} ${fastq[1]}\
        ${sampleId}_trimmed_1.fastq.gz ${sampleId}_trimmed_1_unpaired.fastq.gz\
        ${sampleId}_trimmed_2.fastq.gz ${sampleId}_trimmed_2_unpaired.fastq.gz\
        ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    }    
    else {
        """
        java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar SE\
        -threads $task.cpus\
        $fastq\
        ${sampleId}_trimmed.fastq.gz\
        ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    }
}

process ALIGNMENT_AND_MOUSE_REMOVAL {
    debug true

    publishDir "$params.outdir/star_out/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(fastq)
    path index_human
    path index_mouse

    output:
    tuple val("$sampleId"), path("${sampleId}_Aligned.sortedByCoord.out.bam"), path("${sampleId}_mouse_Aligned.sortedByCoord.out.bam")

    script:
    def single = fastq instanceof Path
    if( !single ) {
        """
        STAR --runThreadN $task.cpus\
        --genomeDir $index_human\
        --readFilesIn ${fastq[0]} ${fastq[2]}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_

        STAR --runThreadN $task.cpus\
        --genomeDir $index_mouse\
        --readFilesIn ${fastq[0]} ${fastq[2]}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_mouse_
        """
    }
    else {
        """
        STAR --runThreadN $task.cpus\
        --genomeDir $index_human\
        --readFilesIn ${fastq}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_

        STAR --runThreadN $task.cpus\
        --genomeDir $index_mouse\
        --readFilesIn ${fastq}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_mouse_
        """
    }
}

process ALIGNMENT {
    debug true

    publishDir "$params.outdir/star_out/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(fastq)
    path index_human

    output:
    tuple val("$sampleId"), path("${sampleId}_Aligned.sortedByCoord.out.bam")

    script:
    def single = fastq instanceof Path
    if( !single ) {
        """
        STAR --runThreadN $task.cpus\
        --genomeDir $index_human\
        --readFilesIn ${fastq[0]} ${fastq[2]}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_
        """
    }
    else {
        """
        STAR --runThreadN $task.cpus\
        --genomeDir $index_human\
        --readFilesIn ${fastq}\
        --outSAMtype BAM SortedByCoordinate\
        --twopassMode Basic\
        --outSAMattributes NM\
        --quantMode GeneCounts\
        --limitBAMsortRAM 35145460689\
        --readFilesCommand zcat\
        --outFileNamePrefix ${sampleId}_
        """
    }
}

process XENOFILTER {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true
  
    input:
    tuple val(sampleId), path(human_bam_file), path(mouse_bam_file)

    output:
    tuple val("$sampleId"), path("Filtered_bams/${sampleId}_Aligned.sortedByCoord.out_Filtered.bam")
    
    if ( file("Filtered_bams") ) {
      file("Filtered_bams").delete()
    } 

    script:
    """
    #!/usr/bin/env Rscript
    library("XenofilteR")
    sample.list <- data.frame(human = "$human_bam_file", mouse="$mouse_bam_file")
    bp.param <- SnowParam(workers = "$task.cpus", type = "SOCK")
    XenofilteR(sample.list, destination.folder = "./", bp.param = bp.param)
    """
}

process FEATURE_COUNT {
    debug true

    publishDir "$params.outdir/feature_count/$sampleId", mode:'copy', enabled: "$params.keepInter" == true
 
    input:
    tuple val(sampleId), path(bam_file)
    path gtf

    output:
    path("${sampleId}_edited_counts.txt")
    
    script:   
    if (params.single == true) {
      """
      featureCounts\
      --extraAttributes "gene_name"\
      -T $task.cpus\
      -a $gtf\
      -s $params.strandness\
      -o ${sampleId}_counts.txt\
      $bam_file

      tail -n +3 ${sampleId}_counts.txt | awk 'BEGIN { print "Gene ID\tChr\tLength\tsymbol\tcounts" } { print \$1 "\t" \$2 "\t" \$6 "\t" \$7 "\t" \$8 }' > ${sampleId}_edited_counts.txt
      """
    }
    else {
      """
      featureCounts\
      --extraAttributes "gene_name"\
      -T $task.cpus\
      -a $gtf\
      -s $params.strandness\
      -o ${sampleId}_counts.txt\
      -p\
       --countReadPairs\
      $bam_file

      tail -n +3 ${sampleId}_counts.txt | awk 'BEGIN { print "Gene ID\tChr\tLength\tsymbol\tcounts" } { print \$1 "\t" \$2 "\t" \$6 "\t" \$7 "\t" \$8 }' > ${sampleId}_edited_counts.txt
      """
    }
}

process ADD_READ_GROUPS {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val("$sampleId"), path("${sampleId}_Aligned.out_rg.bam")

    script:
    """
    java -jar /picard.jar AddOrReplaceReadGroups\
    I=$bam_file\
    O=${sampleId}_Aligned.out_rg.bam\
    SO=coordinate RGID=rnasq RGLB=lb RGPL=illumina RGPU=pu RGSM=${sampleId}
    """
}

process MARK_DUPLICATES {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val("$sampleId"), path("${sampleId}_dedupped.bam")

    script:
    """
    java -jar /picard.jar MarkDuplicates\
    I=$bam_file\
    O=${sampleId}_dedupped.bam\
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
    """
}

process SPLIT_N_CIGAR {
    debug true
    label 'var_call'

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam_file)
    path genome
    path index
    path dict

    output:
    tuple val("$sampleId"), path("${sampleId}_split.bam")

    script:
    """
    gatk SplitNCigarReads\
    -R $genome\
    -I $bam_file\
    -O ${sampleId}_split.bam
    """
}

process BASE_RECALIBRATION {
    debug true
    label 'var_call'

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam_file)
    path genome
    path index
    path dict
    path annotation_intervals
    path annotation_intervals_index
    path dbsnp
    path dbsnp_index

    output:
    tuple val("$sampleId"), path("${sampleId}_recalibrated.bam")

    script:
    """
    gatk BaseRecalibrator\
    -R $genome\
    -L $annotation_intervals\
    --known-sites $dbsnp\
    -ip 100\
    -I $bam_file\
    -O ${sampleId}_recal_data.table

    gatk ApplyBQSR\
    -R $genome\
    -L $annotation_intervals\
    --bqsr-recal-file ${sampleId}_recal_data.table\
    -ip 100\
    -I $bam_file\
    -O ${sampleId}_recalibrated.bam
    """
}

process HAPLOTYPE_CALLER {
    debug true
    label 'var_call'

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam_file)
    path genome
    path index
    path dict
    path annotation_intervals
    path annotation_intervals_index

    output:
    tuple val("$sampleId"), path("${sampleId}_output.vcf*")

    script:
    """
    gatk --java-options "-Xmx8G" HaplotypeCaller\
    -R $genome\
    -L $annotation_intervals\
    --dont-use-soft-clipped-bases\
    -stand-call-conf 20.0\
    -I $bam_file\
    -O ${sampleId}_output.vcf

    bgzip -f ${sampleId}_output.vcf
    tabix -f ${sampleId}_output.vcf.gz
    """
}

process VARIANT_FILTRATION {
    debug true
    label 'var_call'

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)
    path genome
    path index
    path dict
    path annotation_intervals
    path annotation_intervals_index

    output:
    tuple val("$sampleId"), path("${sampleId}_filtered.vcf*")

    script:
    """
    gatk VariantFiltration\
    -R $genome\
    -L $annotation_intervals\
    -window 35 -cluster 3\
    --filter-name FS --filter-expression "FS > 30.0"\
    --filter-name QD --filter-expression "QD < 2.0"\
    -V ${vcf_file[0]}\
    -O ${sampleId}_filtered.vcf
    
    bgzip -f ${sampleId}_filtered.vcf
    tabix -f ${sampleId}_filtered.vcf.gz
    """
}

process ADD_GENE_NAMES {
    debug true

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)
    path annotation_intervals
    path annotation_intervals_index
    path headers

    output:
    tuple val("$sampleId"), path("${sampleId}_names.vcf*")

    script:
    """
    bcftools annotate --threads $task.cpus\
    --annotations $annotation_intervals\
    -h $headers\
    -c CHROM,FROM,TO,Gene\
    --output-type z\
    --output ${sampleId}_names.vcf.gz\
    ${vcf_file[0]}

    tabix -f ${sampleId}_names.vcf.gz
    """
}

process ADD_dbSNP_ID {
    debug true

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)
    path dbSNP
    path dbSNP_index

    output:
    tuple val("$sampleId"), path("${sampleId}_dbSNP.vcf*")

    script:
    """
    bcftools annotate --threads $task.cpus\
    --annotations $dbSNP\
    -c ID,INFO/RS,INFO/COMMON\
    --output-type z\
    --output ${sampleId}_dbSNP.vcf.gz\
    ${vcf_file[0]}

    tabix -f ${sampleId}_dbSNP.vcf.gz
    """
}

process KEEP_ANNOTATED_SNPS {
    debug true

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)

    output:
    tuple path("${sampleId}_dbSNP_filtered.vcf.gz"), path("${sampleId}_dbSNP_filtered.vcf.gz.tbi")

    script:
    """
    bcftools view --threads $task.cpus\
    -i 'INFO/RS!="." & FILTER="PASS" & FMT/DP > 9'\
    --output-type z\
    --output ${sampleId}_dbSNP_filtered.vcf.gz\
    ${vcf_file[0]}

    tabix -f ${sampleId}_dbSNP_filtered.vcf.gz
    """
}
