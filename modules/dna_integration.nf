process TRIM_DNA_FASTQ {
    debug true

    publishDir "$params.outdir/dna_trimmed_fastq", mode:'copy', enabled: "$params.keepInter" == true

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

process ALIGN_DNA {
    time '1d'

    debug true
 
    input:
    tuple val(sampleId), path(fastq)
    path ref    
    path index

    output:
    tuple val("$sampleId"), path("${sampleId}_aligned.bam")  
    
    script:
    def single = fastq instanceof Path
    if( !single ) {
        """
        bwa mem -M -t $task.cpus\
        $ref\
        ${fastq[0]} ${fastq[2]} | samtools view -b -@ $task.cpus -o ${sampleId}_aligned.bam -
        """
    }
    else {
        """
        bwa mem -M -t $task.cpus\
        $ref\
        $fastq | samtools view -b -@ $task.cpus -o ${sampleId}_aligned.bam -
        """
    }
}

process CUT_BAM {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam)
    path intervals

    output:
    tuple val("$sampleId"), path("${sampleId}_cut.bam")

    script:
    """
    samtools view -b -L $intervals -o ${sampleId}_cut.bam -@ $task.cpus $bam
    """
}

process ADD_DNA_READ_GROUPS {
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

process SORT_BAM {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(bam)

    output:
    tuple val("$sampleId"), path("${sampleId}_sorted.bam")

    script:
    """
    java -jar /picard.jar SortSam\
    I=$bam O=${sampleId}_sorted.bam\
    SORT_ORDER=coordinate CREATE_INDEX=true
    """
}

process MARK_DNA_DUPLICATES {
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

process DNA_BASE_RECALIBRATION {
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

process DNA_HAPLOTYPE_CALLER {
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

process DNA_VARIANT_FILTRATION {
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
    gatk --java-options "-Xms3000m" VariantFiltration\
    -R $genome\
    -L $annotation_intervals\
    -V ${vcf_file[0]}\
    -O ${sampleId}_filtered.vcf\
    --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0"\
    --filter-name "HardFiltered"

    bgzip -f ${sampleId}_filtered.vcf
    tabix -f ${sampleId}_filtered.vcf.gz
    """
}

process CNN_SCORE_VARIANTS {
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
    tuple val("$sampleId"), path("${sampleId}_CNN.vcf*")

    script:
    """
    gatk --java-options -Xmx10g CNNScoreVariants\
    -R $genome\
    -L $annotation_intervals\
    -V ${vcf_file[0]}\
    -ip 100\
    -O ${sampleId}_CNN.vcf

    bgzip -f ${sampleId}_CNN.vcf
    tabix -f ${sampleId}_CNN.vcf.gz
    """
}

process FILTER_VARIANT_TRANCHES {
    debug true
    label 'var_call'

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)
    path annotation_intervals
    path annotation_intervals_index
    path dbsnp
    path dbsnp_index

    output:
    tuple val("$sampleId"), path("${sampleId}_tranches_filter.vcf*")

    script:
    """
    gatk --java-options -Xmx6g FilterVariantTranches\
    -L $annotation_intervals\
    -V ${vcf_file[0]}\
    -O ${sampleId}_tranches_filter.vcf\
    -ip 100\
    --resource $dbsnp\
    --info-key CNN_2D\
    --snp-tranche 99.95\
    --indel-tranche 99.4

    bgzip -f ${sampleId}_tranches_filter.vcf
    tabix -f ${sampleId}_tranches_filter.vcf.gz
    """
}

process ADD_dbSNP_ID_DNA {
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

process ADD_DNA_GENE_NAMES {
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


process FINAL_DNA_VCF {
    debug true

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(vcf_file)
    path imprinted_intervals

    output:
    path "dna.vcf.gz"

    script:
    """
    bcftools view --threads $task.cpus\
    -i 'INFO/RS!="." & FILTER="PASS" & FMT/DP > 9'\
    --regions-file $imprinted_intervals\
    --output-type z\
    --output dna.vcf.gz\
    ${vcf_file[0]}
    """
}   

process OUTPUT_DNA_TABLE {
    debug true

    publishDir "$params.outdir/vcf/loi", mode:'copy', enabled: "$params.keepInter" == true

    input:
    path vcf_file

    output:
    path "dna_final_table.txt"

    script:
    """
    bcftools query\
    -H -f '%CHROM\t%POS\t%ID\t%INFO/Gene\t%INFO/COMMON\t%REF\t%ALT[\t%GT:%AD{0}:%AD{1}:%DP]\n'\
    --output dna_final_table.txt\
    $vcf_file
    """
}

process LOI_MATRIX_DNA_INTEGRATION {
    debug true
    publishDir "$params.outdir/results/loi", mode:'copy'
    input:
    path dna_table
    path var_table
    path count_files
    path imprinted_genes_file

    output:
    path "LOI_matrix_with_dna_integration.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    
    #upload imprinted genes
    gene_locations <- read.csv("$imprinted_genes_file")
    
    #read dna variants
    dna.file <- read.delim("$dna_table", stringsAsFactors = F)
    
    #process dataframe and keep variants that are annotated to a gene's exon (Gene),
    #that have a dbSNP ID (RS) and are one base-pair long
    names(dna.file) <- str_replace(names(dna.file), 'X\\\\.+[0-9]+\\\\.', '')
    names(dna.file) <- str_replace(names(dna.file), '\\\\.GT.*', ':genotype;allele1;allele2;sum')
    
    dna.file <- dna.file[(grepl('[A-Z]+',dna.file\$Gene)) & (grepl('^[ACGT]{,2}\$', dna.file\$REF)) & (grepl('^[ACGT]{,2}\$', dna.file\$ALT)), -5]
    
    rownames(dna.file) <- NULL
    
    dna.file <- merge(gene_locations, dna.file, by = 'Gene')
    
    #replace missing snps with nan for each sample, extract information
    #about present values and calulate biallelic ratio (minor allele / major allele).
    #than only keep snps that have a biallelic ratio bigger than 0.2, coverage of at-least 10 reads
    #and are eterozygous (genotype=0/1) 
    return_allelic_freq <- function(Sample) {
      
      has_info <- grep('\\\\./\\\\.:\\\\.:\\\\.:\\\\.', Sample, invert = TRUE)
      no_info <- grep('\\\\./\\\\.:\\\\.:\\\\.:\\\\.', Sample)
      Sample[no_info] <- NaN
      
      snp_info <- data.frame(str_split(Sample[has_info], ':', simplify = T))
      snp_info[, c(2:4)] <- apply(snp_info[, c(2:4)], 2, as.numeric)
      
      allelic_freq <- apply(snp_info[,c(2:3)], 1, function(x){return(min(x) / max(x))})
      snp_info <- cbind(snp_info, allelic_freq)
      
      heterozygot_pos <- grep('^0/1\$', snp_info[,1])
      homozigot_pos <- grep('^0/1\$', snp_info[,1], invert = T)
      snp_info[homozigot_pos, 5] <- 0
      
      snp_info[snp_info[,5] < 0.2, 5] <- 0
      snp_info[snp_info[, 4] < 10, 5] <- 0 
      
      Sample[has_info] <- snp_info[, 5]
      Sample <- as.numeric(Sample)
      Sample[Sample == 0] <- NaN
      
      return(Sample)
    }
    
    dna_allele_table <- cbind.data.frame(dna.file[c(1:7)], apply(dna.file[-c(1:7)], 2, return_allelic_freq))
    
    heterozygote_pos <- dna_allele_table[!is.na(dna_allele_table[8]),-8]
    ###################
    
    #read variant call file
    variants.file <- read.delim("$var_table", stringsAsFactors = F)
    
    #process dataframe and keep variants that are annotated to a gene's exon (Gene),
    #that have a dbSNP ID (RS) and are one base-pair long
    names(variants.file) <- str_replace(names(variants.file), 'X\\\\.+[0-9]+\\\\.', '')
    names(variants.file) <- str_replace(names(variants.file), '\\\\.GT.*', ':genotype;allele1;allele2;sum')
    
    variants.file <- variants.file[(grepl('[A-Z]+',variants.file\$Gene)) & (grepl('^[ACGT]{,2}\$', variants.file\$REF)) & (grepl('^[ACGT]{,2}\$', variants.file\$ALT)), -5]
    
    rownames(variants.file) <- NULL
    
    merged_table <- merge(heterozygote_pos, variants.file, all.x = T)
    
    merged_table[is.na(merged_table)] <- './.:.:.:.'
    
    merged_table <- merged_table[,-c(2:6)]
    
    merged_table_frequencies <- cbind.data.frame(merged_table[c(1:2)], apply(merged_table[-c(1:2)], 2, return_allelic_freq))
    
    colnames(merged_table_frequencies) <- str_replace(colnames(merged_table_frequencies), ':.*', '')
    
    merged_table_frequencies <- as.data.frame(merged_table_frequencies %>%
      group_by(Gene, Location) %>%
      summarise_all(.,mean, na.rm = TRUE))
    
    row.names(merged_table_frequencies) <- merged_table_frequencies\$Gene
    
    merged_table_frequencies <- merged_table_frequencies[,-1]
    
    merged_table_frequencies <- merged_table_frequencies[order(merged_table_frequencies\$Location),]
    
    merged_table_frequencies[-1][is.na(merged_table_frequencies[-1])] <- 0
    
    merged_table_frequencies[-1][merged_table_frequencies[-1] > 0] <- 1
    
    ###count files
    files_vector <- unlist(str_split("$count_files", pattern = ' '))    
    #merging counts
    for (file in files_vector){
      
      # if the merged dataset doesn't exist, create it
      if (!exists("dataset")){
        dataset <- read.table(file, header=TRUE, sep="\t")
        colnames(dataset) <- c('Gene.ID', 'Chr', 'Length', 'gene_name', file)
        dataset\$Chr <- str_replace(dataset\$Chr, ';.*', '')
      }
      
      # if the merged dataset does exist, append to it
      else if (exists("dataset")){
        temp_dataset <-read.table(file, header=TRUE, sep="\t")
        colnames(temp_dataset) <- c('Gene.ID', 'Chr', 'Length', 'gene_name', file)
        temp_dataset\$Chr <- str_replace(temp_dataset\$Chr, ';.*', '')
        dataset<-merge(dataset, temp_dataset, by=c('Gene.ID', 'Chr', 'gene_name', 'Length'))
        rm(temp_dataset)
      }
    }
    
    colnames(dataset) <- str_replace(colnames(dataset), '.*/', '')
    colnames(dataset) <- str_replace(colnames(dataset), '_edited.*', '')
    
    #######end of check######
    counts <- dataset[-c(1:4)]
    rpk= counts / dataset\$Length #first normalize for gene length
    tot_RPK_per_samp=apply(rpk,2,sum)/10^6 #Sum the normalized values
    tpm=t(t(rpk)/tot_RPK_per_samp) # Normalize to library size
    tpm <- cbind(dataset[1:4], tpm)
    
    #subsetting tpm for imprinted genes & changing order of biallelic score matrix to tpm matrix order
    tpm <- tpm[tpm\$gene_name %in% row.names(merged_table_frequencies),]
    
    row.names(tpm) <- tpm\$gene_name
    
    tpm <- tpm[row.names(merged_table_frequencies),]
    tpm <- tpm[colnames(merged_table_frequencies)[-1]]
    
    #keeping allele scores only for expressed genes
    merged_table_frequencies[-1][tpm < 1] <- NaN
    
    output_table_per_gene <- merged_table_frequencies
    
    output_table_per_gene[-1][is.na(output_table_per_gene[-1])] <- 'Not expressed'
    
    write.csv(output_table_per_gene, file = 'LOI_matrix_with_dna_integration.csv', row.names = T, quote=FALSE)
    """
}
