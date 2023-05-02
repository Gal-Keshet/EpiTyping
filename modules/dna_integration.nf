process TRIM_DNA_FASTQ {
    time '1d'

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
        ${fastq[0]} ${fastq[2]} | samtools view -F 4 -b -@ $task.cpus -o ${sampleId}_aligned.bam -
        """
    }
    else {
        """
        bwa mem -M -t $task.cpus\
        $ref\
        $fastq | samtools view -F 4 -b -@ $task.cpus -o ${sampleId}_aligned.bam -
        """
    }
}

process ALIGN_DNA_MOUSE {
    time '1d'

    debug true

    input:
    tuple val(sampleId), path(fastq)
    path ref_mouse
    path index_mouse

    output:
    tuple val("$sampleId"), path("${sampleId}_mouse_aligned.bam")

    script:
    def single = fastq instanceof Path
    if( !single ) {
        """
        bwa mem -M -t $task.cpus\
        $ref_mouse\
        ${fastq[0]} ${fastq[2]} | samtools view -F 4 -b -@ $task.cpus -o ${sampleId}_mouse_aligned.bam -
        """
    }
    else {
        """
        bwa mem -M -t $task.cpus\
        $ref_mouse\
        $fastq | samtools view -F 4 -b -@ $task.cpus -o ${sampleId}_mouse_aligned.bam -
        """
    }
}

process SORT_MOUSE_BAM {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(mouse_bam)

    output:
    path "${sampleId}_mouse_sorted.bam"

    script:
    """
    java -jar /picard.jar SortSam\
    I=$mouse_bam O=${sampleId}_mouse_sorted.bam\
    SORT_ORDER=coordinate CREATE_INDEX=true
    """
}


process CUT_BAM {
    time '1d'

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

process DNA_XENOFILTER {
    debug true

    publishDir "$params.outdir/bam/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(human_bam_file) 
    path mouse_bam_file


    output:
    tuple val("$sampleId"), path("Filtered_bams/${sampleId}_sorted_Filtered.bam")

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
    --info-key CNN_1D\
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
    path "dna.vcf*"

    script:
    """
    bcftools view --threads $task.cpus\
    -i 'INFO/RS!="." & FILTER="PASS" & FMT/DP > 9'\
    --regions-file $imprinted_intervals\
    --output-type z\
    --output dna.vcf.gz\
    ${vcf_file[0]}

    tabix -f dna.vcf.gz 
    """
}   

process OUTPUT_ASE {
    debug true
    label 'var_call'

    publishDir "$params.outdir/vcf/$sampleId", mode:'copy', enabled: "$params.keepInter" == true

    input:
    tuple val(sampleId), path(rna_bam_file)
    path dna_vcf_file
    path genome
    path index
    path dict

    output:
    path "${sampleId}_ASE"

    script:
    """
    gatk ASEReadCounter\
    -R $genome\
    -I $rna_bam_file\
    -V ${dna_vcf_file[0]}\
    -O ${sampleId}_ASE
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
    ${vcf_file[0]}
    """
}

process LOI_MATRIX_DNA_INTEGRATION {
    debug true
    publishDir "$params.outdir/results/loi", mode:'copy'
    input:
    path dna_table
    path ase_samples
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
    
    dna.file <- dna.file[(grepl('[A-Z]+',dna.file\$Gene)) & (grepl('^[ACGT]{,1}\$', dna.file\$REF)) & (grepl('^[ACGT]{,1}\$', dna.file\$ALT)), -5]
    
    rownames(dna.file) <- NULL
    
    dna.file <- dna.file[1:7]

    heterozygote_pos <- inner_join(gene_locations, dna.file, by = 'Gene', multiple = 'all')
        
    ####read rna allele expression####
    ase_files <- unlist(str_split("$ase_samples", pattern = ' '))

    ase_data <- list()

    for(ase_sample in ase_files) {
      ase_data[[str_replace_all(ase_sample, '.*/|_ASE', '')]] = read.delim(ase_sample,  stringsAsFactors = F)
    }

    ase_merged_data <- list()

    ####allelic freq function####
    return_allelic_freq <- function(snp) {
      if(as.numeric(snp['totalCount']) < 10) return(NaN)
      else {
        return(min(as.numeric(snp[c('refCount', 'altCount')])) / max(as.numeric(snp[c('refCount', 'altCount')])))
      }
    }

    ####looping through ase samples####
    for(ase_sample in names(ase_data)) {

      ase_data[[ase_sample]] = ase_data[[ase_sample]][1:8]
      
      colnames(ase_data[[ase_sample]])[1:5] = c('CHROM','POS','ID','REF','ALT')
      
      #merge ase samples with dna
      ase_merged_data[[ase_sample]] = inner_join(heterozygote_pos, ase_data[[ase_sample]], multiple = 'all')
      
      #call allelic freq
      ase_merged_data[[ase_sample]]\$freq <- apply(ase_merged_data[[ase_sample]], 1, return_allelic_freq)
      
      #determine for each snp whether its biallelic, monoallelic or uninformative
      ase_merged_data[[ase_sample]]\$allelic <- sapply(ase_merged_data[[ase_sample]]\$freq, function(x) {
        if(is.na(x)) return(0)
        else if (x < 0.1) return(1)
        else if (x > 0.2) return(2)
        else return(0)
      })
      
      #count the number of monoallelic/biallelic/uninformative SNPs per gene
      ase_merged_data[[ase_sample]] <- as.data.frame(ase_merged_data[[ase_sample]] %>%
                                      group_by(Gene, Location, Imprinting.related.disease) %>%
                                      summarise(n.biallelic = sum(allelic == 2), n.monoallelic = sum(allelic == 1), n.uninformative = sum(allelic == 0)))
      
      #summarize  results
      ase_merged_data[[ase_sample]]\$ase_sample <- apply(ase_merged_data[[ase_sample]][-c(1:3)],1, function(x) {
        if(x['n.monoallelic'] == 0 & x['n.biallelic'] == 0) {return('Uninformative')}
        else if (x['n.biallelic'] == 0) {return(paste0('Monoallelic by ', x['n.monoallelic'], ' SNPs'))}
        else {return(paste0('Biallelic by ', x['n.biallelic'], ' SNPs'))}
      })
      
      colnames(ase_merged_data[[ase_sample]])[7] <- ase_sample
      
      ase_merged_data[[ase_sample]] <- ase_merged_data[[ase_sample]][c(1,2,3,7)]
    }

    aes_joint_df <- Reduce(function(x, y) merge(x, y, all = T), ase_merged_data)

    aes_joint_df <- inner_join(gene_locations, aes_joint_df)

    row.names(aes_joint_df) <- aes_joint_df\$Gene

    aes_joint_df <- aes_joint_df[-1]

    ###count files
    files_vector <- unlist(str_split("$count_files", pattern = ' '))    
    #merging counts
    for (file in files_vector){
      
      # if the merged dataset doesn't exist, create it
      if (!exists("dataset")){
        dataset <- read.table(file, header=TRUE, sep="\t")
        dataset <- dataset[-c(2:3)]
        colnames(dataset) <- c('Gene.ID', 'Length', 'gene_name', file)
      }
      
      # if the merged dataset does exist, append to it
      else if (exists("dataset")){
        temp_dataset <-read.table(file, header=TRUE, sep="\t")
        temp_dataset <- temp_dataset[-c(2:3)]
        colnames(temp_dataset) <- c('Gene.ID', 'Length', 'gene_name', file)
        dataset<-merge(dataset, temp_dataset, by=c('Gene.ID', 'gene_name', 'Length'))
        rm(temp_dataset)
      }
    }
    
    colnames(dataset) <- str_replace(colnames(dataset), '.*/', '')
    colnames(dataset) <- str_replace(colnames(dataset), '_edited.*', '')
    
    #######tpm calculation######
    counts <- dataset[-c(1:3)]
    rpk= counts / dataset\$Length #first normalize for gene length
    tot_RPK_per_samp=apply(rpk,2,sum)/10^6 #Sum the normalized values
    tpm=t(t(rpk)/tot_RPK_per_samp) # Normalize to library size
    tpm <- cbind(dataset[1:3], tpm)

    #subsetting tpm for imprinted genes & changing order of biallelic score matrix to tpm matrix order
    tpm <- tpm[tpm\$gene_name %in% row.names(aes_joint_df),]

    row.names(tpm) <- tpm\$gene_name

    tpm <- tpm[row.names(aes_joint_df),]
    tpm <- tpm[colnames(aes_joint_df)[-c(1:2)]]

    #keeping allele scores only for expressed genes
    #new
    aes_joint_df[-c(1:2)][tpm < 1] <- 'Not expressed'

    aes_joint_df[-c(1:2)][is.na(aes_joint_df[-c(1:2)])] <- 'Uninformative'

    aes_joint_df <- aes_joint_df[apply(aes_joint_df[-c(1:2)], 1, function(x) {sum(x == 'Not expressed') + sum(x == 'Uninformative') != length(x)}),]

    output_table_per_gene <- aes_joint_df
    
    write.csv(output_table_per_gene, file = 'LOI_matrix_with_dna_integration.csv', row.names = T, quote=FALSE)
    """
}
