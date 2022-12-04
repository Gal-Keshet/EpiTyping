process COMBINE_X_VCF {
    debug true

    publishDir "$params.outdir/vcf/xci", mode:'copy', enabled: "$params.keepInter" == true

    input:
    path vcf_files
    path vcf_index_files

    output:
    path "xci_samples_merged.vcf.gz"

    script:
    """
    bcftools merge\
    --threads $task.cpus\
    --output-type z\
    -m none\
    --output xci_samples_merged.vcf.gz\
    $vcf_files
    """
}

process OUTPUT_XCI_TABLE {
    debug true

    publishDir "$params.outdir/vcf/xci", mode:'copy', enabled: "$params.keepInter" == true

    input:
    path vcf_file

    output:
    path "xci_final_table.txt"

    script:
    """
    bcftools query\
    -H -f '%CHROM\t%POS\t%ID\t%INFO/Gene\t%INFO/COMMON\t%REF\t%ALT[\t%GT:%AD{0}:%AD{1}:%DP]\n'\
    --output xci_final_table.txt\
    $vcf_file
    """
}

process XCI_MATRIX {
    debug true
    publishDir "$params.outdir/results/xci", mode:'copy'
    input:
    path var_table
    path count_files

    output:
    path("XCI_status.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("stringr")
    library("dplyr")

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
    
    counts <- dataset[-c(1:4)]
    rpk= counts / dataset\$Length #first normalize for gene length
    tot_RPK_per_samp=apply(rpk,2,sum)/10^6 #Sum the normalized values
    tpm=t(t(rpk)/tot_RPK_per_samp) # Normalize to library size
    tpm[tpm < 1] <- 1
    tpm <- cbind(dataset[1:4], tpm)
    
    #Count n expressed genes per sample 
    per_chrm_n_expressed <- as.data.frame(tpm[-c(1,3,4)] %>%
      group_by(Chr) %>%
      summarise_all(., function(x) { sum(x>1)}))
    
    row.names(per_chrm_n_expressed) <- per_chrm_n_expressed\$Chr
    per_chrm_n_expressed <- per_chrm_n_expressed[-1]
    
    #read variant call file
    variants.file <- read.delim("$var_table", stringsAsFactors = F)
    #process dataframe and keep variants that are annotated to a gene's exon (Gene),
    #that have a dbSNP ID (RS) and are one base-pair long
    names(variants.file) <- str_replace(names(variants.file), 'X\\\\.+[0-9]+\\\\.', '')
    names(variants.file) <- str_replace(names(variants.file), '\\\\.GT.*', ':genotype;allele1;allele2;sum')
    
    variants.file <- variants.file[(grepl('[A-Z]+',variants.file\$Gene)) & 
                                     (grepl('^[ACGT]{,2}\$', variants.file\$REF)) &
                                     (grepl('^[ACGT]{,2}\$', variants.file\$ALT)), -c(3,5,6,7)]
    
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
    
      Sample[has_info] <- snp_info[, 5]
      Sample <- as.numeric(Sample)
      Sample[Sample == 0] <- NaN
      
      return(Sample)
    }
    
    allele_table <- cbind.data.frame(variants.file[1], apply(variants.file[-c(1:3)], 2, return_allelic_freq))
    
    colnames(allele_table) <- str_replace(colnames(allele_table), ':.*', '')
    
    allele_table <- as.data.frame(allele_table %>%
      group_by(CHROM) %>%
      summarise_all(., function(x) {sum(!is.na(x))}))
    
    allele_table <- allele_table[grepl('^chr[0-9]{1,2}\$|chrX', allele_table\$CHROM),]
    
    row.names(allele_table) <- allele_table\$CHROM
    
    allele_table <- allele_table[-1]
    
    per_chrm_n_expressed <- per_chrm_n_expressed[row.names(allele_table), colnames(allele_table), drop = FALSE]
    
    normalised_snp_sum <- allele_table / per_chrm_n_expressed
      
    x_a_ratios <- apply(normalised_snp_sum, 2, function(x) {
      x[row.names(normalised_snp_sum) == 'chrX'] / mean(x[row.names(normalised_snp_sum) != 'chrX'])
    })
    
    scores <- data.frame(list(x_vs_a_allelic_ratios = x_a_ratios))
    
    #Defining the thresholds for each X inactivation status
    output_table <- data.frame(X_status = apply(scores, 1, function(x) {
      if(x['x_vs_a_allelic_ratios'] > 0.3) return('XaXa')
      else if(x['x_vs_a_allelic_ratios'] > 0.12) return('XaXe')
      else return('XaXi') 
    }))

    write.csv(output_table, file = 'XCI_status.csv', row.names = T, quote=FALSE)
    """
}
