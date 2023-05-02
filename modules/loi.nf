process COMBINE_LOI_VCF {
    debug true

    publishDir "$params.outdir/vcf/loi", mode:'copy', enabled: "$params.keepInter" == true

    input:
    path vcf_files
    path vcf_index_files
    path imprinted_intervals

    output:
    path "loi_samples_merged.vcf*"

    script:
    """
    bcftools merge\
    --threads $task.cpus\
    --output-type z\
    --regions-file $imprinted_intervals\
    -m none\
    --output loi_samples_merged.vcf.gz\
    $vcf_files

    tabix -f loi_samples_merged.vcf.gz
    """
}

process OUTPUT_LOI_TABLE {
    debug true

    publishDir "$params.outdir/vcf/loi", mode:'copy', enabled: "$params.keepInter" == true

    input:
    path vcf_file
    path imprinted_intervals    

    output:
    path "loi_final_table.txt"

    script:
    """
    bcftools query\
    -H -f '%CHROM\t%POS\t%ID\t%INFO/Gene\t%INFO/COMMON\t%REF\t%ALT[\t%GT:%AD{0}:%AD{1}:%DP]\n'\
    --regions-file $imprinted_intervals\
    --output loi_final_table.txt\
    ${vcf_file[0]}
    """
}

process LOI_MATRIX {
    debug true
    publishDir "$params.outdir/results/loi", mode:'copy'
    input:
    path var_table
    path count_files
    path imprinted_genes_file

    output:
    tuple path("LOI_matrix.csv"), path("LOI_per_locus_matrix.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library("stringr")
    library("dplyr")
      
    #upload imprinted genes
    gene_locations <- read.csv("$imprinted_genes_file")
    #read variant call file
    variants.file <- read.delim("$var_table", stringsAsFactors = F)
    #process dataframe and keep variants that are annotated to a gene's exon (Gene),
    #that have a dbSNP ID (RS) and are one base-pair long
    names(variants.file) <- str_replace(names(variants.file), 'X\\\\.+[0-9]+\\\\.', '')
    names(variants.file) <- str_replace(names(variants.file), '\\\\.GT.*', ':genotype;allele1;allele2;sum')
    variants.file <- variants.file[(grepl('[A-Z]+',variants.file\$Gene)) & (grepl('^[ACGT]{,1}\$', variants.file\$REF)) & (grepl('^[ACGT]{,1}\$', variants.file\$ALT)), -c(3,5,6,7)]    
    
    rownames(variants.file) <- NULL
    
    variants.file <- left_join(gene_locations, variants.file, by = 'Gene', multiple = 'all')
    
    variants.file[is.na(variants.file)] <- './.:.:.:.'
    
    variants.file <- variants.file[,-c(4:5)]

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
    
    allele_table <- cbind.data.frame(variants.file[c(1:3)], apply(variants.file[-c(1:3)], 2, return_allelic_freq))
    
    colnames(allele_table) <- str_replace(colnames(allele_table), ':.*', '')
        
    allele_table <- as.data.frame(allele_table %>%
      group_by(Gene = factor(Gene, levels = unique(Gene)), Location, Imprinting.related.disease) %>%
      summarise_all(.,function(x) {sum(!is.na(x))}))

    row.names(allele_table) <- allele_table\$Gene
    
    allele_table <- allele_table[,-1]
    
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
    tpm <- tpm[tpm\$gene_name %in% row.names(allele_table),]
    
    row.names(tpm) <- tpm\$gene_name
    
    tpm <- tpm[row.names(allele_table),]
    tpm <- tpm[colnames(allele_table)[-c(1:2)]]

    allele_table[-c(1:2)] <- apply(allele_table[-c(1:2)], 2, function(x) {paste0('Biallelic by ', x, ' SNPs')})
    
    #keeping allele scores only for expressed genes
    allele_table[-c(1:2)][tpm < 1] <- 'Unexpressed'
    
    #changing biallelic score matrix 0 values with 'uninformative'
    allele_table[-c(1:2)][allele_table[-c(1:2)] == 'Biallelic by 0 SNPs'] <- 'Uninformative'

    write.csv(allele_table, file = 'LOI_matrix.csv', row.names = T, quote=FALSE)
    
    #Per regio LOI
    allele_table <- as.data.frame(allele_table[-2] %>%
     group_by(Location = factor(Location, levels = unique(Location))) %>%
     summarise_all(., function(x) {sum(x != 'Uninformative' & x != 'Unexpressed')}))    

    row.names(allele_table) <- allele_table\$Location
    allele_table <- allele_table[-1]

    allele_table[] <- apply(allele_table, 2, function(x) {paste0(x, ' genes with biallelic expression')})
    
    write.csv(allele_table, file = 'LOI_per_locus_matrix.csv', row.names = T, quote=FALSE)
    """
}
  
