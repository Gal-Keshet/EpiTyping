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
    path x_escape

    output:
    path("XCI_status*")

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
        colnames(dataset) <- c('Gene.ID', 'Chr', 'Start', 'Length', 'gene_name', file)
        dataset\$Chr <- str_replace(dataset\$Chr, ';.*', '')
        dataset\$Start <- str_replace(dataset\$Start, ';.*', '')
      }
      
      # if the merged dataset does exist, append to it
      else if (exists("dataset")){
        temp_dataset <-read.table(file, header=TRUE, sep="\t")
        colnames(temp_dataset) <- c('Gene.ID', 'Chr', 'Start', 'Length', 'gene_name', file)
        temp_dataset\$Chr <- str_replace(temp_dataset\$Chr, ';.*', '')
        temp_dataset\$Start <- str_replace(temp_dataset\$Start, ';.*', '')
        dataset<-merge(dataset, temp_dataset, by=c('Gene.ID', 'Chr', 'Start', 'gene_name', 'Length'))
        rm(temp_dataset)
      }
    }
    
    colnames(dataset) <- str_replace(colnames(dataset), '.*/', '')
    colnames(dataset) <- str_replace(colnames(dataset), '_edited.*', '')
    
    counts <- dataset[-c(1:5)]
    rpk= counts / dataset\$Length #first normalize for gene length
    tot_RPK_per_samp=apply(rpk,2,sum)/10^6 #Sum the normalized values
    tpm=t(t(rpk)/tot_RPK_per_samp) # Normalize to library size
    tpm[tpm < 1] <- 1
    tpm <- cbind(dataset[1:5], tpm)

    #Extracting of Y-linked genes

    X_degenerate=c("SRY","RPS4Y1","ZFY","AMELY","TBL1Y","PRKY","USP9Y","DBY","UTY", "TMSB4Y","NLGN4Y","CYorf15A","CYorf15B","SMCY","EIF1AY","RPS4Y2")

    tpm_y = tpm[tpm\$gene_name %in% X_degenerate,]

    #Count n expressed genes per sample 
    per_chrm_n_expressed <- as.data.frame(tpm[-c(1,3,4,5)] %>%
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
                                     (grepl('^[ACGT]{,1}\$', variants.file\$REF)) &
                                     (grepl('^[ACGT]{,1}\$', variants.file\$ALT)), -c(3,5,6,7)]
    

    ##Removing X-escapees from analysis
    X_esacpe=read.csv("$x_escape", stringsAsFactors = F) #includes escaping and variable
    variants.file <- variants.file[!variants.file\$Gene %in% X_esacpe\$Gene.name,]

    #replace missing snps with nan for each sample, extract information about present values and calulate 
    #biallelic ratio (minor allele / major allele). than only keep snps that have a biallelic ratio 
    #bigger than 0.2, coverage of at-least 10 reads and are eterozygous (genotype=0/1)
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


    #reordering y_tpm and calculating mean y expression

    tpm_y <- tpm_y[-c(1:5)][colnames(allele_table)]

    mean_Y = apply(tpm_y, 2, mean, na.rm = T)
        
    scores <- data.frame(list(x_vs_a_allelic_ratios = x_a_ratios, Y_av_exp=mean_Y))

    #Defining the thresholds for each X inactivation status

    output_table <- data.frame(X_status = apply(scores, 1, function(x) {
     if(x['x_vs_a_allelic_ratios'] > 0.25 ) return('XaXa')
     else if(x['x_vs_a_allelic_ratios'] > 0.08) return('XaXe')
     else if(x['x_vs_a_allelic_ratios'] < 0.08 & x['Y_av_exp'] < 5) return('XaXi')
     else if(x['x_vs_a_allelic_ratios'] < 0.08 & x['Y_av_exp'] > 5) return('Male')
    }))

    write.csv(output_table, file = 'XCI_status.csv', row.names = T, quote=FALSE)

    #Biallelic genes and location
    allele_table <- cbind.data.frame(variants.file[c(1,3)], apply(variants.file[-c(1:3)], 2, return_allelic_freq)) 
    colnames(allele_table) <- str_replace(colnames(allele_table), ':.*', '') 

    allele_table <- as.data.frame(allele_table %>%
                                group_by(CHROM,Gene) %>%
                                summarise_all(.,function(x) {sum(!is.na(x))}))

    allele_table <- allele_table%>%filter(CHROM=="chrX")

    row.names(allele_table) <- allele_table\$Gene

    allele_table <- allele_table[-c(1:2)]

    tpm <- tpm[tpm\$gene_name %in% row.names(allele_table),]

    tpm <- tpm[!duplicated(tpm\$gene_name),]

    row.names(tpm) <- tpm\$gene_name

    tpm <- tpm[row.names(allele_table),]

    tpm <- tpm[c('Start', colnames(allele_table))]

    allele_table[] <- apply(allele_table, 2, function(x) {paste0('Biallelic by ', x, ' SNPs')})

    #keeping allele scores only for expressed genes
    allele_table[tpm[-1] == 1] <- 'Unexpressed'

    #changing biallelic score matrix nan values with 0
    allele_table[allele_table == 'Biallelic by 0 SNPs'] <- 'Uninformative'

    allele_table <- merge(tpm[1], allele_table, by = 'row.names')

    row.names(allele_table) <- allele_table\$Row.names

    allele_table <- allele_table[-1]

    allele_table <- allele_table[order(as.numeric(allele_table\$Start)),]

    write.csv(allele_table, file = 'XCI_status_per_gene.csv', row.names = T, quote=FALSE)
    """
}
