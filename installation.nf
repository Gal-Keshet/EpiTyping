params.humanGTF = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz'
params.humanFasta = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz'
params.mouseGTF = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz'
params.mouseFasta = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz'
params.dbSNP = 'https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz'
params.dbSNP_index = 'https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.tbi'
params.readlength=100
params.remap_file="$PWD/genome_files/remapNCBI.txt"

process DOWNLOAD_HUMAN_GTF {
    debug true
    publishDir "$PWD/genome_files/star_ref", mode: 'copy'     
    input:
    val gtf
        
    output:
    path "gencode.v42.primary_assembly.annotation.gtf"

    """
    wget -O gencode.v42.primary_assembly.annotation.gtf.gz $gtf
    gunzip gencode.v42.primary_assembly.annotation.gtf.gz
    """ 
}

process CREATE_INTERVALS {
        debug true
        publishDir "$PWD/genome_files", mode: 'copy'
        
        input:
        path gtf
         
        output: 
        path "GRCh38_exome.bed.gz"

        """
        awk '{if(\$3=="exon") {print \$1"\\t"\$4-1"\\t"\$5+1"\\t"substr(\$16,2,length(\$16)-3)}}' $gtf | grep ^chr | sort -k 1,1 -k2,2n | bgzip > GRCh38_exome.bed.gz
        """ 
}

process INDEX_INTERVALS {
        debug true
        publishDir "$PWD/genome_files", mode: 'copy'
        
        input:
        path intervals
        
        output: 
        path "GRCh38_exome.bed.gz.tbi"

        """
        tabix $intervals
        """ 
}

process DOWNLOAD_HUMAN_FASTA {
    debug true    
    publishDir "$PWD/genome_files/star_ref", mode: 'copy' 
        
    input:
    val fasta
  
    output:
    path "GRCh38.primary_assembly.genome.fa"

    """
    wget -O GRCh38.primary_assembly.genome.fa.gz $fasta
    gunzip GRCh38.primary_assembly.genome.fa.gz
    """ 
}

process INDEX_FASTA {
        debug true
        publishDir "$PWD/genome_files/star_ref", mode: 'copy'
        
        input:
        path fasta 
         
        output: 
        path "${fasta}.fai"

        """
        samtools faidx $fasta 
        """  
}


process CREATE_DICTIONARY {
        debug true
        label 'var_call'
        publishDir "$PWD/genome_files/star_ref", mode:'copy'
        memory '8GB'
        
        input:
        path fasta
        path index
         
        output: 
        path '*'

        """
        gatk CreateSequenceDictionary -R $fasta 
        """  
}

process GENERATE_HUMAN_STAR_INDEX {
    debug true         
    publishDir "$PWD/genome_files/", mode:'copy'
    memory '40GB'

    input:
    path fasta
    path gtf

    output:
    path "star_index/*"
        
    """
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $params.readlength 
    """
}

process DOWNLOAD_MOUSE_GTF {
    debug true
    input:
    val gtf

    output:
    path "gencode.vM31.primary_assembly.annotation.gtf"

    """
    wget -O gencode.vM31.primary_assembly.annotation.gtf.gz $gtf
    gunzip gencode.vM31.primary_assembly.annotation.gtf.gz
    """
}

process DOWNLOAD_MOUSE_FASTA {
    debug true
    input:
    val fasta

    output:
    path "GRCm39.primary_assembly.genome.fa"

    """
    wget -O GRCm39.primary_assembly.genome.fa.gz $fasta
    gunzip GRCm39.primary_assembly.genome.fa.gz
    """
}

process GENERATE_MOUSE_STAR_INDEX {
    debug true
    publishDir "$PWD/genome_files/", mode: 'copy'
    memory '40GB'

    input:
    path fasta
    path gtf

    output:
    path "star_mouse_index/*"

    """
    STAR --runThreadN $task.cpus --runMode genomeGenerate --genomeDir star_mouse_index --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang $params.readlength
    """
}

process DOWNLOAD_DBSNP {
    debug true

    input:
    val dbsnp
    val dbsnp_index

    output:
    tuple path("dbSNP155.vcf.gz"), path("dbSNP155.vcf.gz.tbi")

    """
    wget --no-check-certificate -O dbSNP155.vcf.gz $dbsnp
    wget --no-check-certificate -O dbSNP155.vcf.gz.tbi $dbsnp_index
    """
}

process RENAME_DBSNP {
        debug true
        publishDir "$PWD/genome_files", mode: 'copy' 
        
        input:
        tuple path(dbSNP), path(dbSNPindex)
        path remapNCBI

        output:
        path "dbSNP155Renamed.vcf*"

        """
        bcftools annotate --threads $task.cpus --output-type z --rename-chrs $remapNCBI --output dbSNP155Renamed.vcf.gz $dbSNP
        tabix dbSNP155Renamed.vcf.gz
        """ 
}

process GENERATE_BWA_INDEX {
    debug true
    memory '40GB'

    publishDir "$PWD/genome_files/bwa_ref", mode: 'copy'

    input:
    path fasta

    output:
    path "*"

    """
    bwa index $fasta
    """
}
 
workflow {
    
    human_gtf_ch = DOWNLOAD_HUMAN_GTF(params.humanGTF)

    exon_intervals_ch = CREATE_INTERVALS(human_gtf_ch)
    
    INDEX_INTERVALS(exon_intervals_ch)

    human_fasta_ch = DOWNLOAD_HUMAN_FASTA(params.humanFasta)    

    human_fasta_index_ch = INDEX_FASTA(human_fasta_ch)

    CREATE_DICTIONARY(human_fasta_ch, human_fasta_index_ch)

    GENERATE_HUMAN_STAR_INDEX(human_fasta_ch, human_gtf_ch)

    mouse_gtf_ch = DOWNLOAD_MOUSE_GTF(params.mouseGTF)

    mouse_fasta_ch = DOWNLOAD_MOUSE_FASTA(params.mouseFasta)
 
    GENERATE_MOUSE_STAR_INDEX(mouse_fasta_ch, mouse_gtf_ch)

    dbsnp_ch = DOWNLOAD_DBSNP(params.dbSNP, params.dbSNP_index)

    RENAME_DBSNP(dbsnp_ch, params.remap_file)

    GENERATE_BWA_INDEX(human_fasta_ch) 
    
}
