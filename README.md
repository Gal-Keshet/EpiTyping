# EpiTyping
EpiTyping is a tool for detecting imprinting and X-chromosome inactivation status from RNA-seq
![scheme_cropped](https://user-images.githubusercontent.com/112553439/204306399-e99da88e-0d30-4665-bfdb-e4329b1e8f53.jpg)

## Pre-requisites 
Nextflow and Singularity (or apptainer, the new name for the Singularity project) are the only pre-requisites for the EpiTyping tool. Install both if needed and make sure they are properly running on your system. 

* For nextflow installation refer to https://www.nextflow.io/ 

* for apptainer/singularity installation refer to https://apptainer.org/ 

If the following commands do not generate any error message you are good to go.
```bash
nextflow run hello  # test that nextflow is working

singularity help  # test that singularity is working
```

## Installation

1. Download the project directory:
```bash
git clone https://github.com/Gal-Keshet/EpiTyping.git  # clone the project using git
```
2. Download docker images to run the software:
```bash
cd /path/to/EpiTyping

singularity pull docker://broadinstitute/gatk

singularity pull docker://galclbl/epityping
```
3. Execute the script named installation.nf which is responsible for setting up all the reference data and will complete the installation (this might take a while).
```bash
nextflow run installation.nf [parameters]  # run the installation script
```

### Basic parameters
-profile: Choose the executor profile between a local workstation or usage on a SLURM cluster (standard/cluster, default: standard).

--cpu: The number of threads for multi-threading (int, default 8).

--mem: Memory limit for the run ([int].GB, default: 64.GB)

--readlength: The expected Illumina read length for optimal alignment by STAR (int, default 100).

### Advanced parameters 

* **Use with caution**. Might lead to unexpected results if the alternative gtf/fasta files have a different format than the default files.

--humanFasta: url for a human fasta file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz)

--humanGTF: url for a human gtf file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz)

--mouseFasta: url for a mouse fasta file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz)

--mouseGTF: url for a mouse gtf file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz)

--dbSNP: url for a dbSNP vcf file (default: https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz)

--dbSNP_index: url for a dbSNP index file (default: https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.tbi)

--remap_file: Path to the NCBI chromosome remapping names (default: "$projectDir/genome_files/remapNCBI.txt").

# Usage
```bash
bgzip /path/to/fastq_dir/* # Use if fastq files are not compressed.

nextflow run /path/to/EpiTyping/main.nf --fastq_folder /path/to/fastq_dir --single true [parameters] # for single-end reads (all files in fastq folder need to be single ended)

nextflow run /path/to/EpiTyping/main.nf --fastq_folder /path/to/fastq_dir --single false [parameters] # for paired-ends reads (default. All files in fastq folder need to be paired ended)
```

* **Important!** All RNA (or DNA) fastq files need to be gz-compressed. Single ended fastq files need to have the following namming format: sample-name.fastq.gz, and paired ended fastq files need to have the following format: sample-name_1.fastq.gz/sample-name_2.fastq.gz

### Basic parameters

--fastq_folder: Path to a folder containing the fastq files (default: $PWD/fastq). The folder can contain multiple fastq files but they all need to represent samples which were prepared with the same library composition (e.g all paired and unstranded; all single and unstranded etc.)

--multiple_samples: Whether there is more than one sample (true/false, default: false)

--single: Whether the RNA-seq library is single ended (true/false, default: false)

-profile: (standard/cluster, default: standard).

--cpu: (int, default 12).

--mem: Memory limit for the run ([int].GB, default: 64.GB)

--keepInter: Whether to keep intermediate alignment and VCF files (true/false, default: false). 

--mouse_feeders: Whether to perform mouse contamination cleanup (true/false, default: false).

--strandness: Whether the library is unstranded or strand-specific (int, 0: unstranded; 1: forward strand; 2: reverse strand, default: 0)

--dna_sample_name: Name of DNA-seq sample without the .fastq.gz extension. Use if a DNA-seq for the examined cell line is available (default: "" for no integration).

--dna_fastq_folder: Path to dna fastq files if DNA-seq for the examined cell line is available. Used together with --dna_sample_name parameter (default: $PWD/dna_fastq). 

--dna_single: Whether the DNA-seq library is single ended (true/false, default: false)

* If a DNA-seq sample is given for integration, make sure **all** the RNA-seq samples in the fastq folder are from the same cell line

--outdir: Existing or to-be-created path for output files (default: $PWD/output)

### Advanced parameters (Use with caution. see note on advanced parameters for installation)

--common_adapters: Full path to a fasta file containing ilumina adapter sequences to be used with the Trimmomatic software. We provide such file in the genome_files directory (default: $projectDir/genome_files/CommonAdapters.fa).

--index: Full path to the human STAR index folder (default: $projectDir/genome_files/star_index)

--mouse_index: Full path to the mouse STAR index folder (default: $projectDir/genome_files/star_mouse_index)

--gtf: Full path to a human gtf file (default: $projectDir/genome_files/star_ref/gencode.v42.primary_assembly.annotation.gtf)

--ref: Full path to a human fasta file (default: $projectDir/genome_files/star_ref/GRCh38.primary_assembly.genome.fa)

--ref_dict: Full path to a human fasta dictionary file (created and to be used by gatk) (default: $projectDir/genome_files/star_ref/GRCh38.primary_assembly.genome.dict)

--dbsnp: Full path to a dbsnp vcf file (default: $projectDir/genome_files/dbSNP155Renamed.vcf.gz)

--exon_intervals: Full path to a bed file containing exon intervals (default: $projectDir/genome_files/GRCh38_exome.bed.gz)

--imprinted_intervals: Full path to a sorted bed file containing imprinted genes intervals (default: $projectDir/genome_files/imprinted_genes_sorted.bed)

--imprinted_genes: Full path to a csv file containing two columns (Gene, Location) with imprinted genes name and their cytogenetic band (default: $projectDir/genome_files/imprinted_genes_location.csv)

--headers: Full path to a vcf headers file for vcf annotation with gene names (default: $projectDir/genome_files/headers_for_annotation.txt)

--bwa_index: Full path to the BWA index files (default: $projectDir/genome_files/bwa_ref/*)

# Examples

**Let there be a folder in a path /User/gal/epigenetic_analysis/fastq within a slurm cluster which contains fastq files from 2 paired-end samples that are unstranded: ERR3466738_1.fastq.gz;  ERR3466738_2.fastq.gz; ERR3466740_1.fastq.gz; ERR3466740_2.fastq.gz. Suppose these samples grew on MEFs.**

* To run the analysis for these samples and keep intermediate files, use the following command (after running the installation script):

```bash
nextflow run /path/to/EpiTyping/main.nf\
 --fastq_folder /User/gal/epigenetic_analysis/fastq\
 --multiple_samples true\
 -profile cluster\
 --keepInter true\
 --mouse_feeders true\
 --outdir /User/gal/epigenetic_analysis/output
```
* Suppose you also have a folder /User/gal/epigenetic_analysis/dna_fastq which contains fastq files from  whole genome sequencing of H1 ESC line: SRR2070629_1.fastq.gz; SRR2070629_2.fastq.gz. Considering that both ERR3466740 and ERR3466738 are samples taken from H1 ESC, you can integrate the DNA-seq with the RNA-seq to detect monoallelic as well as biallelic expression from the RNA-seq samples with the following command:

```bash
nextflow run /path/to/EpiTyping/main.nf\
 --fastq_folder /User/gal/epigenetic_analysis/fastq\
 --multiple_samples true\
 --dna_fastq_folder /User/gal/epigenetic_analysis/dna_fastq\
 --dna_reference SRR6377128\
 -profile cluster\
 --keepInter true\
 --mouse_feeders true\
 --strandness 1\
 --outdir /User/gal/epigenetic_analysis/output
```

# Output
Output files are saved to the outdir folder (default: $PWD/output)

Under the outdir folder, a results folder will also be created with two subdirectories: loi and xci.

Under loi, if there is no DNA integration, two files will be created:

 - LOI_matrix.csv: a table where each row correspondes to a gene that exists in the imprinted_genes_location.csv file under the genome_files folder, and each column      corresponds to a sample. In this table, 0 corresponds to an uninformative gene, "not expressed" corresponds to a gene that was not expressed and an int > 0            corresponds to the number of biallelic SNPs detected in the given gene.
 
 - LOI_per_locus_matrix.csv: a table where each row correspondes to an imprinted region that exists in the imprinted_genes_location.csv file under the genome_files        folder, and each column corresponds to a sample. In this table, 0 corresponds to an uninformative region and an int > 0 corresponds to the number of biallelic genes    detected in the given region.
 
If there is DNA integration, one file will be created in the loi folder:

 - LOI_matrix_with_dna_integration.csv: a table where each row correspondes to a gene that exists in the imprinted_genes_location.csv file under the genome_files          folder, and that showed at least one heterozygote position in one of its exons. Each column corresponds to a sample. In this table, 0 corresponds to a gene which      showed no biallelic expression of any heterozygote SNP detected in the DNA-seq, "not expressed" corresponds to a gene that was not expressed and an int > 0            corresponds to the number of biallelic SNPs detected in the given gene (only SNPs which were detected in the DNA).
 
 Under the xci folder, one file will appear:
 
  - XCI_status.csv: a table where each row correspondes to a sample and one column which will either output "XaXa" for no X-chromosome inactivation (XCI), "XaXe" for       XCI erosion, "XaXi" for intact XCI or "Male" in case the RNA-seq sample belongs to a male cell line.
 
* If keepInter paremeter is true, subdirectories under the outdir folder will be created with the intermediate results. 



