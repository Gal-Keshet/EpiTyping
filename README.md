# EpiTyping
EpiTyping is a tool for detecting imprinting and X-chromosome inactivation status from RNA-seq
![](https://github.com/Gal-Keshet/EpiTyping/files/10101966/scheme.pdf)

## Pre-requisites 
Nextflow and Singularity (or apptainer, the newer name for the Singularity project) are the only pre-requisites for the EpiTyping tool. Install both if needed and make sure they are properly running on your system. If the following commands do not generate any error message you are good to go.
```bash
nextflow run hello  # test that nextflow is working
singularity help  # test that singularity is working
```

## Installation

1. Download the project directory:
```bash
git clone https://github.com/Gal-Keshet/EpiTyping.git  # clone the project using git
```
2. Execute the script named installation.nf which is responsible for setting up all the reference data and will complete the installation (this might take a while).
```bash
nextflow run /path/to/EpiTyping/installation.nf [parameters]  # run the installation script
```
### Basic parameters
-profile: Choose the executor profile between a local workstation or usage on a SLURM cluster (standard/cluster, default: standard).

--cpu: The number of threads for multi-threading (int, default 8).

--readlength: The expected Illumina read length for optimal alignment by STAR (int, default 100).
### Advanced parameters (Use with cation. Might lead to unexpected resaults if the alternative gtf/fasta files have a different format than the default files.)
--humanFasta: url for a human fasta file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.primary_assembly.genome.fa.gz)

--humanGTF: url for a human gtf file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz)

--mouseFasta: url for a mouse fasta file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz)

--mouseGTF: url for a mouse gtf file (default: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gtf.gz)

--dbSNP: url for a dbSNP vcf file (default: https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz)

--dbSNP_index: url for a dbSNP index file (default: https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.tbi)

--remap_file: Path to the NCBI chromosome remapping names. (default: "$projectDir/genome_files/remapNCBI.txt")
## Usage

```bash
nextflow run /path/to/EpiTyping/main.nf --fastq_folder /path/to/fastq_dir --single true  # for single-end reads
nextflow run /path/to/EpiTyping/main.nf --fastq_folder /path/to/fastq_dir --single false  # for paired-ends reads (default)
```
### Basic parameters
-profile: (standard/cluster, default: standard).

--cpu: (int, default 8).

--keepInter: Whether to keep intermediate alignment and VCF files (true/false, default: false). 

--filterMouse: Whether to perform mouse contamination cleanup (true/false, default: true).

Example for a paired-ends RNA-seq run, using 4 CPUs, keeping intermediate files:
```bash
nextflow run /path/to/RNA2CM --fastq esc_1.fastq.gz --fastq2 esc_2.fastq.gz --cpu 4 --keepInter true 
```

Example for a single-end RNA-seq run, skipping mouse read filtration and running on a SLURM cluster (nextflow manages batch jobs, so no need to use sbatch):
```bash
nextflow run /path/to/RNA2CM --fastq SRR1234567.fastq.gz --filterMouse false -profile cluster
```
