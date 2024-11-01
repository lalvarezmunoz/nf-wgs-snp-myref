# nf-wgs-snp-myref
Nextflow pipeline for whole genome sequencing (WGS) and variant calling of bacterial genomes, using a custom reference.


## Description

This is a Nextflow pipeline designed to simplify WGS and variant calling analyses of bacterial genomes. It uses Docker containers making installation trivial and results highly reproducible.

The user provides the path to the FASTQ files and the custom reference genome (.fna format).

First, sequences are aligned agains the reference genome using BWA-MEM. The aligned reads are then sorted using Samtools, duplicate reads are marked with Picard and the BAM files are indexed. Variant calling is performed using Freebayes.


## Requirements

- Nextflow version 24.04.4.5917
- Docker version 27.2.0


## Usage

Step 1: create a **results** folder in the working folder

Step 2: create a sample file

Step 3: run standard CLI code:

```bash
nextflow run main.nf --input <samples.csv> --genome_id <acc_number> --genome_fna <path_to_fna_file> --profile prod -resume
```

Step 4: output files will be created in the **results** folder.

*Example:*

```bash
nextflow run main.nf --input mysamples.csv --genome_id myreference --genome_fna </myfolder1/myreferencesfolder/refrence.fna --profile prod -resume
```

---

## Input

### Sample file

CSV file with the following format:
id,R1,R2

- id: sample identification
- R1: path to FASTQ file with sequencing results R1 reads
- R2: path to FASTQ file with sequencing results R2 reads

Example:

```
id,R1,R2
sample1,/mnt/c/myfolder/sample1_S1_L001_R1_001.fastq.gz,/mnt/c/myfolder/sample1_S1_L001_R2_001.fastq.gz
sample2,/mnt/c/myfolder/sample2_S2_L001_R1_001.fastq.gz,/mnt/c/myfolder/sample2_S2_L001_R2_001.fastq.gz
[...]
```

### Genome ID

Name to identifiy the reference genome used.

### Genome FNA

Absolute path to the .fna file of the reference genome.


## Output

### Logs

- BWA-MEM statistics: information about aligned reads
- Picard metrics: information about duplicate reads

### Bam files

Indexed aligned reads in .bam and .bai format.

### Vcf files

Variant calling report files, containing the following information:

```
#CHROM	POS	ID	REF	ALT	QUAL    FILTER	INFO	FORMAT	unknown
```

Information lines describing the *FILTER*, *INFO* and *FORMAT* entries used in the body of the VCF file are included in the meta-information section (file header).


---

## Tools and references

- BWA-MEM (https://bio-bwa.sourceforge.net/)
- samtools/flagstat (http://www.htslib.org/doc/samtools-flagstat.html)
- samtools/sort (http://www.htslib.org/doc/samtools-sort.html)
- Picard/MarkDuplicates (https://broadinstitute.github.io/picard/)
- Freebayes (https://github.com/freebayes/freebayes)
- VCFtools (https://vcftools.github.io/index.html)


## Acknowledgements

Thanks to [loAlon](https://github.com/loalon) for his contribution and guidance.


## Citing this pipeline

Please, refer to the GitHub repository when using this pipeline:

> L. Alvarez, nf-wgs-snp-myref pipeline, (2024), GitHub repository, https://github.com/lalvarezmunoz/nf-wgs-snp-myref
