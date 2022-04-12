# *P*ersonalised *A*SE *C*aller (PAC)
---------------------------------------

Author: Anna Saukkonen

anna.saukkonen@gmail.com

See our preprint [PAC: Highly accurate quantification of allelic gene expression for population and disease genetics](https://www.biorxiv.org/content/10.1101/2021.07.13.452202v1) for additional information

## TABLE OF CONTENTS
1. [Introduction](#INTRODUCTION)  
2. [Installation and running](#INSTALLATION-AND-RUNNING)  
3. [Options](#OPTIONS)
4. [Output](#OUTPUT)
5. [Test Dataset](#TEST-DATASET)

<!-- toc -->

## INTRODUCTION:

Allele-specific expression (ASE) is the imbalanced expression of the two alleles of a gene. While 
many genes are expressed equally from both alleles, gene regulatory differences driven by
genetic changes (i.e. regulatory variants) frequently cause the two alleles to be expressed at
different levels, resulting in allele-specific expression patterns. The detection of ASE events 
relies on accurate alignment of RNA-sequencing reads, where challenges still remain. This pipeline 
has been created to adjust for computational biases associated with allelic counts.
It comprises of the following steps:
1.	Local phasing of genetic data using PHASER
2.	Creation of parental genomes to align sequencing data to
3.	Re-allocation of multimapping reads using RSEM
4.	Selection of the best mapping for each read across the two parental genomes
5.	Outputs haplotype and site level allelic counts





## INSTALLATION AND RUNNING
### 1. Download nextflow

`curl -fsSL get.nextflow.io | bash`

Make sure you have [Java v8+](https://www.oracle.com/java/technologies/javase-downloads.html):

`java -version`

### 2. Install either [Docker]((https://docs.docker.com/get-docker/)) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) if cluster doesn't have them yet


### 3. Run PAC with following command:

* You can either run with this:

`path_to/nextflow run https://github.com/anna-saukkonen/PAC -r main --genome_version GRCh37/38 --reads "path_to_reads_{1,2}.fq.gz" --variants "path_to_variants" --id ID -profile docker/singularity`

-r command specifies the branch


* Or download repository and run with this:

`path_to/nextflow run PAC/main.nf --genome_version GRCh37/38 --reads "path_to_reads_{1,2}.fq.gz" --variants "path_to_variants" --id ID -profile docker/singularity`



## OPTIONS:

### Required
##### --genome_version: GRCh37 *or* GRCh38


##### --reads:  "pathtoreads_**{1,2}.fq.gz**
reads have to be saved in the same directory in the format: *path_to_read_1.fq.gz* and *path_to_read_2.fq.gz*


##### --variants:  "path_to_variants.vcf.gz"
vcf file needs to be phased

##### -profile:  docker *or* singularity
     

##### --id:  "name_of_sample"  
this needs to be same as in the VCF file      




### Optional
##### -N:  name@email_address.com  (To receive email when the pipeline is finished)

##### -outdir:  "name_of_results_file_directory"  
(default:  "/pac_results")
 
##### -cpus:  number  
(default:10  We recommend at least 10 for speed)

Depending on the size of file you might need up to 128000MB, min 64000MB



## OUTPUT

PAC generates 4 output files:

#### * Haplotype level ASE calls:
  - *ID*_gene_level_ae.txt
  
| Haplotype level ASE results columns  | Description          |
| ------------------------------------ |:--------------------:| 
| contig                               | chromosome           | 
| start                                | gene start position  |  
| stop                                 | gene end position    |
| name                                 | gene name            |
| aCount                               | haplotype a coverage |
| bCount                               | haplotype b coverage |
| totalCount                           | total coverage       |

#### * Single nucleotide level ASE calls from PAC: 
  - results_2genomes_*ID*.RSEM.STAR.SOFT.NOTRIM_baq.txt
  - results_2genomes_*ID*.RSEM.STAR.SOFT.NOTRIM.txt
   
#### * Single nucleotide level ASE calls based on standard single genome mapping (for comparison):
  - results_1genome_*ID*.SOFT.NOTRIM_baq.txt
  - results_1genome_*ID*.SOFT.NOTRIM.txt

| Single nucleotide level ASE results columns  | Description                 |
| -------------------------------------------- |:---------------------------:| 
| Chr                                          | chromosome                  | 
| Pos                                          | position along chromosome   |  
| RefAl                                        | reference allele            |
| AltAl                                        | alternative allele          |
| MapRef                                       | reference allele coverage   |
| MapAlt                                       | alternative allele coverage |
| MapRatio                                     | reference allele ratio      |
| Mapcov                                       | total coverage at the site  |





## TEST DATASET

To test PAC on smaller dataset:

load java

load singularity

git clone https://github.com/anna-saukkonen/PAC.git

path_to_nextflow/nextflow run PAC/main.nf --genome_version GRCh37 --reads "PAC/test/NA12890_merged_sample_0.005_{1,2}.fq.gz" --variants "PAC/test/NA12877_output.phased.downsampled.vcf.gz" --id NA12877 -profile singularity

See [this](https://github.com/anna-saukkonen/PAC/tree/main/test/pac_results) folder for output files you should get




```
 Just use
   __    ___      __
 ||__)  /___\\  /   `
 ||    /     \\ \\__, ,

 man ;)
 ```
