# PAC **P**ersonalised **A**SE **C**aller
Author: Anna Saukkonen
anna.saukkonen@gmail.com

## Introduction:

This pipeline has been created to adjust for computational biases associated with allelic counts.
It comprises of the following steps:
1.	Local phasing of genetic data using PHASER
2.	Creation of parental genomes to align sequencing data to
3.	Selection of the best mapping for each read across the two parental genomes
4.	Re-allocation of multimapping reads using RSEM
5.	Outputs haplotype and site level allelic counts

See our [preprint](https://www.biorxiv.org/content/10.1101/2021.07.13.452202v1) for additional information



## Installation
##### 1. Download nextflow

`curl -fsSL get.nextflow.io | bash`

Make sure you have [Java v8+](https://www.oracle.com/java/technologies/javase-downloads.html):

`java -version`


##### 2. Install either [Docker]((https://docs.docker.com/get-docker/)) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) if cluster doesn't have them yet


##### 3. Run PAC with following command:

`nextflow run http://github.com/anna-saukkonen/PAC.git --genome_version GRCh37 --reads "path_to_reads_{1,2}.fq.gz" --variants "path_to_variants" -profile docker`




## Options:

#### Required
--genome_version:  GRCh37 *or* GRCh38


--reads:  "pathtoreads_**{1,2}.fq.gz**
reads have to be saved in the same directory in the format: *path_to_read_1.fq.gz* and *path_to_read_2.fq.gz*


--variants:  "path_to_variants"


-profile:  docker *or* singularity
     

-id:  "name_of_sample"       



#### Additional non-essential
To receive email when the pipeline is finished
-N:  name@email_address.com


Name of output directory
-outdir:  "name_of_results_file_directory"

default:  "/pac_results"
 

Number of threads
-cpus:  number

default:1

We recommend at least 10 for speed



```
 Just use
   __    ___      __
 ||__)  /___\\  /   `
 ||    /     \\ \\__, ,

 man ;)
 ```