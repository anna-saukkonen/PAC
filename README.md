# PAC
```
 Just use
   __    ___      __
 ||__)  /___\\  /   `
 ||    /     \\ \\__, ,

 man ;)
 ```

**P**ersonalised **A**SE **C**aller

Author: Anna Saukkonen
anna.saukkonen@gmail.com

## Introduction:

This pipeline has been created to adjust for reference mapping biases that often cause error in the detection of allele specific expression.

It comprises of the following steps:

1. Local phasing of genetic data using PHASER
2. Creation of parental genomes to align sequencing data to
3. Selection of the best mapping for each read across the two parental genomes


## Quick Start
1. Download nextflow
2. Install either Docker or Singularity
3. Run PAC with following command:

`nextflow run anna-saukkonen/main.nf --genome_version GRCh37 --reads "pathtoreads_{1,2}.fq.gz" --variants "pathtovariants" -profile docker`


### Installation of dependancies
Make sure you have Java v8+:
`java -version`

Install Nextflow
`curl -fsSL get.nextflow.io | bash`

Download docker to your system from
[here](https://docs.docker.com/get-docker/)




## Options:

#### Required
--genome_version GRCh37/GRCh38

--reads "pathtoreads **_{1,2}.fq.gz**"

--variants "pathtovariants"

-profile docker/singularity


#### Additional non-essential

To be emailed when the pipeline is finished
-N name@email_address.com

Name of the sample
-id "name_of_sample"
Default: "default_id"

Name of output directory
-outdir "name_of_results_file_directory"
Default: "/pac_results"




