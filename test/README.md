## Smaller simultaed dataset from our [paper](https://www.biorxiv.org/content/10.1101/2021.07.13.452202v1) to test pipeline


#### To run:

load java

load singularity

nextflow run run https://github.com/anna-saukkonen/PAC -r main --genome_version GRCh37 --reads "/test/NA12890_merged_sample_0.05_{1,2}.fq.gz" --variants "/test/NA12877_output.phased.downsampled.vcf.gz" --id NA12877 -profile singularity

In this folder you can see the results files you should receive
