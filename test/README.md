## Smaller simultaed dataset from our [paper](https://www.biorxiv.org/content/10.1101/2021.07.13.452202v1) to test pipeline


#### To run:

load java

load singularity

git clone https://github.com/anna-saukkonen/PAC.git


path_to_nextflow/nextflow run PAC/main.nf --genome_version GRCh37 --reads "PAC/test/NA12890_merged_sample_0.005_{1,2}.fq.gz" --variants "PAC/test/NA12877_output.phased.downsampled.vcf.gz" --id NA12877 -profile singularity

See [pac_results](https://github.com/anna-saukkonen/PAC/tree/main/test/pac_results) folder for the results files you should receive
