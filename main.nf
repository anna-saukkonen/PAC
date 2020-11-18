#!/usr/bin/env nextflow


/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */

params.variants   = "/away/asaukkonen/gitdir/NA12877_output.phased.vcf.gz"
params.reads       = "/away/asaukkonen/gitdir/NA12890_merged_{1,2}.fq.gz" 

params.genome        = params.genomes[ params.genome_version ]?.genome
params.annot         = params.genomes[ params.genome_version ]?.annot



// Check if genome exists in the config file
if (params.genomes && params.genome_version && !params.genomes.containsKey(params.genome_version)) {
    exit 1, "The provided genome '${params.genome_version}' is not available. Currently the available genomes are ${params.genomes.keySet().join(", ")}. Please check your spelling."
}

if (!params.variants) exit 1, "Path to phased variants has to be specified!"
if (!params.reads) exit 1, "Path to reads has to be specified!"



Channel
  .fromFilePairs(params.reads)
  .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\n"}
  .into {reads_ch1; reads_ch2; reads_ch3}




log.info """\

nextflow run anna-saukkonen/main.nf --genome_version GRCh37 --reads '*_{1,2}.fq.gz' -profile docker


genome        : $params.genome
reads         : $params.reads
variants      : $params.variants
annot         : $params.annot
"""


readLength = ''

process read_length {
  input:
    path variants from params.variants 

  output:
    val readLength into read_len_ch

  script:

  """
  $readLength = \$(gunzip -c ${variants} | sed '2q;d' | wc -m)
  """
}


process prepare_star_genome_index {
  tag "$genome.baseName"


  input:
    path genome from params.genome
    path annot from params.annot
    val readLength from read_len_ch
  output:
    path STARhaploid into genome_dir_ch

  script:

  """
  mkdir STARhaploid



  STAR --runMode genomeGenerate \
       --genomeDir STARhaploid \
       --genomeFastaFiles ${genome} \
       --sjdbGTFfile ${annot} \
       --sjdbOverhang ${readLength} \
       --runThreadN ${task.cpus}
  """
}






process rnaseq_mapping_star {
  tag "$id"

  
  input: 
    path genome from params.genome 
    path STARhaploid from genome_dir_ch
    set val(id), file(reads) from reads_ch1
    val readLength from read_len_ch

  output: 
    tuple \
      val(id), \
      path('NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam'), \
      path('NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam.bai') into aligned_bam_ch

  script: 

  """
  # Align reads to genome
  STAR --genomeDir STARhaploid \
       --readFilesIn $reads \
       --readFilesCommand zcat \
       --runThreadN ${task.cpus} \
       --outSAMstrandField intronMotif \
       --outFilterMultimapNmax 30 \
       --alignIntronMax 1000000 \
       --alignMatesGapMax 1000000 \
       --outMultimapperOrder Random \
       --outSAMunmapped Within \
       --outSAMattrIHstart 0 \
       --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
       --sjdbOverhang ${readLength} \
       --outFilterMismatchNmax 8 \
       --outSAMattributes NH nM NM MD HI \
       --outSAMattrRGline  ID:$id PU:Illumina PL:Illumina LB:NA12877.SOFT.NOTRIM SM:NA12877.SOFT.NOTRIM CN:Seq_centre \
       --outSAMtype BAM SortedByCoordinate \
       --twopassMode Basic \
       --outFileNamePrefix NA12877.SOFT.NOTRIM.STAR.pass2. \
       --outSAMprimaryFlag AllBestScore

  # Index the BAM file
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam
  """
}











process clean_up_reads {

  
  input:
    tuple val(id), path(bam), path(index) from aligned_bam_ch
    path variants from params.variants

  output:
    path ('STAR_original/phaser_version.bam') into phaser_ch
    path ('STAR_original/phaser_version.bam.bai') into phaser_bai_ch
    path ('STAR_original/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') into pp_um_ch

  script:

  """
  mkdir STAR_original
  

  #KEEP ONLY PROPERLY PAIRED READS
  samtools view -@ 10 -f 0x0002 -b -o NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam

  #KEEP UNIQUELY MAPPED READS
  samtools view -h NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam | grep -P "NH:i:1\t|^@" | samtools view -bS - > NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

  #Create BAM compatible with PHASER:
  gunzip -c ${variants} | grep -q 'chr' || (samtools view -h NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | sed -e 's/chr//' >> phaser_version.sam; samtools view -bh phaser_version.sam >> phaser_version.bam; samtools index phaser_version.bam; rm phaser_version.sam)
  gunzip -c  ${variants}  | grep -q 'chr' && (samtools view -bh NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam >> phaser_version.bam; samtools index phaser_version.bam)

  mv phaser_version.bam STAR_original/phaser_version.bam
  mv phaser_version.bam.bai STAR_original/phaser_version.bam.bai
  mv NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam STAR_original/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam
  """
}



process phaser_step {

  input:
  path variants from params.variants
  path ('phaser_version.bam') from phaser_ch
  path ('phaser_version.bam.bai') from phaser_bai_ch

  output:
  path ('NA12877_output_phaser.vcf') into (phaser_out_ch1, phaser_out_ch2, phaser_out_ch3)

  script:

  """
  tabix -p vcf ${variants}


  python2 /phaser/phaser/phaser.py --vcf ${variants} --bam phaser_version.bam --paired_end 1 --mapq 0 --baseq 10 --isize 0 --include_indels 1 --sample NA12877 --id_separator + --pass_only 0 --o NA12877_output_phaser

  gunzip NA12877_output_phaser.vcf.gz
  rm phaser_version.bam
  rm phaser_version.bam.bai
  """
}





process create_parental_genomes {


  input:
    path genome from params.genome
    path annot from params.annot
    path ('NA12877_output_phaser.vcf') from phaser_out_ch1

  output:
    path ('STAR_2Gen_Ref/maternal.chain') into maternal_chain_ch
    path ('STAR_2Gen_Ref/paternal.chain') into paternal_chain_ch
    path ('STAR_2Gen_Ref/chr22_NA12877.map') into map_ch
    path ('STAR_2Gen_Ref/NA12877_maternal.fa') into (mat_fa1, mat_fa2)
    path ('STAR_2Gen_Ref/NA12877_paternal.fa') into (pat_fa1, pat_fa2)
    path ('STAR_2Gen_Ref/mat_annotation.gtf') into (mat_annotation_ch1, mat_annotation_ch2)
    path ('STAR_2Gen_Ref/not_lifted_m.txt') into not_lift_m_ch
    path ('STAR_2Gen_Ref/pat_annotation.gtf') into (pat_annotation_ch1, pat_annotation_ch2)
    path ('STAR_2Gen_Ref/not_lifted_p.txt') into not_lift_p_ch
    path ('STAR_2Gen_Ref/map_over.txt') into (adjusted_ref_ch1, adjusted_ref_ch2)

  
  script:

  """
  mkdir STAR_2Gen_Ref
  java -Xmx10000m -jar /vcf2diploid_v0.2.6a/vcf2diploid.jar -id NA12877 -chr ${genome} -vcf NA12877_output_phaser.vcf -outDir STAR_2Gen_Ref > logfile.txt
  
  
  liftOver -gff ${annot} STAR_2Gen_Ref/maternal.chain STAR_2Gen_Ref/mat_annotation.gtf STAR_2Gen_Ref/not_lifted_m.txt
  liftOver -gff ${annot} STAR_2Gen_Ref/paternal.chain STAR_2Gen_Ref/pat_annotation.gtf STAR_2Gen_Ref/not_lifted_p.txt
  
  
  
  cat STAR_2Gen_Ref/chr1_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr2_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr3_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr4_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr5_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr6_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr7_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr8_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr9_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr10_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr11_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr12_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr13_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr14_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr15_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr16_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr17_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr18_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr19_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr20_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr21_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chr22_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chrX_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chrY_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  cat STAR_2Gen_Ref/chrM_NA12877_maternal.fa >> STAR_2Gen_Ref/NA12877_maternal.fa
  
  

  
  cat STAR_2Gen_Ref/chr1_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr2_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr3_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr4_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr5_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr6_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr7_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr8_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr9_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr10_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr11_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr12_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr13_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr14_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr15_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr16_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr17_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr18_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr19_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr20_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr21_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chr22_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chrX_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chrY_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  cat STAR_2Gen_Ref/chrM_NA12877_paternal.fa >> STAR_2Gen_Ref/NA12877_paternal.fa
  

  sed 's/\\*/N/g' STAR_2Gen_Ref/NA12877_maternal.fa > STAR_2Gen_Ref/NA12877_maternal.hold.fa
  mv STAR_2Gen_Ref/NA12877_maternal.hold.fa STAR_2Gen_Ref/NA12877_maternal.fa
  
  sed 's/\\*/N/g' STAR_2Gen_Ref/NA12877_paternal.fa > STAR_2Gen_Ref/NA12877_paternal.hold.fa
  mv STAR_2Gen_Ref/NA12877_paternal.hold.fa STAR_2Gen_Ref/NA12877_paternal.fa

  mv NA12877_output_phaser.vcf STAR_2Gen_Ref/NA12877_output_phaser.vcf
  cd STAR_2Gen_Ref/
  
  perl ${baseDir}/bin/adjust_reference.pl NA12877_output_phaser.vcf NA12877
  """
} 



process STAR_reference_genomes {
  input:
    path ('STAR_2Gen_Ref/NA12877_maternal.fa') from mat_fa1
      path ('STAR_2Gen_Ref/NA12877_paternal.fa') from pat_fa1
      path ('STAR_2Gen_Ref/mat_annotation.gtf') from mat_annotation_ch1
      path ('STAR_2Gen_Ref/pat_annotation.gtf') from pat_annotation_ch1
      val readLength from read_len_ch

  output:
    path Paternal_STAR into Paternal_STAR_ch
    path Maternal_STAR into Maternal_STAR_ch
    

  script:

  """
  mkdir Maternal_STAR
  mkdir Paternal_STAR

  STAR --runMode genomeGenerate --genomeDir Paternal_STAR --genomeFastaFiles STAR_2Gen_Ref/NA12877_paternal.fa --sjdbGTFfile STAR_2Gen_Ref/pat_annotation.gtf --sjdbOverhang ${readLength} --runThreadN 5 --outTmpDir pat
  STAR --runMode genomeGenerate --genomeDir Maternal_STAR --genomeFastaFiles STAR_2Gen_Ref/NA12877_maternal.fa --sjdbGTFfile STAR_2Gen_Ref/mat_annotation.gtf --sjdbOverhang ${readLength} --runThreadN 5 --outTmpDir mat
  """    

}


process map_paternal_gen_filter {
  tag "$id"
  
  input:
    path Paternal_STAR from Paternal_STAR_ch
    set val(id), file(reads) from reads_ch2
    path ('STAR_2Gen_Ref/pat_annotation.gtf') from pat_annotation_ch2
    path ('STAR_2Gen_Ref/NA12877_paternal.fa') from pat_fa2
    val readLength from read_len_ch

  output:
    path ('STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') into (paternal_mapgen_ch1, paternal_mapgen_ch2, paternal_mapgen_ch3)
    path ('STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam') into pat_rsem_ch

  script:

  """
  STAR --genomeDir Paternal_STAR --runThreadN 10 --quantMode TranscriptomeSAM --readFilesIn $reads --readFilesCommand zcat --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang ${readLength} --outFilterMismatchNmax 8 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:NA12877.SOFT.NOTRIM PU:Illumina PL:Illumina LB:NA12877.SOFT.NOTRIM SM:NA12877.SOFT.NOTRIM CN:Seq_centre --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix NA12877.SOFT.NOTRIM.STAR.pass2. --outSAMprimaryFlag AllBestScore

  
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam

  #KEEP ONLY PROPERLY PAIRED READS
  samtools view -@ 10 -f 0x0002 -b -o NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam

  #KEEP UNIQUELY MAPPED READS
  samtools view -h NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam | grep -P "NH:i:1\t|^@" | samtools view -bS - > NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

  mkdir STAR_Paternal
  mv NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

  ##Create RSEM Files:
  mkdir RSEM_MAT_GEN

  /RSEM/rsem-prepare-reference --gtf STAR_2Gen_Ref/pat_annotation.gtf STAR_2Gen_Ref/NA12877_paternal.fa RSEM_MAT_GEN/RSEM_MAT_GEN

  /RSEM/rsem-calculate-expression --bam --output-genome-bam --sampling-for-bam -p 20 --paired-end  NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.toTranscriptome.out.bam RSEM_MAT_GEN/RSEM_MAT_GEN NA12877.RSEM.TEST

  samtools view -@ 10 -f 0x0002 -b -o NA12877.RSEM.TEST.genome.PP.bam NA12877.RSEM.TEST.genome.bam
  samtools sort -@ 20 -o NA12877.RSEM.TEST.genome.PP.s.bam NA12877.RSEM.TEST.genome.PP.bam
  mv NA12877.RSEM.TEST.genome.PP.s.bam NA12877.RSEM.TEST.genome.PP.bam
  samtools index NA12877.RSEM.TEST.genome.PP.bam
  samtools view -h NA12877.RSEM.TEST.genome.PP.bam | grep -P "ZW:f:1|^@" | samtools view -bS - > NA12877.RSEM.TEST.genome.PP.SM.bam
  samtools index NA12877.RSEM.TEST.genome.PP.SM.bam
  mv NA12877.RSEM.TEST.genome.PP.SM.bam STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam
  """

}


process map_maternal_gen_filter {
  tag "$id"

  input:
    path Maternal_STAR from Maternal_STAR_ch
    set val(id), file(reads) from reads_ch3
    path ('STAR_2Gen_Ref/mat_annotation.gtf') from mat_annotation_ch2
    path ('STAR_2Gen_Ref/NA12877_maternal.fa') from mat_fa2
    val readLength from read_len_ch

  output:
    path ('STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') into (maternal_mapgen_ch1, maternal_mapgen_ch2, maternal_mapgen_ch3) 
    path ('STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam') into mat_rsem_ch

  script:

  """
  STAR --genomeDir Maternal_STAR --runThreadN 10 --quantMode TranscriptomeSAM --readFilesIn $reads --readFilesCommand zcat --outSAMstrandField intronMotif --outFilterMultimapNmax 30 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outMultimapperOrder Random --outSAMunmapped Within --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --sjdbOverhang ${readLength} --outFilterMismatchNmax 8 --outSAMattributes NH nM NM MD HI --outSAMattrRGline  ID:NA12877.SOFT.NOTRIM PU:Illumina PL:Illumina LB:NA12877.SOFT.NOTRIM SM:NA12877.SOFT.NOTRIM CN:Seq_centre --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix NA12877.SOFT.NOTRIM.STAR.pass2. --outSAMprimaryFlag AllBestScore


  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam

  #KEEP ONLY PROPERLY PAIRED READS
  samtools view -@ 10 -f 0x0002 -b -o NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam

  #KEEP UNIQUELY MAPPED READS
  samtools view -h NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.bam | grep -P "NH:i:1\t|^@" | samtools view -bS - > NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam
  samtools index NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

  mkdir STAR_Maternal
  mv NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam

  ##Create RSEM Files:
  mkdir RSEM_MAT_GEN

  /RSEM/rsem-prepare-reference --gtf STAR_2Gen_Ref/mat_annotation.gtf STAR_2Gen_Ref/NA12877_maternal.fa RSEM_MAT_GEN/RSEM_MAT_GEN
  /RSEM/rsem-calculate-expression --bam --output-genome-bam --sampling-for-bam -p 20 --paired-end  NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.toTranscriptome.out.bam RSEM_MAT_GEN/RSEM_MAT_GEN NA12877.RSEM.TEST
  
  samtools view -@ 10 -f 0x0002 -b -o NA12877.RSEM.TEST.genome.PP.bam NA12877.RSEM.TEST.genome.bam
  samtools sort -@ 20 -o NA12877.RSEM.TEST.genome.PP.s.bam NA12877.RSEM.TEST.genome.PP.bam
  mv NA12877.RSEM.TEST.genome.PP.s.bam NA12877.RSEM.TEST.genome.PP.bam
  samtools view -h NA12877.RSEM.TEST.genome.PP.bam | grep -P "ZW:f:1|^@" | samtools view -bS - > NA12877.RSEM.TEST.genome.PP.SM.bam
  samtools index NA12877.RSEM.TEST.genome.PP.SM.bam
  mv NA12877.RSEM.TEST.genome.PP.SM.bam STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam
  """

}


process merge_parental_bam {
  input:
    path ('STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from maternal_mapgen_ch1
    path ('STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from paternal_mapgen_ch1
    path ('STAR_2Gen_Ref/map_over.txt') from adjusted_ref_ch1
    path ('NA12877_output_phaser.vcf') from phaser_out_ch2
    path ('STAR_original/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from pp_um_ch
  output:
    path ('results/results*.txt') into results_ch
  script:

  """
  samtools view STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | cut -f1 | sort | uniq >> maternal_tags.txt
  samtools view STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | cut -f1 | sort | uniq >> paternal_tags.txt
  cat maternal_tags.txt paternal_tags.txt | sort | uniq -u >> unique_tags.txt
  cat maternal_tags.txt paternal_tags.txt | sort | uniq -d >> duplicate_tags.txt
  samtools view STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | grep -Fwf duplicate_tags.txt >> tempout_mat.sam
  samtools view STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | grep -Fwf duplicate_tags.txt >> tempout_pat.sam
  sort -k 1,1 tempout_mat.sam > tempout_mat.sort.sam
  sort -k 1,1 tempout_pat.sam > tempout_pat.sort.sam

  perl ${baseDir}/bin/filter_2genomes.pl tempout_mat.sort.sam tempout_pat.sort.sam

  cat maternal_wins.txt unique_tags.txt > maternal_wins_final.txt
  cat paternal_wins.txt unique_tags.txt > paternal_wins_final.txt

  samtools view -H STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam > final_mat.sam
  samtools view -H STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam > final_pat.sam
  samtools view STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | grep -Fwf maternal_wins_final.txt >> final_mat.sam
  samtools view STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | grep -Fwf paternal_wins_final.txt >> final_pat.sam
  samtools view -bS final_mat.sam -o final_mat.bam
  samtools index final_mat.bam
  samtools view -bS final_pat.sam -o final_pat.bam
  samtools index final_pat.bam

  perl ${baseDir}/bin/compare_2genomes.pl STAR_2Gen_Ref/map_over.txt NA12877_output_phaser.vcf final_mat.bam final_pat.bam NA12877 results_2genomes_NA12877.SOFT.NOTRIM_baq.txt results_2genomes_NA12877.SOFT.NOTRIM.txt

  perl ${baseDir}/bin/compare_basic_map.pl NA12877_output_phaser.vcf STAR_original/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam NA12877 results_1genome_NA12877.SOFT.NOTRIM_baq.txt results_1genome_NA12877.SOFT.NOTRIM.txt

  mkdir results
  mv results_2genomes_NA12877.SOFT.NOTRIM_baq.txt results_2genomes_NA12877.SOFT.NOTRIM.txt results_1genome_NA12877.SOFT.NOTRIM_baq.txt results_1genome_NA12877.SOFT.NOTRIM.txt results/
  """

}



process extra_reads_rsem {

  input:
    path ('STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from maternal_mapgen_ch2
    path ('STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam') from mat_rsem_ch
    path ('STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from paternal_mapgen_ch2
    path ('STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam') from pat_rsem_ch

  output:
    path ('Maternal.RSEM.bam') into mat_rsembam
    path ('Paternal.RSEM.bam') into pat_rsembam


  script:

  """
  samtools view STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | cut -f1 | sort | uniq >> maternal_tags_UM.txt
  samtools view STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam | cut -f1 | sort | uniq > maternal_tags_UM.RSEM.txt

  perl ${baseDir}/bin/filter_rsem.pl maternal

  samtools view STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam | cut -f1 | sort | uniq >> paternal_tags_UM.txt
  samtools view STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam | cut -f1 | sort | uniq > paternal_tags_UM.RSEM.txt

  perl ${baseDir}/bin/filter_rsem.pl paternal


  samtools view -H STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam > Maternal.RSEM.sam
  samtools view STAR_Maternal/NA12877.RSEM.TEST.genome.PP.SM.bam | grep -Fwf extra.rsem.maternal.txt | sed -e 's/339\tchr/83\tchr/' | sed -e 's/355\tchr/99\tchr/' | sed -e 's/403\tchr/147\tchr/' | sed -e 's/419\tchr/163\tchr/' >> Maternal.RSEM.sam
  samtools view -bS Maternal.RSEM.sam -o Maternal.RSEM.bam

  samtools view -H STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam > Paternal.RSEM.sam
  samtools view STAR_Paternal/NA12877.RSEM.TEST.genome.PP.SM.bam | grep -Fwf extra.rsem.paternal.txt | sed -e 's/339\tchr/83\tchr/' | sed -e 's/355\tchr/99\tchr/' | sed -e 's/403\tchr/147\tchr/' | sed -e 's/419\tchr/163\tchr/' >> Paternal.RSEM.sam
  samtools view -bS Paternal.RSEM.sam -o Paternal.RSEM.bam
  """
}



process add_rsemreads_bam {
  publishDir "$params.outdir/result", mode: 'copy'

  input:
    path ('Maternal.RSEM.bam') from mat_rsembam
    path ('Paternal.RSEM.bam') from pat_rsembam
    path ('STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from paternal_mapgen_ch3
    path ('STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam') from maternal_mapgen_ch3
    path ('STAR_2Gen_Ref/map_over.txt') from adjusted_ref_ch2
    path ('NA12877_output_phaser.vcf') from phaser_out_ch3

  output:
    path ('results_2genomes_NA12877.RSEM.STAR.SOFT.NOTRIM.txt')
    
  script:

  """
  samtools merge Maternal.RSEM.STAR.bam STAR_Maternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam Maternal.RSEM.bam
  samtools merge Paternal.RSEM.STAR.bam STAR_Paternal/NA12877.SOFT.NOTRIM.STAR.pass2.Aligned.sortedByCoord.out.PP.UM.bam Paternal.RSEM.bam
  samtools view Maternal.RSEM.STAR.bam | cut -f1 | sort | uniq >> maternal_tags.txt
  samtools view Paternal.RSEM.STAR.bam | cut -f1 | sort | uniq >> paternal_tags.txt
  cat maternal_tags.txt paternal_tags.txt | sort | uniq -u >> unique_tags.txt
  cat maternal_tags.txt paternal_tags.txt | sort | uniq -d >> duplicate_tags.txt
  samtools view Maternal.RSEM.STAR.bam | grep -Fwf duplicate_tags.txt >> tempout_mat.sam
  samtools view Paternal.RSEM.STAR.bam | grep -Fwf duplicate_tags.txt >> tempout_pat.sam
  sort -k 1,1 tempout_mat.sam > tempout_mat.sort.sam
  sort -k 1,1 tempout_pat.sam > tempout_pat.sort.sam

  perl ${baseDir}/bin/filter_2genomes.pl tempout_mat.sort.sam tempout_pat.sort.sam

  cat maternal_wins.txt unique_tags.txt > maternal_wins_final.txt
  cat paternal_wins.txt unique_tags.txt > paternal_wins_final.txt
  samtools view -H Maternal.RSEM.STAR.bam > final_mat.sam
  samtools view -H Paternal.RSEM.STAR.bam > final_pat.sam
  samtools view Maternal.RSEM.STAR.bam | grep -Fwf maternal_wins_final.txt >> final_mat.sam
  samtools view Paternal.RSEM.STAR.bam | grep -Fwf paternal_wins_final.txt >> final_pat.sam
  samtools view -bS final_mat.sam -o final_mat.bam
  samtools sort -@ 20 -o final_mat.sorted.bam final_mat.bam
  samtools index final_mat.sorted.bam
  samtools view -bS final_pat.sam -o final_pat.bam
  samtools sort -@ 20 -o final_pat.sorted.bam final_pat.bam
  samtools index final_pat.sorted.bam

  perl ${baseDir}/bin/compare_2genomes.pl STAR_2Gen_Ref/map_over.txt NA12877_output_phaser.vcf final_mat.sorted.bam final_pat.sorted.bam NA12877 results_2genomes_NA12877.RSEM.STAR.SOFT.NOTRIM_baq.txt results_2genomes_NA12877.RSEM.STAR.SOFT.NOTRIM.txt



  """

}



