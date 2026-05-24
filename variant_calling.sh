#!/bin/bash


# Create directories
mkdir -p Sequences Reference QC Alignment Variants

# Step 1: Download sample sequences

cd Sequences

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/034/SRR22352834/SRR22352834_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR223/034/SRR22352834/SRR22352834_2.fastq.gz
 
# Step 2: Download chr11 Reference Genome

cd Reference

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chr11.fa.gz

gunzip chr11.fa.gz

cd ..

# Step 3: Reference Genome Indexing

bwa index Reference/chr11.fa

samtools faidx Reference/chr11.fa

# Step 4: Quality Control

fastqc SRR22352834_1.fastq.gz SRR22352834_2.fastq.gz -o QC

# Step 5: Read Alignment

bwa mem -t 4 Reference/chr11.fa SRR22352834_1.fastq.gz SRR22352834_2.fastq.gz > Alignment/SRR22352834.sam

# Step 6: SAM to BAM Conversion

samtools view -Sb Alignment/SRR22352834.sam > Alignment/SRR22352834.bam

# Step 7: Sort BAM

samtools sort Alignment/SRR22352834.bam -o Alignment/SRR22352834.sorted.bam

# Step 8: Index BAM

samtools index Alignment/SRR22352834.sorted.bam

# Step 9: Alignment Statistics

samtools flagstat Alignment/SRR22352834.sorted.bam > Alignment/alignment_stats.txt

# Step 10: Variant Calling
# Generate bcf file
 
bcftools mpileup -f Reference/chr11.fa Alignment/SRR22352834.sorted.bam -Ou -o Variants/output.bcf

# Generate vcf from bcf file
bcftools call -mv -Oz -o Variants/variants.vcf Variants/output.bcf

# Step 11: Variant Filtering

bcftools filter -i 'QUAL>30 && DP>10' Variants/variants.vcf -Oz -o Variants/filtered_variants.vcf

# Filter SNPs only

bcftools view -v snps Variants/filtered_variants.vcf -Oz -o Variants/filtered_snps.vcf

# Filter INDELs only

bcftools view -v indels Variants/filtered_variants.vcf -Oz -o Variants/filtered_indels.vcf

# Step 12: Variant Annotation

# Download database once:
# snpEff download GRCm39.115

java -Xmx8g -jar snpEff.jar GRCm39.115 Variants/filtered_variants.vcf.gz > Variants/annotated_variants.vcf

# Step 13: Extract HIGH Impact Variants

grep "HIGH" Variants/annotated_variants.vcf > Variants/high_impact_variants.vcf

# Pipeline Completed

echo "Variant calling pipeline completed successfully."