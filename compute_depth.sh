#!/bin/bash
faidx /repositorio/exomas/datos/datos_gatk/referencia/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -i  chromsizes >  chrom.sizes
awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' /repositorio/exomas/datos/datos_gatk/referencia/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai > bed_split.bed
samtools view -L bed_split.bed -o out.bam $1
cat chrom.sizes|sort -V > sizes.genome.sort
bedtools coverage -a /repositorio/exomas/datos/coverage/hg38/xgen-exome.bed  -b out.bam -g sizes.genome.sort -sorted -hist 
