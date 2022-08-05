# blobtools_tutorial

A practical guide to blobtools

## what is blobtools?

Blobtools was created by [Laetsch DR and Blaxter ML, 2017.](https://f1000research.com/articles/6-1287/v1) \
The current release is found [here.](https://zenodo.org/badge/latestdoi/23453/DRL/blobtools) \
This guide will follow "Workflow A" in the [manual](https://blobtools.readme.io/docs/what-is-blobtools), providing example code at each step.

Blobtools allows for taxonomic partitioning of metagenomes or metatranscriptomes, without requiring reference genomes.

## what input files does it require?

First, we will create a blobplot using just one set of trimmed paired reads: \

coral_1.fq.gz
coral_2.fq.gz

After the base name, 1 refers to forward reads, and 2 refers to reverse reads. \

Run Trinity to assemble these trimmed paired reads. The resulting metatranscriptome represents an octocoral holobiont, and is expected to contain genetic information both from the _Sympodium sp._ host and its endosymbiotic zooxanthellae.

Trinity_holobiont_metatranscriptome.fasta

## run bwa-mem to align the trimmed paired reads to the holobiont transcriptome

### index the metatranscriptome assembly
bwa index Trinity_holobiont_metatranscriptome.fasta

### map reads to the assembly as an array

#!/bin/bash
#PBS -V \
#PBS -N bwa_mem_array \
#PBS -q batch \
#PBS -S /bin/bash \
#PBS -l nodes=1:ppn=8 \
#PBS -l mem=10G \
#PBS -l walltime=48:00:00 \
#PBS -o /error_logs/bwa_mem_array.out \
#PBS -e /error_logs/bwa_mem_array.err \
#PBS -J 1-2 \

mkdir -p /results/alignments
cd /results/alignments

### first, align the reads
bwa mem -t 8 -HM -k 19 -w 100 -d 100 -R "@RG\tID:${PBS_ARRAY_INDEX}\tSM:${PBS_ARRAY_INDEX}\tLB:${PBS_ARRAY_INDEX}\tPL:ILLUMINA" \
/nas5/dpierro/data/assemblies/Trinity_holobiont_metatranscriptome.fasta \
/nas5/dpierro/data/trimmed_paired_reads_renamed/coral_${PBS_ARRAY_INDEX}_1.fq.gz \
/nas5/dpierro/data/trimmed_paired_reads_renamed/coral_${PBS_ARRAY_INDEX}_2.fq.gz \
> coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam &&

samtools view -@ 8 -bS -o coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam_bam \
coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam &&

samtools sort -m 1G -@ 8 coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam_bam -o coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam &&

samtools index coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam &&

rm *sam* &&

samtools flagstat coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam > coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam.flagstat

