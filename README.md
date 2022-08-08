# blobtools_tutorial

A practical guide to BlobTools

## why use BlobTools?

BlobTools allows taxonomic partitioning of metagenome or metatranscriptome assemblies, without requiring reference genomes. Using GC content, read coverage, and taxonomy, it can visualize how much genetic information in the assembly comes from a host species versus members of its microbiome. 

The example dataset here includes the metatranscriptome of an octocoral host and its associated *Symbiodinium* endosymbionts, during a heat stress experiment. 

BlobTools was created by [Laetsch DR and Blaxter ML, 2017.](https://f1000research.com/articles/6-1287/v1) The current release is found [here.](https://zenodo.org/badge/latestdoi/23453/DRL/blobtools)

This guide follows "Workflow A" in the [manual](https://blobtools.readme.io/docs/what-is-blobtools), providing example code and explanations at each step.

## what input? trimmed paired reads

After running Trimmomatic to trim the reads, rename the paired reads as:

After the base name, list an array number, followed by a 1 (forward read) or 2 (reverse read).

This example will use two sets of trimmed paired reads:
```
coral_1_1.fq.gz
coral_1_2.fq.gz

coral_2_1.fq.gz
coral_2_2.fq.gz
```
Record which array numbers correspond to which experimental conditions. Here, coral_1 experiences control conditions, while coral_2 experiences thermal stress. 

Then, run Trinity to assemble these trimmed paired reads, generating this assembly:
```
holobiont.fasta
```
The resulting holobiont metatranscriptome assembly is expected to contain genetic information from the coral host and its endosymbionts. 

We will use PBS scripts to allow for PBS arrays. Many of these commands can also be run as loops. Unless otherwise indicated, all PBS scripts below will have these flags:
```
#!/bin/bash
#PBS -V
#PBS -N job_name
#PBS -q batch
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=10G
#PBS -l walltime=48:00:00
#PBS -o /error_logs/job_name.out
#PBS -e /error_logs/job_name.err
#PBS -J 1-2
```

## use bwa-mem to align trimmed paired reads to the assembly

### index the assembly
```
bwa index holobiont.fasta
```
### read mapping
Map reads as an array in a PBS script that includes the following steps:

create an output folder
```
mkdir -p /alignments
cd /alignments
```
align the reads
```
bwa mem -t 8 -HM -k 19 -w 100 -d 100 -R "@RG\tID:${PBS_ARRAY_INDEX}\tSM:${PBS_ARRAY_INDEX}\tLB:${PBS_ARRAY_INDEX}\tPL:ILLUMINA"
holobiont.fasta
coral_${PBS_ARRAY_INDEX}_1.fq.gz
coral_${PBS_ARRAY_INDEX}_2.fq.gz
> coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam &&
```
convert sam to bam
```
samtools view -@ 8 -bS -o coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam_bam \
coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam &&
```
sort the bam files
```
samtools sort -m 1G -@ 8 coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.sam_bam 
-o coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam &&
```
index the bam files
```
samtools index coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam &&
```
remove all intermediate output files
```
rm *sam* &&
```
get statistics on the alignment
```
samtools flagstat coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam 
> coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam.flagstat
```
## install BlobTools

It can be helpful to create a conda environment for blobtools.
