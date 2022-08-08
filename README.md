# blobtools_tutorial

A practical guide to BlobTools

## why use BlobTools?

BlobTools allows taxonomic partitioning of metagenome or metatranscriptome assemblies, without requiring reference genomes. Using GC content, read coverage, and taxonomy, it can visualize how much genetic information in the assembly comes from a host species versus members of its microbiome. 

The example dataset here includes the metatranscriptome of an octocoral host and its associated *Symbiodinium* endosymbionts, during a heat stress experiment. 

BlobTools was created by [Laetsch DR and Blaxter ML, 2017.](https://f1000research.com/articles/6-1287/v1) The current release is found [here.](https://zenodo.org/badge/latestdoi/23453/DRL/blobtools)

This guide follows "Workflow A" in the [manual](https://blobtools.readme.io/docs/what-is-blobtools), providing example code and explanations at each step.

## what input? 
### trimmed paired reads

After running Trimmomatic to trim the reads, rename the paired reads as follows.

First a base name, then an array number, then 1 (forward read) or 2 (reverse read).

This example will use two sets of trimmed paired reads:
```
coral_1_1.fq.gz
coral_1_2.fq.gz

coral_2_1.fq.gz
coral_2_2.fq.gz
```
Record which array numbers correspond to which experimental conditions. Here, coral_1 experiences control conditions, while coral_2 experiences thermal stress. 

Then, run Trinity to assemble these trimmed paired reads, generating the assembly ```holobiont.fasta```

This holobiont metatranscriptome assembly is expected to contain genetic information from the coral host and its endosymbionts. 

This tutorial will use PBS scripts to allow for PBS arrays. Many of these commands can also be run as loops. Unless otherwise indicated, all PBS scripts below will have these flags:
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

## run bwa-mem
### to align trimmed paired reads to the assembly

### index the assembly
```
bwa index holobiont.fasta
```
### read mapping
map reads as an array in a PBS script that includes the following steps:

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

Several methods for installing BlobTools can be found [here](https://github.com/DRL/blobtools)

Upload ```blob.txt``` and use it to create a conda environment
```
conda create --name blobtools --file blob.txt
```
Install the dependencies
```
conda activate blobtools
conda install -n blobtools matplotlib 
conda install -n blobtools docopt 
conda install -n blobtools tqdm 
conda install -n blobtools wget 
conda install -n blobtools pyyaml 
conda install -n blobtools git
conda install -n blobtools pysam --update-deps 
conda install -n blobtools blobtools
```
you should recieve a message that all requested packages are installed.

note where the blobtools environment is located, for instance:
```/envs/blobtools```

cd into the data directory, for instance:
```/envs/blobtools/lib/python3.9/site-packages/data```

then get the taxdump
```
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz
```
create nodesDB
```
blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
```
the output file nodesDB.txt should be in the data/ directory

move the /data directory and all of its contents (including nodesDB) into the directory /blobtools
```
cp -r /home/conda/envs/blobtools/lib/python3.9/site-packages/data/ /home/conda/envs/blobtools/
```
## coverage
### generate coverage files from mapping individual bwamem alignments to the holobiont reference

```
cd /alignments/

for i in coral_*_toHolobiont_bwamem.bam;
do blobtools map2cov --infile holobiont.fasta
--bam $i --output $i; done
```
## taxonomic annotation
### blastn holobiont contigs to the nt database to find taxonomic hits

split the reference assembly ```holobiont.fasta``` into 10+ fasta files, using the script ```splitFASTA.pl```

make the directory /split_holobiont

```
perl splitFASTA.pl -i holobiont.fasta -o /split_holobiont/ -s 10000
```
Run blastn as an array job. 

Download the blast nt database locally (shown here) or if possible, use a remote search

Troubleshoot:\
Specify ```#PBS -J 1-n``` such that n = the number of files that splitFASTA split the assembly into.\
Make sure ```#PBS -l nodes=1:ppn=8``` matches the number of threads in the blastn command ```-num_threads 8``` \
Ensure that blast is updated to the most recent version with ```conda update blast```

```
cd /fullpath/split_holobiont/

blastn -query holobiont-${PBS_ARRAY_INDEX}.fasta
-db /database/nt.oct.19.2021/nt -outfmt '6 qseqid staxids bitscore std'
-max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 8
-out ${PBS_ARRAY_INDEX}_to_nt.tab
```
In the split_holobiont directory, combine all blastn output files to one hits directory
```
cat *_to_nt.tab > holobiont_to_nt.tab
```
## create a BlobTools database

run as an array
make the output directory /blob_db

```
source ~/.bash_profile
conda activate blobtools

cd /fullpath/blob_db

blobtools create --infile holobiont.fasta
--hitsfile /fullpath/split_holobiont/holobiont_to_nt.tab
--cov /fullpath/alignments/coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam.cov
--bam /fullpath/alignments/coral_${PBS_ARRAY_INDEX}_toHolobiont_bwamem.bam
--out coral_${PBS_ARRAY_INDEX}_toHolobiont_blob.DB

```
the output databases are named coral_{array number}_toHolobiont_blob.DB.blobDB.json

## create a summary table
change {1..n} depending on how many array numbers you have in total

```
for i in {1..4}; do blobtools view -i /fullpath/blob_db/coral_${i}_toHolobiont_blob.DB.blobDB.json 
-o coral_$i -x bestsum -r phylum; done
```

## create a blobplot

create a blobplot of one individual file (for instance, array number 1)
```
blobtools plot -i coral_1_toHolobiont_blob.DB.blobDB.json --nohit --rank phylum --taxrule bestsum --format pdf --out coral1.phylum
```
create a blobplot for all files. 
change {1..n} depending on how many array numbers you have in total
since this may take several hours, running it in a screen is highly recommended
```
for i in {1..4}; do blobtools plot -i /home/dpierro/nas5/results/blob_db/coral_${i}_toHolobiont_blob.DB.blobDB.json --nohit --rank phylum --taxrule bestsum --format pdf --out phylum_$i; done
```
In specifying taxonomic rank (--rank), note that BlobTools supports superkingdom, phylum, order, family, genus, and species, but does *not* support class, kingdom, or domain.

## good luck!



