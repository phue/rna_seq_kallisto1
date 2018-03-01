# How to run the pipeline
## Introduction

The first time you atempt to run the pipeline it might seem a bit complicated. But even though the instructions are quite long, the whole procedure is really very simple once you get the hang of it. Much easier than e.g. Galaxy!

And I am very happy to help out at any time, so don't spend too much time getting frustrating if something is not clear or if you don't manage to get something to work!

## Requirements:

### Computational requirements:
The **one and only** thing you'll need is a Mendel account. If you do not have one talk to the hpc office.

### Data requirements:
You need to have bam files of your sequencing runs (unaligned, demultiplexed). If you sequenced at the vbcf your data should be availaible in the group folder under /lab/Raw/demultiplexed/. If you can't find your data or have any other problem, then just contact me (Elin).

## Setting up the pipeline:

### Computational setup
**This step is needed ONLY the very first time you run the pipeline on Mendel**

On the login node execute the script called setup_r_packages.sh. This is done by simply writing in the command line:

./setup_r_packages.sh

The script will take a few minutes to run so please be patient.

### Data setup

**You will need to do every time you want to run a *new* analysis**

The text files you need are the following:

* info.tab : information about the samples
* contrasts.tab : defining the contrasts you want to test

**info.tab**

This file contains three columns (run_accession,condition,sample), where the first is the name of the bam file (without the .bam extension), the second one defines the sample condition (e.g. wildtype, knockout etc) and the third indictes the seperate replicates (e.g 1,2,3... or A,B,C)

**Tip1:** To get the right structure of this file you can simply open the included file info.tab file and edit it. Keep the first line (containing the headers) as is is and on line 2 and onwards insert your bam names, conditions and sample info.

**Tip2:** The condition can (should?)  be set to be more specific then the rather generic "wt" and "ko" used in the templete. E.g. one could consider using the naming guidelines for sequenicing submission.

**Example:** Let's say you have 4 samples; two Col WT replicates and two Col clf-29 replicates. Your bam files are called 12345_barcode1_extra_info.bam, 12346_barcode2_extra_info.bam, 12347_barcode3_extra_info.bam, 12348_barcode4_extra_info.bam  You want to name the replicates 1 and 2. The your info.tab would look like this:

run_accession condition sample<br/>
12345_barcode1_extra_info WT 1<br/>
12346_barcode2_extra_info WT 2<br/>
12347_barcode3_extra_info clf-29 1<br/>
12348_barcode4_extra_info clf-29 2<br/>


**contrasts.tab**

This is a very simple file, that defines which conditions you want to compare with each other. The file should have one column for each condition of interest. The name of the conditons should be the header (first line of the file). Each additional row is then one contrast (pair-wise comparison). If you want to compare condition A with condition B, so that a positive fold change means that a gene is higher expressed in A then in B, then  you define this contrast by setting a 1 in the column of the condition A and a -1 in the column of condition B.

**Tip1:** Again it is possible to use the existing contrasts.tab as a templete. Now you have to change the first row so that the names there are the same as the condition column of your info.tab. On the next line you define the contrast of your choice.

**Example:** Given the info.tab example above, let's say that you want to compare the clf-29 KO with the WT. That is genes that are upregulated in the ko will have a positive fold change, whereas downregulated genes have negative fold changes. Then the contrasts.tab should look like this:

elf-29 WT<br/>
1 -1<br/>


**Extra:**
Described above is the example where one has two condition and is interested in a comparison of the two. It is however possible to used the pipeline in more complex situations too. If there is for example three conditions (A,B,C) and one wants to for example compare A with B and B with C than this can be done too.


## Running the pipeline:

### Nextflow paramters:

If you open the file called rna_seq1.nf you will find on the very top a section called Parameters. This is a set of input information that is given to the pipeline. The parameters are run specific, meaning that one give different parameters for different datasets. However, some of the parameters do not need to be changed, so the ones you (may) need to provide are the following:

**params.bam:** this is the path to the folder where you have your raw bam files (note the folder should contain ONLY the bam files your are interested in) followed by  '/\*.bam' which tells the pipeline to take all files with the .bam extension as input files. If you follow the 'recommended project set up', then the default path will work for you and you can leave it as it is.<br/>
**params.seqtype:** this should be SR (Single read) or PR  (Paired read). As most RNA-seq is single read, SR is set to be the default, meaning that if you have single read data you do not need to change this parameter.<br/>
**params.strand:** this parameter tells the pipeline what type of strand specificity your data has. This is decided by the kit used for library preparation. The most common type is "reverse first" (RF-strand), this means that the read (or the first read if paired) comes from the reverse strand. This is the default setting and if you have this type of data you do not need to care about setting this parameter. The other options are: fr-stand (read, or first read in paired data, comes from forward starnd) and NULL for unstranded data (e.g from SMART2)<br/>
**params.anno_set:** which annotations you want to use. Options so far (more can be added on demand) are "tair10" and "araport11" where "tair10" is the default.

In addition, if you have single read data, you might consider changing the following two parameters:

**params.fragment_len**: this is an estimate of the mean fragment length. The default setting is 180 bp.<br/>
**params.fragment_sd**: this is an estimate of the fragment length standard deviasion. The default is set to 20.<br/>

In my experience those two parameters are not super critical to get 100% right. But if you know the mean and standard deviation it makes sense to provide it. 

There are two ways to provide the parameters: either you edit the rna_seq1.nf file

To acually run the pipeline you need to do the following:
in the command line write
ml Nextflow
nextflow run rna_seq1.nf

the pipeline will submit jobs to the cluster so this you do not need to take care of. However you should not close the mendel terminal while the jobs are still running. If you think that your jobs will take a long time you can use the screen command

The pipeline will generate a output folder (by default called 'results') This folder would contain everything you need and is the only thing you need to save somewhere safe. Inside this folder is another folder called report and there is a file called report.html. That is the file you should. 

