# How to run the pipeline
## Introduction

The first time you atempt to run the pipeline it might seem a bit complicated. But even though the instructions are quite long, the whole procedure is really very simple once you get the hang of it. Much easier than e.g. Galaxy!


And I am very happy to help out at any time, so don't spend too much time getting frustrating if something is not clear or if you don't manage to get something to work!

## Requirements

### Computational requirements
The **one and only** thing you'll need is a Mendel account. If you do not have one talk to the hpc office.

### Data requirements
You need to have bam files of your sequencing runs (unaligned, demultiplexed). If you sequenced at the vbcf your data should be availaible in the group folder under /lab/Raw/demultiplexed/. If you can't find your data or have any other problem, then just contact me (Elin).

## Get pipeline

### Using git
If you have a github account that is added to the GMI github, you can use git to clone the pipeline resposity by simply type the following:

git clone

### Without git
If you don't have access to the GMI github and have no inteters in getting this set up for you, there is an alternative way. A copy of the respository is located at ..... . You can cp this folder:

cp -r .../.... .

## Recommended project set up
In the work directory on mendel ($WORK), create a folder called something that fits your project. Inside this folder create a subfolder called bams. Copy (using the data moving node) your bam files (see [Data requirements](#Data-requirements)) into this folder. Next get the pipeline code (see section [get pipeline](#Get-pipeline). Now you will have two subfolders, the bams from before and a new folder called rna_seq_kallisto1. Move into rna_seq_kallisto1, here you will now have some folders and files. The ones you need to care about are: info.tab, contrasts.tab and to some extent rna_seq1.nf. What you have to do with those three files is decribed in this documentention in the sections [Data setup](#Data-setup) and [Nextflow parameters](#Nextflow-parameters).

The following lines of code is an example where "my_cool_project" is the name you want to give to your project and /full/path/to/bam/file is the path to the location of the bam file (e.g. for a file named myBAM that is in the lab folder: /net/gmi.oeaw.ac.at/berger/lab/Raw/demultiplexed/myBAM)

/# all but the data copying step can be done on the login or data moving node. The copying step should be done on a data moving node.
cd $WORK /# changing directory to the work directory
mkdir my_cool_project /# creating a new folder for project
cd my_cool_project /# moving into the new folder
mkdir bams /# creating a new folder for the bam files
 \# this step should be done on a data moving node e.g. dmn0
 cp /full/path/to/bam/file bams \#  copying the bam files to the new folder - you will need to do this for each bam file 
\# with git

\# without git
cp -r ////rna_seq_kallisto1 .
cd rna_seq_kallisto1



## Data Setup

**You will need to do every time you want to run a *new* analysis**

The text files you need are the following:

* info.tab : information about the samples
* contrasts.tab : defining the contrasts you want to test

### info.tab

This file contains three columns (run_accession,condition,sample), where the first is the name of the bam file (without the .bam extension), the second one defines the sample condition (e.g. wildtype, knockout etc) and the third indictes the seperate replicates (e.g 1,2,3... or A,B,C)

**Tip1:** To get the right structure of this file you can simply open the included file info.tab file and edit it. Keep the first line (containing the headers) as it is and on line 2 and onwards insert your bam names, conditions and sample info.

**Tip2:** The condition can (should?)  be set to be more specific then the rather generic "wt" and "ko" used in the templete. E.g. one could consider using the naming guidelines for sequenicing submission.

**Example:** Let's say you have 4 samples; two Col WT replicates and two Col clf-29 replicates. Your bam files are called 12345_barcode1_extra_info.bam, 12346_barcode2_extra_info.bam, 12347_barcode3_extra_info.bam, 12348_barcode4_extra_info.bam  You want to name the replicates 1 and 2. The your info.tab would look like this:

run_accession condition sample<br/>
12345_barcode1_extra_info WT 1<br/>
12346_barcode2_extra_info WT 2<br/>
12347_barcode3_extra_info clf-29 1<br/>
12348_barcode4_extra_info clf-29 2<br/>


### contrasts.tab

This is a very simple file, that defines which conditions you want to compare with each other. The file should have one column for each condition of interest. The name of the conditons should be the header (first line of the file). Each additional row is then one contrast (pair-wise comparison). If you want to compare condition A with condition B, so that a positive fold change means that a gene is higher expressed in A then in B, then you define this contrast by setting a 1 in the column of the condition A and a -1 in the column of condition B.

**Tip1:** Again it is possible to use the existing contrasts.tab as a templete. Now you have to change the first row so that the names there are the same as the condition column of your info.tab. On the next line you define the contrast of your choice.

**Example:** Given the info.tab example above, let's say that you want to compare the clf-29 KO with the WT. That is genes that are upregulated in the ko will have a positive fold change, whereas downregulated genes have negative fold changes. Then the contrasts.tab should look like this:

elf-29 WT<br/>
1 -1<br/>


**Extra:**
Described above is the example where one has two condition and is interested in a comparison of the two. It is however possible to used the pipeline in more complex situations too. If there is for example three conditions (A,B,C) and one wants to for example compare A with B and B with C than this can be done too.:

A B C</br>
1 -1 0</br>
0 1 -1</br>

This will result in two different comparisons: First one where A is compared to B (line 2) and one where B is compared to C (line3).


## Running the pipeline

### Nextflow parameters

If you open the file called rna_seq1.nf you will find on the very top a section called Parameters. This is a set of input information that is given to the pipeline. The parameters are run specific, meaning that one give different parameters for different datasets. However, some of the parameters do not need to be changed, so the ones you (may) need to provide are the following:

**params.bam:** this is the path to the folder where you have your raw bam files (note the folder should contain ONLY the bam files your are interested in) followed by  '/\*.bam' which tells the pipeline to take all files with the .bam extension as input files. If you follow the 'recommended project set up', then the default path will work for you and you can leave it as it is.<br/>
**params.seqtype:** this should be SR (Single read) or PR  (Paired read). As most RNA-seq is single read, SR is set to be the default, meaning that if you have single read data you do not need to change this parameter.<br/>
**params.strand:** this parameter tells the pipeline what type of strand specificity your data has. This is decided by the kit used for library preparation. The most common type is "reverse first" (RF-strand), this means that the read (or the first read if paired) comes from the reverse strand. This is the default setting and if you have this type of data you do not need to care about setting this parameter. The other options are: fr-stand (read, or first read in paired data, comes from forward starnd) and NULL for unstranded data (e.g from SMART2)<br/>
**params.anno_set:** which annotations you want to use. Options so far (more can be added on demand) are "tair10" and "araport11" where "tair10" is the default.

In addition, if you have single read data, you might consider changing the following two parameters:

**params.fragment_len**: this is an estimate of the mean fragment length. The default setting is 180 bp.<br/>
**params.fragment_sd**: this is an estimate of the fragment length standard deviasion. The default is set to 20.<br/>

In my experience those two parameters are not super critical to get 100% right. But if you know the mean and standard deviation it makes sense to provide it. 

There are two ways to provide the parameters: either you edit the rna_seq1.nf file so that the file defined parameters suit you. By the way anything behind // is a comment, here in the parameters section the possible options are listed in this comment. So if you have paired data you would change the line:<br/>

params.seqtype = 'SR' // 'PR' <br/>
to <br/> 
params.seqtype = 'PR'

The other way is to provide the parameters when calling the pipeline by adding --paraName Value to the call, e.g for above example --seqtype 'PR'.

### Starting the pipeline

To acually run the pipeline you first need to do load the Nextflow module by typing in the command line of the login node:<br/>
ml Nextflow


Now you are ready to call the nextflow pipeline by typing, again on the command line of the login node:<br/>
nextflow run rna_seq1.nf

please remember that if you did not edit the .nf file with your parameters you will need to add any non-default parameters, e.g.:<br/>
nextflow run rna_seq1.nf --seqtype 'PR'

The pipeline will submit jobs to the cluster so this you do not need to take care of. However you should not close the mendel terminal while the jobs are still running. If you think that your jobs will take a long time you can use the screen command. If you don't know it google 'screen command' or ask someone who might now, for example me (Elin).

### Output

The pipeline will generate a output folder (by default called 'results') This folder would contain everything you need and is the only thing you need to save somewhere safe. Inside this folder is another folder called report and there is a file called report.html. That file you should open, it will explain the outputs and direct you to the different files. 




