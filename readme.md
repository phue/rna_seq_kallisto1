# How to run the pipeline

If you at any time have **any problems**; with the manual, the pipeline or anything else related to this, **just contact me (Elin)**

**NB** This document is focused on the case where the starting point is an unaligned bam file. It is also possible to use fastq files, see [Using fastq files](#using-fastq-files).

## Introduction

The first time you attempt to run the pipeline it might seem a bit complicated. But even though the instructions are quite long, the whole procedure is really very simple once you get the hang of it. Much easier than e.g. Galaxy.

Using this pipeline you can start your job with one (okay two) command(s) and then go on and do other stuff while the analysis is running.

Another advantage is that the output of the pipeline contains all information about software, versions etc. It also provides a full documentation of every single step that has been done. This means that
1.  your results can be reproduced, by you or by others
2.  it's easy to write this methods part of your paper

**Tip:** I am happy to help out at any time, so don't spend too much time getting frustrating if something is not clear or if you don't manage to get something to work!

## Requirements

### Computational requirements
The **one and only** thing you'll need is a Mendel account. If you do not have one talk to the hpc office. Oh, and access to a computer of course...

### Data requirements
You need to have bam files of your sequencing runs (unaligned, demultiplexed). If you have sequenced at the vbcf, then your data should be available in the group folder under /lab/Raw/demultiplexed/. If you can't find your data or have any other problem, then just contact me (Elin).

## Get pipeline 

### Using git (recommended)

If you have a github account that is added to the GMI github, you can use git to clone the pipeline repository by simply type the following:

git clone https://github.com/Gregor-Mendel-Institute/rna_seq_kallisto1.git

Using git comes with a lot of benefits, espacially if you want to make any changes/additions to the code. It is also useful if I have to fix bugs or if you ask me to implement any new features. You can sign up for github here: <https://github.com>. If you email Ümit <uemit.seren@imp.ac.at> he will then add you to the GMI account. If you have any questions about github - feel free to ask me.

### Without git
If you don't have access to the GMI github and have no interest in getting this set up for you, there is an alternative way. A copy of the repository is located at /lustre/scratch/projects/berger_common/pipelines/. You can copy this folder to your current folder by typing (**NB the dot at the end of the line**):

cp -r /lustre/scratch/projects/berger_common/pipelines/rna_seq_kallisto1 .

## Recommended project setup
In the work directory on mendel ($WORK), create a folder called something that fits your project. Inside this folder create a subfolder called bams. Copy (using the data moving node) your bam files (see [Data requirements](#data-requirements)) into this folder. **NB the bams folder should only contain the bam files you want to include in the analysis!!**. Next get the pipeline code (see section [Get pipeline](#get-pipeline)). Now you will have two subfolders, the bams folder that you created and a new folder called rna_seq_kallisto1.

The following lines shows how to do what was described above. Here "my_cool_project" is the name you want to give to your project and "/full/path/to/bam/file" is the path to the location of the bam file (e.g. for a file named myBAM that is in the lab folder: /net/gmi.oeaw.ac.at/berger/lab/Raw/demultiplexed/myBAM). All but the data copying step can be done on the login or data moving node. The copying step should be done on a data moving node. </br>

cd $WORK</br>
mkdir my_cool_project </br>
cd my_cool_project </br>
mkdir bams </br>

  \# Use **data moving node** e.g. dmn0 for the next step </br>
cp /full/path/to/bam/file bams  </br>
\# repeat the above line for for each bam file in your project</br>

\# With git **as recommended**- use this step if you have access to GMI github, then skip the next step. </br>
git clone https://github.com/Gregor-Mendel-Institute/rna_seq_kallisto1.git </br>

\# Without git - use this step if you do not have access to GMI github, then skip the previous step. (**NB the dot at the end**) </br>
cp -r /lustre/scratch/projects/berger_common/pipelines/rna_seq_kallisto1 . </br>

If you now type:</br>
ls</br>
you will see that you have the two subfolders mentioned above (bams and rna_seq_kallisto1).

Now you are ready for the next step. First move into rna_seq_kallisto1:</br>
cd rna_seq_kallisto1</br>
here you will now have some folders and files. Type:</br>
ls</br>
to list them. The ones you need to care about are: *info.tab*, *contrasts.tab* and to some extent *rna_seq1.nf*. What you have to do with those three files is described in this documentation in the sections [Data setup](#data-setup) and [Nextflow parameters of interest](#nextflow-parameters-of-interest).



## Data setup

**You will need to do this every time you want to run a *new* analysis**

The text files you need to edit are the following:

* info.tab : information about the samples
* contrasts.tab : defining the contrasts you want to test

### info.tab

This file contains three columns separated by a ",". The first is named "run_accession" and is the name of the bam file (**without the .bam extension**), the second is called "condition" and defines the sample condition (e.g. wildtype, knockout etc) and the last is called "sample" and defines the separate replicates (e.g 1,2,3... or A,B,C)

**Tip1:** To get the right structure of this file you can simply open the included file info.tab file and edit it. Keep the first line (containing the headers) as it is and on line 2 and onwards insert your bam names, conditions and sample info.

**Tip2:** The condition can (should?)  be set to be more specific then the rather generic "wt" and "ko" used in the templete. E.g. one could consider using the naming guidelines for sequencing submission.

**Example:** Let's say you have 4 samples; two Col WT replicates and two Col clf-29 replicates. Your bam files are called 12345_barcode1_extra_info.bam, 12346_barcode2_extra_info.bam, 12347_barcode3_extra_info.bam, 12348_barcode4_extra_info.bam. You want to name the replicates 1 and 2. Then your info.tab would look like this:

run_accession,condition,sample<br/>
12345_barcode1_extra_info,WT,1<br/>
12346_barcode2_extra_info,WT,2<br/>
12347_barcode3_extra_info,clf-29,1<br/>
12348_barcode4_extra_info,clf-29,2<br/>


### contrasts.tab

This is a very simple file, it only defines which conditions you want to compare with each other. The file should have one column for each condition. The name of the conditions should be the header (first line) of the file. Each additional row is then one contrast (pair-wise comparison). If you want to compare condition A with condition B, so that a positive fold change means that a gene is higher expressed in A then in B, then you define this contrast by setting a 1 in the column of the condition A and a -1 in the column of condition B.

**Tip1:** Again it is possible to use the existing contrasts.tab as a template. Now you have to change the first row so that the names there are the same as the condition column of your info.tab. On the next line you define the contrast of your choice.

**Example:** Given the info.tab example above, let's say that you want to compare the clf-29 KO with the WT. That is genes that are up-regulated in the ko will have a positive fold change, whereas down-regulated genes have negative fold changes. Then the contrasts.tab should look like this:

elf-29,WT<br/>
1,-1<br/>


**Extra:**
Described above is the example where one has two condition and is interested in a comparison of the two. It is however possible to used the pipeline in more complex situations too. If there is for example three conditions (A,B,C) and one wants to for example compare A with B and B with C than this can be done too.:

A,B,C</br>
1,-1,0</br>
0,1,-1</br>

This will result in two different comparisons: First one where A is compared to B (line 2) and one where B is compared to C (line3).


## Running the pipeline

### Nextflow parameters of interest

If you open the file called rna_seq1.nf in the rna_seq_kallisto1 folder you will find, on the very top, a section called "Parameters". This section contains is a set of input information that is given to the pipeline. The parameters are run-specific, meaning that one can give different parameters for different datasets. However, some of the parameters here do not need to be changed, so the ones you (may) need to provide are the following:

**params.files:** this is the path to the folder where you have your raw bam files followed by  '/\*.bam' which tells the pipeline to take all files with the .bam extension as input files (NB the folder should contain ONLY the bam files your are interested in). **If you follow the 'recommended project setup', then the default path will work for you and you can leave it as it is.**<br/>

**params.seqtype:** this should be set to SR (Single read) or PR  (Paired read). As most RNA-seq is single read, SR is set to be the default, meaning that **if you have single read data you do not need to change this parameter.**<br/>

**params.strand:** this parameter tells the pipeline what type of strand specificity your data has. This is decided by the kit used for library preparation. The most common type is "reverse first" (RF-strand), this means that the read (or the first read if paired) comes from the reverse strand. RF-strand is the default setting so **if you have "reverse first" strand specificity you do not need to change this parameter.** The other options are: fr-stand (read, or first read in paired data, comes from forward strand) and NULL for un-stranded data (e.g from SMART2)<br/>

**params.anno_set:** which annotations you want to use. Options so far (more can be added on demand) are "tair10" and "araport11" and "tair10_TE" where "tair10" is the default and the "tair10_TE" was added on request by Maggie and Bhagyshree (contact them for details on what is included).

In addition, if you have single read data, you might consider changing the following two parameters:

**params.fragment_len**: this is an estimate of the mean fragment length. The default setting is 180 bp.<br/>
**params.fragment_sd**: this is an estimate of the fragment length standard deviation. The default is set to 20.<br/>

In my experience it's not super critical to get those two parameters 100% right. But if you know the mean and standard deviation it makes sense to provide this information.

### How to use non-default parameters

If you want to use non-default parameters there are two ways to pass your parameters to the pipeline:

**OPTION 1: Editing the rna_seq1.nf file**

You can edit the rna_seq1.nf file so that the file defined parameters suit you. In the rna_seq1.nf file anything after a // is a comment. In the parameters section the possible parameter options are listed as a comment.

**Example:**
If you have paired data you would change the line:<br/>
params.seqtype = 'SR' // 'PR' <br/>
to <br/> 
params.seqtype = 'PR'

**OPTION 2: Include parameters in pipeline run**

The other way is to provide the parameters when calling the pipeline by adding --paraName Value to the call, e.g. for above example add --seqtype 'PR' to the run command. See also [Starting the pipeline](#starting-the-pipeline) for more information and an example.

### Starting the pipeline

To actually run the pipeline you first need to load the Nextflow module by typing in the command line of the login node:<br/>
ml Nextflow

Then make sure you are in the subfolder rna_seq_kallisto1.
Start the nextflow pipeline by typing, again on the command line of the login node:<br/>
nextflow run rna_seq1.nf

Please remember that if you did not edit the rna_seq1.nf file with your parameters you will need to add any non-default parameters to the run command, e.g.:<br/>
nextflow run rna_seq1.nf --seqtype 'PR'

The pipeline will now do the rest, it will for example load modules and submit jobs to the cluster without you having to do anything. However *you should not close the mendel terminal while the jobs are still running*. If you think that your jobs will take a long time you can use the screen command. If you don't know it google 'screen command' or ask someone who might know, for example me (Elin).

### Advanced run options

This section is for advanced users. In most cases you will not need to read this (unless you are interested!)

**Run in background**</br>
If you want to continue working in the terminal where you launch the pipeline, then you can use the -bg option (bg = background). Nextflow will still send messages to the terminal when new processes are being submitted.

nextflow run -bg  rna_seq1.nf<br/>


**Restart failed/interrupted run**</br>
If your run fails for reasons that you know*, e.g. you killed it, you can use the resume option. This means that any process that had finished before the pipeline failed will be "reused" and the pipeline will not spend time on re-running those.

nextflow run rna_seq1.nf -resume

*If the pipeline fails and it's not obvious why - contact me (Elin)


**More options**</br>
The pipeline uses Nextflow, a software developed for reproducible scientific workflows. If you type the following commands in login node terminal:</br>

ml Nextflow</br>
nextflow -h</br>
and/or</br>
nextflow run -h</br>

Then all available options and commands will be listed in your terminal.

For more information about Nextflow see https://www.nextflow.io/ and read the paper: P. Di Tommaso, et al. *Nextflow enables reproducible computational workflows.* Nature Biotechnology 35, 316–319 (2017) [doi:10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820)


### Output

The pipeline will generate a output folder (by default called 'results') This folder would contain everything you need (with the potential exeption of the aligned bam files - see below) and is the only thing you need to save somewhere safe.
To do this you need to use the **data moving node** again. E.g.</br>

cp -r $WORK/my_cool_project/rna_seq_kallisto1/results to/my/storage </br>

where "to/my/storage" is the full path location where you want to keep the results. </br>
Inside this folder is another folder called "report" and there is a file called "report.html". That file you should open, it will explain the outputs and direct you to the different files.

**Aligned BAM files**</br>
As most users do not need the aligned bam files and since they are usually quite big, those are not included in the 'results' output folder. Instead they are copied into a folder called 'result_bams' where you easily can access them should you be interested.

## Using fastq files

The pipeline can also take fastq files as input. The easiest way to get this to work is described below:

Follow the instructions in [Recommended project setup](#recommended-project-setup). But instead of creating a "bams" folder create a "fastqs" folder. In this folder you put all fastq files you want include in the analysis. The fastq files should be have names ending with _1.fastq if you have single read data. If you have paired end data you should have one _1.fastq and one _2.fastq file for each sample.

In the rna_seq1.nf file you then have to change the parameter **params.type** from "bam" to "fastq". All other steps are the same as when using bam files.

