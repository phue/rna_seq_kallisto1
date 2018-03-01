How to run the pipeline

The first time you atempt to run the pipeline it might seem a bit complicated. But even though the instructions are quite long, the whole procedure is really very simple once you get the hang of it. Much easier than e.g. Galaxy!

And I am very happy to help out at any time, so don't spend too much time getting frustrating if something is not clear or if you don't manage to get something to work!

Requirements:

Computational requirements:
The only thing you'll need is a Mendel account. If you do not have one talk to the hpc office.

Data requirements:
You need to have bam files of your sequencing runs ( unaligned, demultiplexed ). If you sequenced at the vbcf your data should be availaible in the group folder under /lab/Raw/demultiplexed. If you can't find your data or have any other problem, then just contact me (Elin)

Setting up the pipeline:

Computational setup -- this step is needed ONLY the very first time you run the pipeline on Mendel:
On the login node execute the script called setup_r_packages.sh. This is done by simply writing: ./setup_r_packages.sh in the command line. The script will take a few minutes to run so please be patient.

Data setup -- this you will need to do every time you want to run a NEW analysis

The text files you need are the following:

info.tab: To get the right structure of this file you can simply open (e.g. using vi ) the file included file called info.tab and then edit it
so instead of the dummy names you have your data:
Line one: keep as is
Line two and onwards:
First column is the name of your bam file without the .bam in the end e.g. 45566_GCCAAT_CA80KANXX_6_20161111B_20161111 instead of 45566_GCCAAT_CA80KANXX_6_20161111B_20161111.bam
The second column should be the condition of the sample (e.g, wildtype, atrx_ko, what ever as long as the all replicates are given the same condition)
The third condition should be an indicator of the replicate. Eg 1,2,3... or A,B,C.

contrasts.tab: This is a very simple file, again you can use the included file as a template.
First line should list all conditions used in the info.tab
The second row and onwords represent one contrast (pair-wise comparison). So if you just want to compare one condition with another you just need one line. Here you should put a -1 for one of the conditons and a 1 for the other.

Running the pipeline:
Nextflow paramters:
If you open the file called rna_seq1.nf you will find on the very top a section called Parameters. This is a set of input information that is given to the pipeline. The parameters are run specific, meaning that one give different parameters for datasets. However, some of the parameters do not need to be changed, so the ones you need to provide are the following:

params.bam: this is the path to the folder where you have your raw bam files (note the folder should contain ONLY the bam files your are interested in)
params.seqtype: Single read (SR) or  Paired read  (PR). As most often the type is single read this is set to default, meaning that if you have single read data you do not need to change this parameter
params.strand: What type of strand specificity your data has. This is decided by the kit used for library preparation. The most common type is "reverse for" (RF-strand), this means that the read (or the first read if paired) comes from the reverse strand. This is the default setting and if you have this type of data you do not need to care about setting this parameter. The other options are: fr-stand (read, or first read in paired data, comes from forward starnd) and NULL for unstranded data (e.g from SMART2)

annotation

To acually run the pipeline you need to do the following:
in the command line write
ml Nextflow
nextflow run rna_seq1.nf

the pipeline will submit jobs to the cluster so this you do not need to take care of. However you should not close the mendel terminal while the jobs are still running. If you think that your jobs will take a long time you can use the screen command

The pipeline will generate a output folder (by default called 'results') This folder would contain everything you need and is the only thing you need to save somewhere safe. Inside this folder is another folder called report and there is a file called report.html. That is the file you should. 

