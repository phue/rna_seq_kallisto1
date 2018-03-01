How to run the pipeline

To run this pipeline you need two things:

- A mendel account, if you do not have one talk to the hpc office
- Bam files of your sequencing runs ( unaligned, demultiplexed ). If you sequenced at the vbcf your data should be availaible in the group folder under /lab/Raw/demultiplexed

If you have the above things than you still have a few steps that you need to do.

The first time only, that you run the pipeline you need to, on the login node, execute the script called setup_r_packages.sh. This is done by simply writing: ./setup_r_packages.sh in the command line. It might take some time but please patient.

the input files you need are the following:

info.tab: To get the right structure of this file you can simply open (e.g. using vi ) the file included file called info.tab and then edit it
so instead of the dummy names you have your data:
Line one: keep as is
Line two and onwards:
First column is the name of your bam file without the .bam in the end e.g. 45566_GCCAAT_CA80KANXX_6_20161111B_20161111 instead of 45566_GCCAAT_CA80KANXX_6_20161111B_20161111.bam
The second column should be the condition of the sample (e.g, wildtype, atrx_ko, what ever as long as the all replicates are given the same condition)
The third condition should be an indicator of the replicate. Eg 1,2,3... or A,B,C.

BAM FILES

contrasts.tab: This is a very simple file, again you can use the included file as a template.
First line should list all conditions used in the info.tab
The second row and onwords represent one contrast (pair-wise comparison). So if you just want to compare one condition with another you just need one line. Here you should put a -1 for one of the conditons and a 1 for the other.

To acually run the pipeline you need to do the following:
in the command line write
ml Nextflow
nextflow run rna_seq1.nf

the pipeline will submit jobs to the cluster so this you do not need to take care of. However you should not close the mendel terminal while the jobs are still running. If you think that your jobs will take a long time you can use the screen command

The pipeline will generate a output folder (by default called 'results') This folder would contain everything you need and is the only thing you need to save somewhere safe. Inside this folder is another folder called report and there is a file called report.html. That is the file you should. 

