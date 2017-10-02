#!/usr/bin/env nextflow

params.in = "../merge_raw_bams/results/*.bam"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.output        = "results_araport11/"
params.organism      = "Arabidopsis"
params.anno	     = "Araport11_genes_and_non_coding" // "TAIR10" "Araport11_genes"
params.design        = "exp.txt"
params.contrast      = "contrast.txt"


if(params.organism == "Arabidopsis") {
	effSize=119146348
	binSize=10
}

if(params.anno == "Araport11_genes" || params.anno == "Araport11_genes_and_non_coding"){
	gtf=file("/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Araport11_GFF3_genes_transposons.201606.gtf")
	starDir="star_araport" // not affect kallisto
	fasta_dna=file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.Araport.dna.toplevel.fa")
}

if(params.anno == "Araport11_genes_and_non_coding"){
	fasta=file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Araport11_genes_and_non_coding.201606.cdna.fasta.gz")
	kallisto_index="araport11_ganc_transcripts.idx"
}

if(params.anno == "Araport11_genes"){
        fasta=file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Araport11_genes.201606.cdna.fasta.gz")
        kallisto_index="araport11_transcripts.idx"
}


design=file(params.design)
contrast=file(params.contrast)

//START

bams = Channel
	.fromPath(params.in)
	.map { file -> tuple(file.baseName, file) } 


process bam2fastq {
publishDir "$params.output/$name", mode: 'copy'
tag "bam: $name"

    input:
    set name, file(reads) from bams

    output:
    set name, file("${name}.fastq") into fastqs

    script: 
    """ 
    samtools bam2fq $reads > ${name}.fastq
    """
     
}
fastqs.into { fastqs_kallisto; fastqs_star }

process kallistoIndex {
storeDir '/lustre/scratch/projects/berger_common/backup_berger_common'
    input:
    file fasta

    output:
    file "${kallisto_index}" into transcriptome_index

    script:
    """
    kallisto index -i ${kallisto_index} ${fasta} 
    """
}

process quantKallisto {
publishDir "$params.output/$name", mode: 'copy'
tag "fq: $name"

    input:
    file index from transcriptome_index
    set name, file(fq) from fastqs_kallisto

    output:
    file "kallisto_${name}" into kallisto_dirs 

    script:
    """
 	kallisto quant -i ${index} -o kallisto_${name} --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}
    """ 
}
process STARindex {
storeDir '/lustre/scratch/projects/berger_common/backup_berger_common/'

    input:
    file fasta_dna
    file gtf

    output: 
    file "${starDir}" into star_index

    script:
    """
    mkdir ${starDir}
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${starDir} --genomeFastaFiles ${fasta_dna} --sjdbGTFfile ${gtf}
    """
}


process STAR {
publishDir "$params.output/$name", mode: 'copy'
tag "star: $name"

    input:
    file index from star_index
    set name, file(fq) from fastqs_star
    
    output:
    set name, file("star_${name}/${name}Aligned.sortedByCoord.out.bam") into sort_bam    

    script:
    """
    mkdir -p star_${name}
    STAR --genomeDir $index --readFilesIn $fq --runThreadN 4  --outFileNamePrefix ./star_${name}/${name} --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
    """
}
/*process bigWig {
publishDir "$params.output/$name", mode: 'copy'
tag "bw: $name"
    
    input:
    set name. file(bam) from sort_bam
    
    output:
    file ${name}.coverage.bw   

    script:
    """
    bamCoverage -b ${bam} -o ${name}.coverage.bw --normalizeTo1x ${effSize} --binSize=${binSize} 

*/
process deseq2 {
publishDir "$params.output/deseq", mode: 'copy'

  input:
  file 'kallisto/*' from kallisto_dirs.collect()
  file design
  
  output:
  file 'pairs.pdf'
  file 'dds.Rdata'
  file 'pca.pdf'
  file 'maplots.pdf'
  file 'contrast_*'

  script:
  """
  $baseDir/bin/deseq2.R kallisto ${design} ${contrast} 
  """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

  	
	
