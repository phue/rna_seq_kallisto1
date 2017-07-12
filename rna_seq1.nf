#!/usr/bin/env nextflow

params.in = "../bams/*.bam"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '100'
params.output        = "results/"
params.fasta 	     = "/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
params.dna_fasta     = "/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
params.gtf 	     = "/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Arabidopsis_thaliana.TAIR10.35.gtf"
params.design        = "exp.txt"

fasta=file(params.fasta)
fasta_dna=file(params.dna_fasta)
gtf=file(params.gtf)
design=file(params.design)

/*Channel
    .fromFilePairs( params.in, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.in}" }
    .set { bams }*/

bams = Channel
	.fromPath(params.in)
	.map { file -> tuple(file.baseName, file) } 


process bam2fastq {
publishDir "$params.output/$name"
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
    file "tair10_transcripts.idx" into transcriptome_index

    script:
    """
    kallisto index -i tair10_transcripts.idx ${fasta} 
    """
}

process quantKallisto {
publishDir "$params.output/$name"
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
    file 'star' into star_index

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star --genomeFastaFiles ${fasta_dna} --sjdbGTFfile ${gtf} 
    """
}

process STAR {
publishDir "$params.output/$name"
tag "star: $name"

    input:
    file index from star_index
    set name, file(fq) from fastqs_star
    
    output:
    file "star_${name}"    

    script:
    """
    mkdir -p star_${name}
    STAR --genomeDir $index --readFilesIn $fq --runThreadN 4 --quantMode GeneCounts --outFileNamePrefix ./star_${name}/
    """
}
process deseq2 {
publishDir "$params.output/deseq"

  input:
  file 'kallisto/*' from kallisto_dirs.collect()
  file design
  
  output:
  file 'pairs.pdf'

script:
"""
$baseDir/bin/deseq2.R kallisto ${design} 
"""
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

  	
	
