#!/usr/bin/env nextflow

params.in = "../../testfiles/*.bam"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '10'
params.output        = "results/"

/*fasta = Channel.fromPath( '../testfiles/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz' )*/
fasta=file('../../testfiles/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz')
design=file('exp.txt')

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


process kallistoIndex {

    input:
    file fasta

    output:
    file "transcriptome.index" into transcriptome_index

    script:
    """
    kallisto index -i transcriptome.index ${fasta} 
    """
}

process quantKallisto {
publishDir "$params.output/$name"
tag "fq: $name"

    input:
    file index from transcriptome_index
    set name, file(fq) from fastqs

    output:
    file "kallisto_${name}" into kallisto_dirs 

    script:
    """
 	kallisto quant -i ${index} -o kallisto_${name} --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}

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

  	
	
