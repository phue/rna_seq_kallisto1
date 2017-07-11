#!/usr/bin/env nextflow

params.in = "$WORK/testfiles/*.bam"
params.fragment_len  = '180'
params.fragment_sd   = '20'
params.bootstrap     = '10'
params.output        = "results/"
params.fasta 	     = "/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

fasta=file(params.fasta)

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
module 'SAMtools/1.3.1-foss-2016a'

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
module 'kallisto/0.42.4-linux-x86_64'

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
module 'kallisto/0.42.4-linux-x86_64'
publishDir "$params.output/$name"
tag "fq: $name"

    input:
    file index from transcriptome_index
    set name, file(fq) from fastqs

    output:
    file "kallisto_${name}" 

    script:
    """
 	kallisto quant -i ${index} -o kallisto_${name} --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}

""" 


}



workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

  	
	
