#!/usr/bin/env nextflow

/*####################################
  parameters
#####################################*/

params.bam 		= "../bams/*.bam"
params.fragment_len  	= '180'
params.fragment_sd  	= '20'
params.bootstrap     	= '100'
params.seqtype 		= 'SR' // 'PR'
params.strand 		= 'rf-stranded'//  fr-stranded,  NULL
params.output        	= "results/"
params.info 		= 'info.tab' // name, type, condition  
params.anno_set 	= "araport_genes" // "atair10" // "araport_genes" "araport_genes_non"
params.deseq_type 	= "kallisto" // "star" // "NULL
params.contrast         = "contrasts.tab"  

/*##################################
  annotation set selection
####################################*/

if(params.anno_set == "tair10"){
	fasta_dna = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
	gtf = file("/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Arabidopsis_thaliana.TAIR10.35.gtf") 
	fasta = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.cdna.all.fa")
	starDir = "star_tair10"
	kallistoDir = "tair10_transcripts.idx"
	params.normtosize = '119146348'
	txdb="tair10"
}
if(params.anno_set == "araport_genes"){
	fasta_dna = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
	gtf = file("/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Araport11_GFF3_genes_transposons.201606.gtf") 
	fasta = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Araport11_genes.201606.cdna.fasta.gz")
	starDir = "star_araport"
	kallistoDir = "araport_genes.idx"
	params.normtosize = '119146348'
	txdb=file("/lustre/scratch/projects/berger_common/backup_berger_common/araport11.txdb")
}
if(params.anno_set == "araport_genes_non"){
	fasta_dna = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
	gtf = file("/lustre/scratch/projects/berger_common/backup_berger_common/gtf/Araport11_GFF3_genes_transposons.201606.gtf")
	fasta = file("/lustre/scratch/projects/berger_common/backup_berger_common/fasta/Araport11_genes_and_non_coding.201606.cdna.fasta.gz")
	starDir = "star_araport"
	kallistoDir = "araport_gene_non.idx"
	params.normtosize = '119146348'
	txdb=file("/lustre/scratch/projects/berger_common/backup_berger_common/araport11.txdb")
}


/*###################################
  Analysis start
####################################*/

design = file(params.info)

bam_files = Channel
          .fromPath(params.bam)
          .map { file -> [ id:file.baseName,file:file] }

// SORT BAM

process sortBam {
tag "sort: $id"

        input:
        set id, file(bam) from bam_files

        output:
        set id, file("${id}.sort.bam") into bam_sorted
	
        script:
        """
        samtools sort -n $bam -o ${id}.sort.bam
        """
}

// BAM TO FASTQ

process generateFastq {
tag "bam : $name, type:$params.seqtype"
//publishDir "${params.output}/fastq", mode: 'copy'

        input:
        set  name, file(bam) from bam_sorted

        output:
        set name, file('*.fastq') into fastqs

        script:
        if (params.seqtype=='SR'){
        """
        bedtools bamtofastq -i ${bam} -fq ${name}_1.fastq
        """
        }
        else {
        """
        bedtools bamtofastq -i ${bam} -fq ${name}_1.fastq -fq2 ${name}_2.fastq
        """
        }
}

// COPY CHANNEL

fastqs.into { fastqs_kallisto; fastqs_star }

// KALLISTO INDEX IF NEEDED

process kallistoIndex {
tag "dir: $kallistoDir"
storeDir '/lustre/scratch/projects/berger_common/backup_berger_common'

   	input:
    	file fasta

    	output:
    	file "${kallistoDir}" into transcriptome_index

    	script:
    	"""
    	kallisto index -i ${kallistoDir} ${fasta} 
    	"""
}

// KALLIST QUANT

process quantKallisto {
tag "fq: $name "

	input:
    	file index from transcriptome_index
    	set name, file(fq) from fastqs_kallisto

    	output:
    	file "kallisto_${name}" into kallisto_dirs

    	script:
    	def single = fq instanceof Path
    	if( single  && params.strand ==null) {
    	"""
	mkdir kallisto_${name}
    	kallisto quant -i ${index} -o kallisto_${name}  --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}  
    	"""
    	}
    	else if( single ){
        """
        mkdir kallisto_${name}
        kallisto quant -i ${index} -o kallisto_${name} --${params.strand} --single -l ${params.fragment_len} -s ${params.fragment_sd} -b ${params.bootstrap} ${fq}
        """
	}
	else if (params.strand ==null) {
    	"""
	mkdir kallisto_${name}
    	kallisto quant -i ${index} -o kallisto_${name} -b ${params.bootstrap} ${fq}
    	"""
    	}
	else {
	"""
        mkdir kallisto_${name}
        kallisto quant -i ${index} -o kallisto_${name} -b --${params.strand} ${params.bootstrap} ${fq}
	"""
	}
}

// COMBINE KALLISTO OUTPUT

process kallistoCountMatrix {
	tag "anno: ${params.anno_set}"
	publishDir "$params.output/kallisto_data" , mode: 'copy'	

	input:
	file 'kallisto/*' from kallisto_dirs.collect()
	

	output:
	file 'kallisto_counts.tab' 
	file 'kallistoData.Rdata' into kallistodata 	

	script:
	"""
	$baseDir/bin/sumkallisto.R kallisto ${txdb} ${design}
	"""
}

// STAR INDEX IF NEEDED

process STARindex {
	tag "dir: $starDir"
	storeDir '/lustre/scratch/projects/berger_common/backup_berger_common/'

   	input:
    	file fasta_dna
    	file gtf

    	output: 
    	file "${starDir}" into star_index

    	script:
    	"""
    	mkdir -p  ${starDir}
    	STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${starDir} --genomeFastaFiles ${fasta_dna} --sjdbGTFfile ${gtf} 
   	 """
}

// STAR ALIGN

process STAR {
	tag "star: $name"

   	input:
    	file index from star_index
    	set name, file(fq) from fastqs_star
    
    	output:
	set name, file("star_${name}/${name}Aligned.sortedByCoord.out.bam") into sort_bam    
	set name, file("star_${name}/${name}Log.final.out") into final_log    
    	file "star_${name}/${name}ReadsPerGene.out.tab" into starcount

    	script:
    	"""
	mkdir -p star_${name}
    	STAR --genomeDir $index --outFileNamePrefix ./star_${name}/${name} --readFilesIn  $fq --runThreadN 4 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate 
    	"""
}

// COMBINE STAR COUNTS

process starCountMatrix {
	tag "strand: ${params.strand}"
  	publishDir "$params.output/star_data" , mode: 'copy'

	input:
	file 'star/*' from starcount.collect()
	
	output:
	file 'star_counts.tab'
	file 'starData.Rdata' into stardata
	
	script:
	"""
	$baseDir/bin/sumstar.R star ${params.strand} ${design} 
	"""
}

// BAM 2 BW

process bam2bw {
	publishDir "$params.output/$name/bam_bw", mode: 'copy'
        tag "bw: $name"

	input:
	set name, file(bam) from sort_bam


	output:
	file("${name}.bw")
	
	script:
	"""
	export TMPDIR=\$(pwd)
   	samtools index ${bam}
   	bamCoverage -b ${bam} -o ${name}.bw --normalizeTo1x ${params.normtosize} --binSize=10	
	"""
}



/*
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
*/
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

  	
	
