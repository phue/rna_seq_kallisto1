#!/usr/bin/env nextflow

/**************
 * Parameters
 **************/

params.bam 		= "../bams/*.bam"
params.fragment_len  	= '180'
params.fragment_sd  	= '20'
params.bootstrap     	= '100'
params.seqtype 		= 'SR' // 'PR'
params.strand 		= 'rf-stranded'//  fr-stranded,  NULL
params.output        	= "results/"
params.info 		= 'info.tab' // name, type, condition  
params.anno_set 	= "araport_genes" // "atair10"  
//params.deseq_type 	= "kallisto" // "star" // "NULL"
params.contrast         = "contrasts.tab"  
params.pvalue		= 0.1
params.binsize		= 10
//params.storage		= "/lustre/scratch/projects/berger_common/backup_berger_common/"

//fasta_dna, fasta, gtf, params.normtosize, txdb, 

/***************
 *  annotation set selection
 ***************/

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

report = file("report/deseq2.Rmd")

log.info "RNA-SEQ N F  ~  version 0.1"
log.info "====================================="
log.info "bam files         	: ${params.bam}"
log.info "fragment length 	: ${params.fragment_len}"
log.info "fragment sd		: ${params.fragment_sd}"
log.info "bootstrap		: ${params.bootstrap}"
log.info "seq type 		: ${params.seqtype}"
log.info "strandness		: ${params.strand}"
log.info "output		: ${params.output}"
log.info "sample info 		: ${params.info}"
log.info "annotations		: ${params.anno_set}"
//log.info "DESeq2 data		: ${params.deseq_type}"
log.info "contrasts		: ${params.contrast}"
log.info "p-value 		: ${params.pvalue}"
log.info "norm. size		: ${params.normtosize}"
log.info "binsize		: ${params.binsize}"
log.info "txdb 			: ${txdb}"
log.info "fasta dna		: ${fasta_dna}"
log.info "fasta			: ${fasta}"
log.info "\n"


/*********************************************
**********************************************
ANALYSIS START
**********************************************
*********************************************/

/*
 * Input parameters validation
 */

design = file(params.info)
contrasts = file(params.contrast) 

/*
 * validate input files
 */

if( !design.exists() ) exit 1, "Missing sample info file: ${design}"
if( !contrasts.exists() ) exit 1, "Missing contrast file: ${contrasts}"
// contrasts

/* 
 * Channel for bam files
 */

bam_files = Channel
          .fromPath(params.bam)
          .map { file -> [ id:file.baseName,file:file] }

/* 
 * SORT BAM
 */

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

/*
 * BAM TO FASTQ
 */

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

/*
 * COPY CHANNEL
 */

fastqs.into { fastqs_kallisto; fastqs_star }

/*
 * KALLISTO INDEX IF NEEDED
 */

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

/*
 *  KALLIST QUANT
 */

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

/*
 *  COMBINE KALLISTO OUTPUT
 */

kallisto_dirs.into{kallisto_dirs; kallisto_dirs_deseq2}

process kallistoCountMatrix {
	tag "anno: ${params.anno_set}"
//	publishDir "$params.output/kallisto_data" , mode: 'copy'	

	input:
	file 'kallisto/*' from kallisto_dirs.collect() 
	

	output:
	file 'kallisto_counts.tab' 
	file 'kallistoData.Rdata' //into kallistodata 	

	script:
	"""
	sumkallisto.R kallisto ${txdb} ${design}
	"""
}

/* 
 * STAR INDEX IF NEEDED
 */
	
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

/*
 * STAR ALIGN
 */

process STAR {
	tag "star: $name"

   	input:
    	file index from star_index
    	set name, file(fq) from fastqs_star
    
    	output:
	set name, file("star_${name}/${name}Aligned.sortedByCoord.out.bam") into sort_bam    
	file("star_${name}/${name}Log.final.out") into final_log    
    	file "star_${name}/${name}ReadsPerGene.out.tab" into starcount

    	script:
    	"""
	mkdir -p star_${name}
    	STAR --genomeDir $index --outFileNamePrefix ./star_${name}/${name} --readFilesIn  $fq --runThreadN 4 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate 
    	"""
}

process STAR_log {
 	publishDir "$params.output/star_data" , mode: 'copy'
	input:
	file 'logs/*' from final_log.collect()

	output:
	file "star_stats.tab" into stats	
	file "star_stats.pdf"
	script:
	"""
	bash star_stats.sh
	plot_star_stats.R ${design}
	"""
}

/*
 * COMBINE STAR COUNTS
 */ 

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
	sumstar.R star ${params.strand} ${design} 
	"""
}

/*
 * BAM 2 BW
 */

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
   	bamCoverage -b ${bam} -o ${name}.bw --normalizeTo1x ${params.normtosize} --binSize=${params.binsize}	
	"""
}




process deseq2 {
publishDir "$params.output/deseq", mode: 'copy'

  input:
  file 'kallisto/*' from kallisto_dirs_deseq2.collect()
  file design
  file contrasts
  
  output:
  file 'pairs.pdf'
  file 'dds.Rdata'
  file 'pca.pdf'
  file 'maplot_*'
  file 'contrast_*'

  script:
  """
  deseq2.R kallisto ${design} ${contrasts} ${params.pvalue}
  """
}

process report {
publishDir "$params.output/deseq", mode: 'copy'

        input:
	file stats from stats
	file report	

	output:
 	file 'report.html'

	script:
 	"""
	createReport.R 1 ${design} ${params.pvalue} ${stats} ${report}
        """
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

  	
	
