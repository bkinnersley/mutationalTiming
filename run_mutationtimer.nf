#!/usr/bin/env nextflow

/*
===================================================
								run_mutationtimer
===================================================
Ben Kinnersley (ben.kinnersley@icr.ac.uk)
------------------------------------------------
run_mutationtimer:
	An analysis pipeline to create input files and carry out clonality/subclonality and mutational timing analysis using MutationTimeR (https://github.com/gerstung-lab/MutationTimeR)
------------------------------------------------
*/

def helpMessage() {
	log.info nfcoreHeader()
	log.info"""
	Usage:
	The typical command for running the pipeline is as follows:
		nextflow run run_mutationtimer.nf --input_samples <INPUT_SAMPLES> --cohort_name <COHORT_NAME> --output_dir /path/to/output -c <CONF_FILE>
		
	Mandatory arguments:
		--input_samples		Tab-separated text file containing all samples to be analysed (<PARTICIPANT_ID> <PLATE_TUMOUR> <PLATE_GERMLINE> <SEX> <PURITY> <VCF_PATH> <CN_PATH> <DPCLUST_PATH>)
									Assumes VEP annotations are in INFO fields of VCF. All chr1-22,X mutations with "PASS" filter overlapping copy number segments will be included
		--cohort_name		String to name_cohort to form basis of output files and directories
		--output_dir			The output directory where the results will be saved
		""".stripIndent()
}

// Show help message
if (params.help) exit 0, helpMessage()

params.lsfProject="re_gecip_cancer_colorectal"
input_samples = file(params.input_samples)
cohort_name = params.cohort_name
output_dir = params.output_dir

initial_samples = Channel
	.from(input_samples)
	.splitCsv(sep: '\t')
	.map { row ->
			def participant_id = row[0]
			def plate_tumour = row[1]
			def plate_germline = row[2]
			def sex = row[3]
			def purity = row[4]
			def wgd = row[5]
			def vcf_file = file(row[6])
			def cn_file = file(row[7])
			def dpclust_file = file(row[8])
			return [ participant_id, plate_tumour, plate_germline, sex, purity, wgd, vcf_file, cn_file, dpclust_file ]
	}
	.set { INITIAL_SAMPLE_FILES }

// in case we want to branch into other processes	
INITIAL_SAMPLE_FILES.into { input_mutationtimer }

process create_mutationtimer_input{
	tag "create mutationtimer input ${plate_tumour}"
	
	publishDir "${params.output_dir}/per_sample/${participant_id}", mode:'copy'
	
	input:
		tuple val(participant_id), val(plate_tumour), val(plate_germline), val(sex), val(purity), val(wgd), file(vcf_file), file(cn_file), file(dpclust_file) from input_mutationtimer
		
	output:
		tuple file("$plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf.gz"), file("${plate_tumour}.GRCh38.PASS.MutationTimeR.segments.tsv"),
			file("${plate_tumour}.GRCh38.PASS.MutationTimeR.clusters.tsv"), val(participant_id), val(plate_tumour), val(plate_germline), val(sex), val(purity), val(wgd) into mutationtimer_mut
		
	beforeScript 'module purge'

	script:
	"""
	#!/bin/bash
	
	FILE_VCF-${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf.gz
	FILE_CNV=${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.segments.tsv
	FILE_CLUSTERS=${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.cluster.tsv
	
	if [ -f "\$FILE_VCF" ] && [ -f "\$FILE_CNV" ] && [ -f "\$FILE_CLUSTERS" ]; then
		ln -s \$FILE_VCF
		ln -s \$FILE_CNV
		ln -s \$FILE_CLUSTERS
		
	else
	
	python3 ${params.mutationtimer_input_script} \
	${vcf_file} \
	${cn_file} \
	${purity} \
	${dpclust_file} \
	${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf \
	${plate_tumour}.GRCh38.PASS.MutationTimeR.segments.tsv \
	${plate_tumour}.GRCh38.PASS.MutationTimeR.clusters.tsv
	
	gzip ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf
	
	fi
	"""
}

process run_mutationtimer{
	tag "run mutationtimer ${plate_tumour}"
	
	publishDir "${params.output_dir}/per_sample/${participant_id}", mode:'copy'
	
	input:
		tuple file("$plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf.gz"), file("${plate_tumour}.GRCh38.PASS.MutationTimeR.segments.tsv"),
			file("${plate_tumour}.GRCh38.PASS.MutationTimeR.clusters.tsv"), val(participant_id), val(plate_tumour), val(plate_germline), val(sex), val(purity), val(wgd) from mutationtimer_mut
	
	output:
		tuple file("${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf.gz"), file("${plate_tumour}_MutationTimeR.pdf"), val(participant_id), val(plate_tumour), val(plate_germline),
			val(sex), val(purity), val(wgd) into mutationtimer_extract_driver_clonal_state
		file("${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf.gz.tbi")
		file("${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.rds")
		file("${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.CNtiming.tsv")
		
	beforeScript 'module purge; module load lang/R/3.6.2-foss-2019b bio/BCFtools/1.11-GCC-8.3.0 bio/tabix/0.2.6-GCCcore-7.3.0'
	
	script:
	"""
	#!Rscript
	
	output_pdf = '${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}_MutationTimeR.pdf'
	output_vcf = '${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf.gz'
	output_vcf_tbi = '${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf.gz.tbi'
	output_rds = '${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.rds'
	output_bb = '${params.output_dir}/${params.cohort_name}/MutationTimeR/per_sample/${participant_id}/${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.CNtiming.tsv'
	
	if(file.exists(output_pdf) & file.exists(output_vcf) & file.exists(output_vcf_tbi) & file.exists(output_rds) & file.exists(output_bb)) {
		file.symlink(output_pdf, ".")
		file.symlink(output_vcf, ".")
		file.symlink(output_vcf_tbi, ".")
		file.symlink(output_rds, ".")
		file.symlink(output_bb, ".")
		
	} else {
	
	.libPaths("${params.R_library_3_6}")
	
	library(MutationTimeR)
	
	system('zcat ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf.gz > ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf')
	
	if("${sex}"=="MALE") {
		mt_gender = 'male'
	} else {
		mt_gender = 'female'
	}
	
	input_vcf = '${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.vcf'
	input_cn = '${plate_tumour}.GRCh38.PASS.MutationTimeR.segments.tsv'
	input_clusters = '${plate_tumour}.GRCh38.PASS.MutationTimeR.clusters.tsv'
	
	cn <- read.delim(file=input_cn, sep="\\t", header=TRUE)
	
	bb <- makeGRangesFromDataFrame(cn,
							keep.extra.columns=TRUE,
							ignore.strand=TRUE,
							seqinfo=NULL,
							seqnames.field=c("seqnames", "seqname",
												"chromosome", "chrom",
												"chr", "chromosome_name", "seqid"),
							start.field="start"
							end.field=c("end", "stop"),
							strand.field="strand",
							starts.in.df.are.0based=FALSE)
							
	vcf <- readVcf(input_vcf)
	
	clusters <- read.delim(file=input_clusters, sep="\\t", header=TRUE)
	
	mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=1000, gender=mt_gender, isWgd=${wgd})
	
	info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
	info(vcf) <- cbind(info(vcf), mt\$V)
	
	mcols(bb) <- cbind(mcols(bb),mt\$T)
	
	bb_file <- '${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.CNtiming.tsv'
	
	bb_df <- as.data.frame(bb)
	bb_df <- subset(bb_df, select=-c(timing_param))
	
	write.table(bb_df, file=bb_file, append=FALSE, quote=FALSE, sep="\\t", eol="\\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE, qmethod="escape", fileEncoding="")
	
	pdf('${plate_tumour}_MutationTimeR.pdf')
	plotSample(vcf,bb)
	dev.off()
	
	vcf_out = '${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf'
	
	writeVcf(vcf, vcf_out)
	
	rds_file = '${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.rds'
	
	saveRDS(mt, file=rds_file)
	
	system('bgzip ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf')
	system('tabix -p vcf ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf.gz')
	system('rm ${plate_tumour}.GRCh38.PASS.MutationTimeR.SNV.INDEL.anno.vcf')
	
	}
	
	"""	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	