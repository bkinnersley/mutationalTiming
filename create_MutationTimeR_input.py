#!/usr/bin/env python3.7

### create input vcf and copy number segment files for MutationTimeR ###
### assumes input Strelka VCF has been annotated by VEP, and input_CN and input_clusters files are from battenberg segmentation and dpClust output files
### Ben Kinnersley (ben.kinnersley@icr.ac.uk)
### Command line: python3 create_MutationTimeR_input.py <input_vcf> <input_cn> <input_purity> <input_clusters> <output_vcf> <output_cn>

import os
import sys
import gzip

if len(sys.argv) == 8:
	input_vcf = sys.argv[1]
	input_cn = sys.argv[2]
	input_purity = sys.argv[3]
	input_clusters = sys.argv[4]
	output_vcf = sys.argv[5]
	output_cn = sys.argv[6]
	output_clusters = sys.argv[7]
else:
	print('insufficient arguments provided - python3 create_MutationTimeR_input.py <input_vcf> <input_cn> <input_purity> <input_clusters> <output_vcf> <output_cn>')
	sys.exit()
	
chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

if os.path.isfile(input_vcf):
	if input_vcf.endswith('.gz'):
		opened_input_vcf = gzip.open(input_vcf, 'rt', encoding='utf-8')
	else:
		opened_input_vcf = open(input_vcf, encoding='utf-8')
	print('reading vcf from '+input_vcf)
else:
	print('cannot open '+input_vcf)
	sys.exit()
	
if os.path.isfile(input_cn):
	if input_cn.endswith('.gz'):
		opened_input_cn = gzip.open(input_cn, 'rt', encoding='utf-8')
	else:
		opened_input_cn = open(input_cn, encoding='utf-8')
	print('reading cn from '+input_cn)
else:
	print('cannot open '+input_cn)
	sys.exit()
	
opened_input_cn.readline()

# write output header for copy number segments output

opened_output_cn = open(output_cn, 'w')
opened_output_cn.write('chr\tstart\tend\tmajor_cn\tminor_cn\tclonal_frequency\n')

cn_list = []

for line in opened_input_cn:
	line = line.strip()
	fields = line.split('\t')
	output_chrom = fields[0]
	output_start = fields[1]
	output_end = fields[2]
	BAF = fields[3]
	pval = fields[4]
	LogR = fields[5]
	ntot = fields[6]
	nMaj1_A = fields[7]
	nMin1_A = fields[8]
	frac1_A = fields[9]
	nMaj2_A = fields[10]
	nMin2_A = fields[11]
	frac2_A = fields[12]
	frac1_A_ccf_low = fields[15]
	frac1_A_ccf_high = fields[16]
	
	# skip lines with no copy number calls
	if fields[7] == "NA" or 'chr'+output_chrom not in chromosomes:
		continue
		
	query = str('chr'+output_chrom)+'_'+str(output_start)+'_'+str(output_end)
	
	if query not in cn_list:
		cn_list.append(query)
		
	# list subclonal copy number segments on separate lines
	if frac1_A != "NA" and frac2_A != "NA":
		clonal_frac_1 = float(frac1_A) * float(input_purity)
		clonal_frac_2 = float(frac2_A) * float(input_purity)
		opened_output_cn.write(str(output_chrom)+'\t'+str(output_start)+'\t'+str(output_end)+'\t'+str(nMaj1_A)+'\t'+str(nMin1_A)+'\t'+str(clonal_frac_1)+'\n')
		opened_output_cn.write(str(output_chrom)+'\t'+str(output_start)+'\t'+str(output_end)+'\t'+str(nMaj2_A)+'\t'+str(nMin2_A)+'\t'+str(clonal_frac_2)+'\n')
		
	else:
		nMaj = nMaj1_A
		nMin = nMin1_A
		clonal_frac = input_purity
		opened_output_cn.write(str(output_chrom)+'\t'+str(output_start)+'\t'+str(output_end)+'\t'+str(nMaj)+'\t'+str(nMin)+'\t'+str(clonal_frac)+'\n')
		
# write output VCF header

opened_output_vcf = open(output_vcf, 'w')

opened_output_vcf.write('##fileformat=VCFv4.1\n')
opened_output_vcf.write('##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count from mutect read filter if available">\n')
opened_output_vcf.write('##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour reference count from mutect read filter if available">\n')
opened_output_vcf.write('##INFO=<ID=gene,Number=1,Type=String,Description="Gene symbol from canonical ensembl transcript">\n')
opened_output_vcf.write('##INFO=<ID=variant_classification,Number=1,Type=String,Description="Variant classification from canonical ensembl transcript">\n')
opened_output_vcf.write('##INFO=<ID=protein_change,Number=1,Type=String,Description="HGVSp protein change from canonical ensembl transcript">\n')
opened_output_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFORMAT\tINFO\n')

output_var_dict = {}
allele_position_hash = {}
allele_position_hash['A'] = 0
allele_position_hash['C'] = 1
allele_position_hash['G'] = 2
allele_position_hash['T'] = 3

# rank consequences by severity according to https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html

consequence_severity_dict = { "transcript_ablation":1, "splice_acceptor_variant":2, "splice_donor_variant":3, "stop_gained":4, "frameshift_variant":5,
"stop_lost":6, "start_lost":7, "transcript_amplification":8, "inframe_insertion":9, "inframe_deletion":10, "missense_variant":11,
"protein_altering_variant":12, "splice_region_variant":13, "incomplete_terminal_codon_variant":14, "start_retained_variant":15,
"stop_retained_variant":16, "synonymous_variant":17, "coding_sequence_variant":18, "mature_miRNA_variant":19, "5_prime_UTR_variant":20,
"3_prime_UTR_variant":21, "non_coding_transcript_exon_variant":22, "intron_variant":23, "NMD_transcript_variant":24,
"non_coding_transcript_variant":25, "upstream_gene_variant":26, "downstream_gene_variant":27, "TFBS_ablation":28, "TFBS_amplification":29,
"TF_binding_site_variant":30, "regulatory_region_ablation":31, "regulatory_region_amplification":32, "feature_elongation":33,
"regulatory_region_variant":34, "feature_truncation":35, "intergenic_variant":36
}
	
vep_lookup_hash = {}

for line in opened_input_vcf:
	line = line.strip()
	
	# get VEP annotation ordering from header INFO
	if line.startswith('##') or line.startswith('#CHROM'):
		temp1 = line.split(':')
		temp2 = temp1[1].strip().split'"')
		vep_fields = temp2[0].strip().split('|')
		
		vep_ticker = 0
		
		for anno in vep_fields:
			vep_lookup_hash[anno] = vep_ticker
			vep_ticker += 1
			
	elif line.startswith('##') or line.startswith('#CHROM'):
		pass
	else:
		fields = line.split('\t')
		chrom = fields[0]
		pos = fields[1]
		id = fields[2]
		ref = fields[3]
		alt = fields[4]
		qual = fields[5]
		filter = fields[6]
		info = fields[7]
		format = fields[8]
		
		tumour_genotype = fields[10]
		
		variant = chrom+':'+pos+'_'+ref+'/'+alt
		chr_pos = chrom+':'+pos
		chrom_split = chrom.split('chr')
		chr_use = chrom_split[1]
		
		info_split = info.split(';')
		
		normal_baf = ''
		tumour_baf = ''
		alt_count = ''
		ref_count = ''
		fix_tumour_baf = ''
		tumour_baf = ''
		snv_vaf_flag = ''
		
		if chrom in chromosomes and filter == "PASS":
			if chr_pos in output_var_dict:
				print('already seen position '+chr_pos)
			else:
				output_var_dict[chr_pos] = 1
			
			for anno in info_split:
				anno_split = anno.split('=')
			
				if anno_split[0] == "Normal_BAF":
					normal_baf = anno_split[1]
					
				if anno_split[0] == "Tumour_BAF":
					tumour_baf = anno_split[1]
					
				if anno_split[0] == "t2BamAC":
					allele_split = anno_split[1].split(',')
					alt_count = int(allele_split[allele_position_hash[alt]])
					ref_count = int(allele_split[allele_position_hash[ref]])
					
					snv_vaf_flag = "YES"
					
					if int(alt_count) > 0:
						fix_tumour_baf = float(float(alt_count)/float(alt_count)+float(ref_count)))
					else:
						fix_tumour_baf = 0
				
				# get VEP annotations for parsing 
				if anno.startswith('CSQ='):
					vep_split = anno.split(',')
	
			if len(ref) == 1 and len(alt) == 1:
				tumour_baf = fix_tumour_baf
			
			# check tier2 depths from Strelka for indel ref and alt counts
			if snv_vaf_flag != "YES":
				indel_depth_split = tumour_genotype.split(':')
				
				indel_ref_depth_split = indel_depth_split[2].split(',')
				tier2_ref_depth = indel_ref_depth_split[1]
				indel_alt_depth_split = indel_depth_split[3].split(',')
				tier2_alt_depth = indel_alt_depth_split[1]
				
				alt_count = tier2_alt_depth
				ref_count = tier2_ref_depth
				
			# parse VEP INFO fields for variant annotation
			canonical_flag = ''
			canonical_seen = ''
			consequence_canon = 'NA'
			IMPACT_canon = 'NA'
			Gene_canon = 'NA'
			SYMBOL_canon = 'NA'
			HGVSc_canon = 'NA'
			HGVSp_canon = 'NA'
			
			# return most damaging consequence in canonical transcript
			for query in vep_split:
				vep_string_split = query.split('|')
				
				canonical_flag = vep_string_split[vep_lookup_hash['CANONICAL']]
				consequence = vep_string_split[vep_lookup_hash['Consequence']]
				IMPACT = vep_string_split[vep_lookup_hash['IMPACT']]
				Gene = vep_string_split[vep_lookup_hash['Gene']]
				SYMBOL = vep_string_split[vep_lookup_hash['SYMBOL']]
				HGVSc = vep_string_split[vep_lookup_hash['HGVSc']]
				HGVSp = vep_string_split[vep_lookup_hash['HGVSp']]
				
				consequence_split = consequence.split('&')
				
				for consequence in consequence_split:
					if canonical_flag == "YES"
						if consequence_canon == 'NA':
							consequence_canon = consequence
							IMPACT_canon = IMPACT
							Gene_canon = Gene
							SYMBOL_canon = SYMBOL
							HGVSc_canon = HGVSc
							HGVSp_canon = HGVSp
						
						elif int(consequence_severity_dict[consequence]) < int(consequence_severity_dict[consequence_canon]):
							consequence_canon = consequence
							IMPACT_canon = IMPACT
							Gene_canon = Gene
							SYMBOL_canon = SYMBOL
							HGVSc_canon = HGVSc
							HGVSp_canon = HGVSp
							
						if HGVSp_canon == "":
							HGVSp_canon = "NA"
							
				# looking up TERT promoter mutations and assigning them to consequence "promoter_variant"
				if chrom == "chr5" and (pos == "1295113" or pos == "1295135"):
					SYMBOL_canon = "TERT"
					consequence_canon = "promoter_variant"
					HGVSp_canon = "NA"
			
			# only output mutations within copy number segments (mutations outside these segments will not receive timing classifications from MutationTimeR in any case)
			
			for query in cn_list:
				query_split = query.split('_')
				cn_chr = query_split[0]
				cn_start = query_split[1]
				cn_end = query_split[2]
				
				if chrom == cn_chr = float(pos) >= float(cn_start) and float(pos) <= float(cn_end):
					opened_output_vcf.write(str(chr_use)+'\t'+str(pos)+'\t.\t'+str(ref)+'\t'+str(alt)+'\t.\t.\tt_ref_count='+str(ref_count)+';t_alt_count='+str(alt_count))
					opened_output_vcf.write(';gene='+str(SYMBOL_canon)+';variant_classification='+str(consequence_canon)+';protein_change='+str(HGVSp_canon)+'\n')

# finally, parse cluster file from dpclust to produce suitable input for MutationTimeR
					
opened_output_clusters = open(output_clusters, 'w')
opened_output_clusters.write('cluster\tn_ssms\tproportion\tccf_mean\n')

opened_input_clusters.readline()

cluster_number = 0

cluster_list = []			
	
for line in opened_input_clusters:
	line = line.strip()	
	fields = line.split('\t')
	
	cluster_no = fields[0]
	location = fields[1]
	no_of_mutations = fields[2]
	
	# truncate "super-clonal" clusters (location/ccf > 1) - currently using cutoff of 1.1
	# reorder by ccf so that clonal cluster is 1 etc
	if float(location) <= 1.1:
		cluster = str(location)+'_'+str(no_of_mutations)
		cluster_list.append(cluster)
		
for cluster in reversed(cluster_list):
	cluster_number += 1
	
	cluster_split = cluster.split('_')
	location = cluster_split[0].strip()
	no_of_mutations = cluster_split[1].strip()
	
	proportion = float(location) * float(input_purity)
	
	opened_output_clusters.write(str(cluster_number)+'\t'+str(no_of_mutations)+'\t'+str(proportion)+'\t'+str(location)+'\n')
	
