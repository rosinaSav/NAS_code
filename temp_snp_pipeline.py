'''
Author: Liam Abrahams
'''

import bed_ops as bed
import bam_ops as bam
import SNP_ops as snp
import generic as gen
import re
import numpy as np
import shutil
import os

#i have used this just to get the top 2000 lines of the gtf as dont need all for testing
gen.extract_head_of_file("./source_data/gtfs/Homo_sapiens.GRCh37.87.gtf", 2000)

#extract cds
# returns: bed file containing coordinates, fasta file containing intervals, fasta file containing whole cds
def extract_cds():
	# I have just used extracted files here to make them more usable
	gtf_file = "./source_data/gtfs/Homo_sapiens.GRCh37.87.extracted.2000.gtf"
	genome_fasta = "./source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	output_bed_file = "./temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds.bed"
	output_cds_fasta = "./temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds.fasta"
	#create temp directory if it doesnt already exist
	gen.create_directory("./temp_files/")
	bed.extract_cds(gtf_file, output_bed_file, output_cds_fasta, genome_fasta, all_checks=True)

#extract snps
def extract_snps():
	#want to create a bed file with the full chromosome names for intersecting snps
	gtf_file = "./source_data/gtfs/Homo_sapiens.GRCh37.87.extracted.2000.gtf"
	cds_full_chr_name_name = "./temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds.full_chr_name.bed"
	bed.extract_features(gtf_file, cds_full_chr_name_name, ["CDS", "stop_codon"], full_chr_name=True)

	#ensure this only contains the featuers that passed filtering
	cds_fasta = "./temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds.fasta"
	bed.filter_bed_from_fasta(cds_full_chr_name_name, cds_fasta, cds_full_chr_name_name)

	#extact snps that overlap cds features
	vcf_file = "./source_data/per_sample_vcfs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
	snps_bed = "./temp_files/Homo_sapiens.GRCh37.87.extracted.snps.bed"
	snp.tabix(cds_full_chr_name_name, snps_bed, vcf_file)

	#create an intersect file
	intersect_file = "./temp_files/cds_snp_intersect.bed"
	bam.intersect_bed(cds_full_chr_name_name, snps_bed, output_file=intersect_file, write_both=True,  no_dups=False)

def get_snp_positions():
	#get the positions of the snps relative to the cds
	intersect_file = "./temp_files/cds_snp_intersect.bed"
	output_file = "./temp_files/snp_relative_cds_positions.bed"
	cds_intervals = "temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds_intervals.fasta"
	snp_exon_realtive_positions = snp.get_snp_relative_exon_position(intersect_file)
	snp.get_snp_relative_cds_position(snp_exon_realtive_positions, output_file, cds_intervals)

def get_snp_change_status():
	#get what type of mutation the snp causes
	snp_relative_cds_positions = "./temp_files/snp_relative_cds_positions.bed"
	cds_fasta = "./temp_files/Homo_sapiens.GRCh37.87.extracted.2000.cds.fasta"
	ptc_outputs = "./temp_files/ptc_causing_snps.bed"
	other_outputs = "./temp_files/other_causing_snps.bed"
	snp.get_snp_change_status(snp_relative_cds_positions, cds_fasta, ptc_outputs, other_outputs)


extract_cds()
extract_snps()
get_snp_positions()
get_snp_change_status()
