'''
Author: Rosina Savisaar and Liam Abrahams
Check whether PTCs are associated with greater rates of exon skipping.
'''

import bed_ops as beo
import bam_ops as bao
import disease_ops as do
import SNP_ops as so
import generic as gen
import os
import random
import SNP_ops as so
import time
import numpy as np
import shutil

def main():

    description = "Analysis of disease mutations from firehose broad."
    args = gen.parse_arguments(description,  ["output_dir", "mutations_dir", "rna_dir", "exon_reads_dir", "results_prefix", "subset", "subset_no", "clean_run", "snp_ops", "process_mutations", "intersect_mutations", "convert_intersect", "get_relative_positions", "get_snp_status", "filter_junctions", "process_exon_reads", "process_reads", "check_ptcs"], flags = [5,7,8,9,10,11,12,13,14,15,16,17], ints = [6])
    output_dir, mutations_dir, rna_dir, exon_reads_dir, results_prefix, subset, subset_no, clean_run, snp_ops, process_mutations, intersect_mutations, convert_intersect, get_relative_positions, get_snp_status, filter_junctions, process_exon_reads, process_reads, check_ptcs = args.output_dir, args.mutations_dir, args.rna_dir, args.exon_reads_dir, args.results_prefix, args.subset, args.subset_no, args.clean_run, args.snp_ops, args.process_mutations, args.intersect_mutations, args.convert_intersect, args.get_relative_positions, args.get_snp_status, args.filter_junctions, args.process_exon_reads, args.process_reads, args.check_ptcs

    start = time.time()
    gen.create_directory('temp_data/')
    if clean_run:
        shutil.rmtree(output_dir)
    gen.create_output_directories(output_dir)

    # create folder containing processed mutation files
    processed_dir = "{0}/processed_mutation_files".format(output_dir)
    filename_prefix = "processed_mutations"
    full_mutation_file = "{0}/mutations.vcf".format(output_dir)
    gen.create_output_directories(processed_dir)

    if process_mutations or snp_ops or clean_run:
        print("Processing mutation files...")
        entrylimit = 10000
        do.refactor_files(mutations_dir, processed_dir, filename_prefix, full_mutation_file, entrylimit, subset, subset_no, clean_directory = True)

    full_mutation_file_zip = "{0}.gz".format(full_mutation_file)
    coding_exons_file = "{0}_coding_exons.bed".format(results_prefix)
    intersect_file = "{0}/mutations_coding_exons_intersect.vcf".format(output_dir)
    if intersect_mutations or snp_ops or clean_run:
        print("Intersecting mutations with coding exons...")
        gen.run_process(["tabix", "-p", "vcf", full_mutation_file_zip])
        so.intersect_snps_parallel(coding_exons_file, full_mutation_file_zip, intersect_file)

    bed_intersect_file = "{0}/mutations_coding_exons_intersect.bed".format(output_dir)
    if convert_intersect or snp_ops or clean_run or not os.path.isfile(bed_intersect_file):
        print("Converting intersect file to bed format...")
        so.intersect_vcf_to_bed(coding_exons_file, intersect_file, bed_intersect_file, change_names = True)

    # get relative positions of the snps in cds and exons
    full_bed = "{0}_CDS.bed".format(results_prefix)
    disease_snps_relative_exon_positions = "{0}/disease_snp_relative_exon_positions.bed".format(output_dir)
    disease_snps_relative_cds_positions = "{0}/disease_snp_relative_cds_positions.bed".format(output_dir)
    if get_relative_positions or snp_ops or clean_run or not os.path.isfile(disease_snps_relative_exon_positions) or not os.path.isfile(disease_snps_relative_cds_positions):
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(bed_intersect_file, disease_snps_relative_exon_positions, broad_snps_shift=True)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, disease_snps_relative_cds_positions, full_bed, broad_snps_shift=True)

    # get the change status of the snps to check them, generate two files
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    disease_ptcs_file = "{0}/disease_ptcs.txt".format(output_dir)
    disease_other_file = "{0}/disease_other_snps.txt".format(output_dir)
    if get_snp_status or snp_ops or clean_run or not os.path.isfile(disease_ptcs_file) or not os.path.isfile(disease_other_file):
        print("Getting snp status...")
        so.get_snp_change_status(disease_snps_relative_cds_positions, cds_fasta, disease_ptcs_file, disease_other_file, broad_snps_shift=True)

    exon_junctions_file = "{0}_exon_junctions.bed".format(results_prefix)
    ptc_exon_junctions_file = "{0}/ptc_filtered_exon_junctions.bed".format(output_dir)
    if filter_junctions or clean_run or not os.path.isfile(ptc_exon_junctions_file):
        #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
        print("Filtering exon-exon junctions to only leave those that flank exons with a PTC variant...")
        beo.filter_exon_junctions(exon_junctions_file, disease_ptcs_file, ptc_exon_junctions_file)

    # junctions = bao.read_exon_junctions(ptc_exon_junctions_file)

    exon_reads_output_dir = "{0}/processed_exon_reads".format(output_dir)
    if process_exon_reads:
        do.process_raw_reads(exon_reads_dir, exon_reads_output_dir)

    processed_rna_dir = "{0}/processed_reads".format(output_dir)
    processed_junction_suffix = "processed_junctions"
    if process_reads:
        do.process_counts(rna_dir, processed_rna_dir, processed_junction_suffix, exon_junctions_file, junctions, results_prefix)


    if check_ptcs:
        do.check_ptcs(disease_ptcs_file, processed_rna_dir, processed_junction_suffix)

if __name__ == "__main__":
    main()
