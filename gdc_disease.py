'''
Author: Rosina Savisaar and Liam Abrahams
'''

import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import generic as gen
import gdc_ops as gdco
import os
import random
import SNP_ops as so
import time
import numpy as np
import shutil

def main():

    description = "Analysis of disease mutations from GDC."
    args = gen.parse_arguments(description,  ["output_dir", "mutations_dir", "results_prefix", "full_run", "snp_run", "process_mutations", "intersect_mutations", "get_relative_positions", "get_snp_status"], flags = [3,4,5,6,7,8])
    output_dir, mutations_dir, results_prefix, full_run, snp_run, process_mutations, intersect_mutations, get_relative_positions, get_snp_status = args.output_dir, args.mutations_dir, args.results_prefix, args.full_run, args.snp_run, args.process_mutations, args.intersect_mutations, args.get_relative_positions, args.get_snp_status

    start = time.time()
    # create the output directories
    gen.create_output_directories('temp_data/')
    gen.create_output_directories(output_dir)

    # prcoess the mutation files
    processed_mutations = "{0}/processed_mutations.vcf".format(output_dir)
    if process_mutations or full_run or snp_run:
        print("Processing mutation files...")
        gdco.process_mutation_files(mutations_dir, output_dir, processed_mutations)

    # intersec the snps with coding exons
    processed_mutations_gz = "{0}.gz".format(processed_mutations)
    coding_exons_file = "{0}_coding_exons.bed".format(results_prefix)
    intersect_file = "{0}/processed_mutations_coding_exons_intersect.vcf".format(output_dir)
    bed_intersect_file = "{0}/processed_mutations_coding_exons_intersect.bed".format(output_dir)
    if intersect_mutations or full_run or snp_run:
        print("Intersecting mutations with coding exons...")
        gen.run_process(["tabix", "-p", "vcf", processed_mutations_gz])
        so.intersect_snps_parallel(coding_exons_file, processed_mutations_gz, intersect_file)
        print("Converting to bed...")
        so.intersect_vcf_to_bed(coding_exons_file, intersect_file, bed_intersect_file, change_names = True)

    # get relative positions of the snps in cds and exons
    full_cds_bed = "{0}_CDS.bed".format(results_prefix)
    relative_exon_positions = "{0}/snp_relative_exon_positions.bed".format(output_dir)
    relative_cds_positions = "{0}/snp_relative_cds_positions.bed".format(output_dir)
    if get_relative_positions or full_run or snp_run:
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(bed_intersect_file, relative_exon_positions, broad_snps_shift=True)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, relative_cds_positions, full_cds_bed, broad_snps_shift=True)

    # get the change status of the snps to check them, generate two files
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    ptc_file = "{0}/gdc_ptcs.txt".format(output_dir)
    other_file = "{0}/gdc_other_snps.txt".format(output_dir)
    if get_snp_status or full_run or snp_run:
        print("Getting snp status...")
        so.get_snp_change_status(relative_cds_positions, cds_fasta, ptc_file, other_file, broad_snps_shift=True)


if __name__ == "__main__":
    main()
