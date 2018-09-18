'''
Auhor: Liam Abrahams.
Look at disease SNPs from ClinVar.
'''

import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import disease_snps_ops as dso
import os
import collections
import re
import copy
import numpy as np
import random

def main():

    description = "Look at disease snps."
    arguments = ["disease_snps_file", "output_directory", "results_prefix", "simulations", "ese_file", "intersect_snps", "get_relative_positions", "get_snp_status", "get_info", "simulate_ptc_location", "get_possible_ptc_locations", "required_simulations", "get_overlaps", "intersect_ptcs", "compare_ptcs" ,"get_introns", "compare_distances", "clinvar_ptc_locations", "location_simulation", "exclude_cpg", "ese_hit_simulation", "only_disease", "only_kgenomes", "only_ese", "get_unique_ptcs", "get_unique_rel_pos", "excess_test", "disease_locations_chisquare"]
    args = gen.parse_arguments(description, arguments, flags = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, 21, 22, 23,24,25,26,27], ints=[3])
    disease_snps_file, output_directory, results_prefix, simulations, ese_file, intersect_snps, get_relative_positions, get_snp_status, get_info, simulate_ptc_location, get_possible_ptc_locations, required_simulations, get_overlaps, intersect_ptcs, compare_ptcs, get_introns, compare_distances, clinvar_ptc_locations, location_simulation, exclude_cpg, ese_hit_simulation, only_disease, only_kgenomes, only_ese, get_unique_ptcs, get_unique_rel_pos, excess_test, disease_locations_chisquare = args.disease_snps_file, args.output_directory, args.results_prefix, args.simulations, args.ese_file, args.intersect_snps, args.get_relative_positions, args.get_snp_status, args.get_info, args.simulate_ptc_location, args.get_possible_ptc_locations, args.required_simulations, args.get_overlaps, args.intersect_ptcs, args.compare_ptcs, args.get_introns, args.compare_distances, args.clinvar_ptc_locations, args.location_simulation, args.exclude_cpg, args.ese_hit_simulation, args.only_disease, args.only_kgenomes, args.only_ese, args.get_unique_ptcs, args.get_unique_rel_pos, args.excess_test, args.disease_locations_chisquare

    if simulations and not isinstance(simulations, int):
        print("\nERROR: Please provide the correct number for simulations.\n")
        raise Exception

    # create the output directory if it doesnt already exist
    gen.create_output_directories(output_directory)

    # disease_snps_file = "./source_data/clinvar_20180429.vcf.gz"
    disease_snps_index_file = "{0}.tbi".format(disease_snps_file)

    if not os.path.isfile(disease_snps_file) or not os.path.isfile(disease_snps_index_file):
        print("\nERROR: Please provide the required disease SNPs file(s).\n")
        raise Exception

    # intersect the coding exons with the disease snps
    exon_bed = "{0}_coding_exons.bed".format(results_prefix)
    disease_snp_intersect_file_vcf = "{0}/disease_snp_intersect.vcf".format(output_directory)
    disease_snp_intersect_file_bed = "{0}/disease_snp_intersect.bed".format(output_directory)
    if intersect_snps:
        print("Intersecting snps with exons")
        so.intersect_snps_parallel(exon_bed, disease_snps_file, disease_snp_intersect_file_vcf)
        so.intersect_vcf_to_bed(exon_bed, disease_snp_intersect_file_vcf, disease_snp_intersect_file_bed, change_names = True)

    # get relative positions of the snps in cds and exons
    full_bed = "{0}_CDS.bed".format(results_prefix)
    disease_snps_relative_exon_positions = "{0}/disease_snp_relative_exon_positions.bed".format(output_directory)
    disease_snps_relative_cds_positions = "{0}/disease_snp_relative_cds_positions.bed".format(output_directory)
    if get_relative_positions:
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(disease_snp_intersect_file_bed, disease_snps_relative_exon_positions)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, disease_snps_relative_cds_positions, full_bed)

    # get the change status of the snps to check them
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    disease_ptcs_file = "{0}/disease_ptcs.txt".format(output_directory)
    disease_other_file = "{0}/disease_other_snps.txt".format(output_directory)
    if get_snp_status:
        print("Getting snp status...")
        so.get_snp_change_status(disease_snps_relative_cds_positions, cds_fasta, disease_ptcs_file, disease_other_file)

    # get intersect between the clinvar ptcs and 1000 genomes ptcs
    ptc_file = "{0}_ptc_file.txt".format(results_prefix)
    ptc_intersect_file = "{0}/ptc_intersect.bed".format(output_directory)
    if intersect_ptcs:
        temp_disease_ptc_file = "temp_data/{0}".format(random.random())
        dso.refactor_ptc_file(disease_ptcs_file, temp_disease_ptc_file)
        temp_k_genomes_ptc_file = "temp_data/{0}".format(random.random())
        dso.refactor_ptc_file(ptc_file, temp_k_genomes_ptc_file, header=True)
        bao.intersect_bed(temp_k_genomes_ptc_file, temp_disease_ptc_file, write_both = True, no_dups=False, output_file = ptc_intersect_file)
        gen.remove_file(temp_disease_ptc_file)
        gen.remove_file(temp_k_genomes_ptc_file)

    # get a list of ptcs unique to each dataset
    unique_ptcs = "{0}/disease_ptcs_no_intersect.bed".format(output_directory)
    unique_ptcs_kgenomes = "{0}/kgenomes_ptcs_no_intersect.bed".format(output_directory)
    if get_unique_ptcs:
        dso.get_unique_ptcs(disease_ptcs_file, ptc_file, ptc_intersect_file, unique_ptcs, unique_ptcs_kgenomes)

    # get the relative positions of the ptcs unique to each dataset
    unique_ptcs_rel_pos_file = "{0}/disease_ptcs_no_intersect_rel_pos.bed".format(output_directory)
    kgenomes_relative_positions = "{0}_PTC_relative_exon_positions.bed".format(results_prefix)
    kgenomes_unique_ptcs_rel_pos_file = "{0}/kgenomes_ptcs_no_intersect_rel_pos.bed".format(output_directory)
    if get_unique_rel_pos:
        dso.get_unique_rel_pos(unique_ptcs, disease_snps_relative_exon_positions, unique_ptcs_kgenomes, kgenomes_relative_positions, unique_ptcs_rel_pos_file, kgenomes_unique_ptcs_rel_pos_file)


    # get the ese file name
    ese_file_name = ese_file.split('/')[-1].split('.')[0]
    # get the coding exons fasta file path
    coding_exons_fasta = "{0}_coding_exons.fasta".format(results_prefix)

    # snp_relative_positions_file = "{0}_SNP_relative_exon_position.bed".format(results_prefix)

    # simulation picking random reference allele matched simulants
    clinvar_location_simulation_file = "{0}/clinvar_ptc_location_simulation.csv".format(output_directory)
    clinvar_location_simulation_ese_overlap_file = "{0}/clinvar_ptc_location_simulation_{1}_ese_overlaps.csv".format(output_directory, ese_file_name)
    kgenomes_location_simulation_file = "{0}/1000_genomes_simulations.csv".format(output_directory)
    kgenomes_location_simulation_ese_overlap_file = "{0}/1000_genomes_simulations_ese_overlaps.csv".format(output_directory)

    if location_simulation:
        if not only_kgenomes:
            print('Running ptc location simulation on disease PTCs...')
            dso.ptc_location_simulation(unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, clinvar_location_simulation_file, clinvar_location_simulation_ese_overlap_file, ese_file, only_ese, exclude_cpg)
        if not only_disease:
            print('Running ptc location simulation on 1000 genomes PTCs...')
            dso.ptc_location_simulation(kgenomes_unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, clinvar_location_simulation_file, clinvar_location_simulation_ese_overlap_file, ese_file, only_ese, exclude_cpg)


    window_start = 3
    window_end = 69
    clinvar_ese_hit_simulation_file = "{0}/clinvar_{1}_ese_hit_simulation_{2}_{3}.csv".format(output_directory, ese_file_name, window_start, window_end)
    kgenomes_ese_hit_simulation_file = "{0}/1000_genomes_ese_hit_simulation_{1}_{2}.csv".format(output_directory, window_start, window_end)

    # do a simulation picking only sites from within the region
    if ese_hit_simulation:
        if not only_kgenomes:
            print("Simulating ESE hits on the {0}-{1} region for disease PTCs...".format(window_start, window_end))
            dso.ese_hit_simulation(unique_ptcs, disease_snps_relative_exon_positions, ptc_file, coding_exons_fasta, simulations, clinvar_ese_hit_simulation_file, ese_file, window_start, window_end, exclude_cpg, clinvar=True)
        if not only_disease:
            print("Simulating ESE hits on the {0}-{1} region for 1000 genomes PTCs...".format(window_start, window_end))
            dso.ese_hit_simulation(unique_ptcs, snp_relative_positions_file, ptc_file, coding_exons_fasta, simulations, kgenomes_ese_hit_simulation_file, ese_file, window_start, window_end, exclude_cpg)


    excess_test_file = "{0}/clinvar_ptc_{1}_{2}_excesses.csv".format(output_directory, window_start, window_end)
    if excess_test:
        dso.excess_test(unique_ptcs_rel_pos_file, coding_exons_fasta, excess_test_file)

    location_test_file = "{0}/clinvar_locations_chisquare.csv".format(output_directory)
    if disease_locations_chisquare:
        dso.disease_ptc_location_test(unique_ptcs_rel_pos_file, coding_exons_fasta, location_test_file)


if __name__ == "__main__":
    main()
