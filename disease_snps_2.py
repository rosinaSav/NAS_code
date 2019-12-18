'''
Auhor: Liam Abrahams.
Look at disease SNPs from ClinVar.
'''

import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import disease_snps_ops2 as dso2
import os
import collections
import re
import copy
import numpy as np
import random

def main():

    description = "Look at disease snps."
    arguments = ["region_chi", "region_sim"]
    args = gen.parse_arguments(description, arguments, flags = [0,1])
    region_chi, region_sim = args.region_chi, args.region_sim

    output_directory = "results/clinvar3"
    gen.create_output_directories(output_directory)

    coding_exons_file = "results/clean_run_2/clean_run_coding_exons.fasta"
    # ptc_file = "results/clean_run_2/.bed"
    clinvar_ptc_file = "results/clinvar2/disease_ptcs_no_intersect_rel_pos_pathogenicity.bed"


    if region_chi:
        output_file = "{0}/locations_chi.csv".format(output_directory)
        file1 = "results/clinvar3/disease_ptcs_no_intersect_rel_pos_pathogenic.bed".format(output_directory)
        file2 = "results/clinvar3/disease_ptcs_no_intersect_rel_pos_likely_pathogenic.bed".format(output_directory)
        dso2.ptc_locations_test(coding_exons_file, file1, file2, output_file)

    if region_sim:
        output_file = "{0}/locations_sim.csv".format(output_directory)
        dso2.ptc_location_sim(coding_exons_file, clinvar_ptc_file, output_file)

    # # simulation picking random reference allele matched simulants
    # clinvar_location_simulation_file = "{0}/clinvar_ptc_location_simulation.csv".format(output_directory)
    # clinvar_location_simulation_ese_overlap_file = "{0}/clinvar_ptc_location_simulation_{1}_ese_overlaps.csv".format(output_directory, ese_file_name)
    # kgenomes_location_simulation_file = "{0}/1000_genomes_simulations.csv".format(output_directory)
    # kgenomes_location_simulation_ese_overlap_file = "{0}/1000_genomes_simulations_ese_overlaps.csv".format(output_directory)
    #
    # if location_simulation:
    #     if not only_kgenomes:
    #         print('Running ptc location simulation on disease PTCs...')
    #         dso.ptc_location_simulation(unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, clinvar_location_simulation_file, clinvar_location_simulation_ese_overlap_file, ese_file, only_ese, exclude_cpg)
    #     if not only_disease:
    #         print('Running ptc location simulation on 1000 genomes PTCs...')
    #         dso.ptc_location_simulation(kgenomes_unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, kgenomes_location_simulation_file, kgenomes_location_simulation_ese_overlap_file, ese_file, only_ese, exclude_cpg)
    #
    #
    # window_start = 3
    # window_end = 69
    # clinvar_ese_hit_simulation_file = "{0}/clinvar_ese_hit_simulation_{1}_{2}_{3}.csv".format(output_directory, window_start, window_end, ese_file_name)
    # kgenomes_ese_hit_simulation_file = "{0}/1000_genomes_ese_hit_simulation_{1}_{2}_{3}.csv".format(output_directory, window_start, window_end, ese_file_name)
    #
    # # do a simulation picking only sites from within the region
    # if ese_hit_simulation:
    #     if not only_kgenomes:
    #         print("Simulating ESE hits on the {0}-{1} region for disease PTCs...".format(window_start, window_end))
    #         dso.ese_hit_simulation(unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, clinvar_ese_hit_simulation_file, ese_file, window_start, window_end, exclude_cpg)
    #     if not only_disease:
    #         print("Simulating ESE hits on the {0}-{1} region for 1000 genomes PTCs...".format(window_start, window_end))
    #         dso.ese_hit_simulation(kgenomes_unique_ptcs_rel_pos_file, coding_exons_fasta, simulations, kgenomes_ese_hit_simulation_file, ese_file, window_start, window_end, exclude_cpg)
    #
    #
    # excess_test_file = "{0}/clinvar_ptc_{1}_{2}_excesses.csv".format(output_directory, window_start, window_end)
    # if excess_test:
    #     dso.excess_test(unique_ptcs_rel_pos_file, coding_exons_fasta, excess_test_file)




if __name__ == "__main__":
    main()
