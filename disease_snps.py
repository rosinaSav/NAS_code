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
    arguments = ["disease_snps_file", "output_directory", "results_prefix", "intersect_snps", "get_relative_positions", "get_snp_status", "get_info", "simulate_ptc_location", "get_possible_ptc_locations", "required_simulations", "get_overlaps", "intersect_ptcs", "compare_ptcs" ,"get_introns", "compare_distances", "clinvar_ptc_locations"]
    args = gen.parse_arguments(description, arguments, flags = [3,4,5,6,7,8,9,10,11,12,13,14,15])
    disease_snps_file, output_directory, results_prefix, intersect_snps, get_relative_positions, get_snp_status, get_info, simulate_ptc_location, get_possible_ptc_locations, required_simulations, get_overlaps, intersect_ptcs, compare_ptcs, get_introns, compare_distances, clinvar_ptc_locations = args.disease_snps_file, args.output_directory, args.results_prefix, args.intersect_snps, args.get_relative_positions, args.get_snp_status, args.get_info, args.simulate_ptc_location, args.get_possible_ptc_locations, args.required_simulations, args.get_overlaps, args.intersect_ptcs, args.compare_ptcs, args.get_introns, args.compare_distances, args.clinvar_ptc_locations

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
    if intersect_snps or not os.path.isfile(disease_snp_intersect_file_vcf) or not os.path.isfile(disease_snp_intersect_file_bed):
        print("Intersecting snps with exons")
        so.intersect_snps_parallel(exon_bed, disease_snps_file, disease_snp_intersect_file_vcf)
        so.intersect_vcf_to_bed(exon_bed, disease_snp_intersect_file_vcf, disease_snp_intersect_file_bed, change_names = True)

    # get relative positions of the snps in cds and exons
    full_bed = "{0}_CDS.bed".format(results_prefix)
    disease_snps_relative_exon_positions = "{0}/disease_snp_relative_exon_positions.bed".format(output_directory)
    disease_snps_relative_cds_positions = "{0}/disease_snp_relative_cds_positions.bed".format(output_directory)
    if get_relative_positions or not os.path.isfile(disease_snps_relative_exon_positions) or not os.path.isfile(disease_snps_relative_cds_positions):
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(disease_snp_intersect_file_bed, disease_snps_relative_exon_positions)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, disease_snps_relative_cds_positions, full_bed)

    # get the change status of the snps to check them
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    disease_ptcs_file = "{0}/disease_ptcs.txt".format(output_directory)
    disease_other_file = "{0}/disease_other_snps.txt".format(output_directory)
    if get_snp_status or not os.path.isfile(disease_ptcs_file) or not os.path.isfile(disease_other_file) or not os.path.isfile(cds_fasta):
        print("Getting snp status...")
        so.get_snp_change_status(disease_snps_relative_cds_positions, cds_fasta, disease_ptcs_file, disease_other_file)

    # get the info about the ptcs
    output_file_ptc_info = "{0}/disease__analysis_ptc_info.txt".format(output_directory)
    output_file_other_info = "{0}/disease__analysis_other_info.txt".format(output_directory)
    if get_info:
        print("Getting PTC information...")
        dso.get_ptc_info(disease_ptcs_file, disease_snps_relative_exon_positions, output_file_ptc_info)
        dso.get_ptc_info(disease_other_file, disease_snps_relative_exon_positions, output_file_other_info)

    location_simulation_output_directory = "{0}/possible_ptc_locations".format(output_directory)
    coding_exons_file = "{0}_exons.bed".format(results_prefix)
    if get_possible_ptc_locations:
        print("Getting possible PTC mutation locations...")
        dso.generate_possible_ptc_locations(full_bed, cds_fasta, location_simulation_output_directory)

    # simulation to see if disease ptcs occur at exon ends more conmonly than by chance
    if simulate_ptc_location:
        if not required_simulations:
            print("\nERROR: please specify the number of simulations required.\n")
        dso.ptc_location_simulation(disease_ptcs_file, full_bed, cds_fasta, location_simulation_output_directory, required_simulations, coding_exons_file)

    ptc_file = "{0}_ptc_file.txt".format(results_prefix)
    overlap_file = "{0}/disease__analysis_overlaps.txt".format(output_directory)
    if get_overlaps:
        print("Getting overlap between disease PTCs and 1000 genomes PTCs...")
        dso.get_ptc_overlaps(disease_ptcs_file, ptc_file, overlap_file)

    # get intersect between the clinvar ptcs and 1000 genomes ptcs
    ptc_intersect_file = "{0}/ptc_intersect.bed".format(output_directory)
    if intersect_ptcs:
        temp_disease_ptc_file = "temp_data/{0}".format(random.random())
        dso.refactor_ptc_file(disease_ptcs_file, temp_disease_ptc_file)
        temp_k_genomes_ptc_file = "temp_data/{0}".format(random.random())
        dso.refactor_ptc_file(ptc_file, temp_k_genomes_ptc_file, header=True)
        bao.intersect_bed(temp_k_genomes_ptc_file, temp_disease_ptc_file, write_both = True, no_dups=False, output_file = ptc_intersect_file)
        gen.remove_file(temp_disease_ptc_file)
        gen.remove_file(temp_k_genomes_ptc_file)

    # get the introns
    cds_bed_file = "{0}_CDS.bed".format(results_prefix)
    intron_bed = "{0}/internal_intron_list.bed".format(output_directory)
    # intron_fasta = "{0}/internal_intron_list.fasta".format(output_directory)
    if get_introns:
        print("Extracting internal introns from exons...")
        dso.get_introns(cds_bed_file, intron_bed)
        # beo.fasta_from_intervals(intron_bed, intron_fasta, genome_file, force_strand = True, names = True)

    compare_file = "{0}/ptc_distances_disease.txt".format(output_directory)
    relative_exon_positions_file = "{0}_SNP_relative_exon_position.bed".format(results_prefix)
    exon_fasta = "{0}_CDS_intervals.fasta".format(results_prefix)
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    if compare_ptcs:
        dso.compare_ptcs(ptc_intersect_file, ptc_file, relative_exon_positions_file, exon_fasta, cds_fasta, cds_bed_file, intron_bed, compare_file)

    kgenomes_relative_exon_positions_file = "{0}_SNP_relative_exon_position.bed".format(results_prefix)
    compare_distances_file = "{0}/compare_distances.txt".format(output_directory)
    if compare_distances:
        dso.compare_distances(disease_ptcs_file, disease_snps_relative_exon_positions, ptc_file, kgenomes_relative_exon_positions_file, compare_distances_file)

    clinvar_ptc_locations_file = "{0}/clinvar_ptc_locations.csv".format(output_directory)
    if clinvar_ptc_locations:
        dso.clinvar_ptc_locations(disease_ptcs_file, disease_snps_relative_exon_positions, ptc_file, clinvar_ptc_locations_file)

    # coding_exons_fasta = "{0}_coding_exons.fasta".format(results_prefix)
    # clinvar_location_sim = "{0}/clinvar_simulations.csv".format(output_directory)
    # if clinvar_simulation:
    #     dso.clinvar_simulation(disease_ptcs_file, disease_snps_relative_exon_positions, ptc_file, coding_exons_fasta, )


if __name__ == "__main__":
    main()
