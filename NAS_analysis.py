'''
Author: Rosina Savisaar and Liam Abrahams
Check whether PTCs are associated with greater rates of exon skipping.
'''

import bed_ops as bo
import bam_ops as bmo
import NAS_ops as nao
import generic as gen
import os
import random
import SNP_ops as so
import time
import numpy as np

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix", "bam_analysis_folder", "number_of_simulations", "simulation_output_folder", "motif_file", "filter_genome_data", "get_SNPs", "process_bams", "simulate_ptc_snps", "motif_complement", "overwrite_intersect", "use_old_sims", "out_of_frame", "simulate_ptcs_with_monomorphic", "generate_monomorphic_indices", "ignore_determine_snp_type", "ignore_psi_calculation", "ptc_location_analysis"], flags = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], ints = [7])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix, bam_analysis_folder, number_of_simulations, simulation_output_folder, motif_file, filter_genome_data, get_SNPs, process_bams, simulate_ptc_snps, motif_complement, overwrite_intersect, use_old_sims, out_of_frame, simulate_ptcs_with_monomorphic, generate_monomorphic_indices, ignore_determine_snp_type, ignore_psi_calculation, ptc_location_analysis = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix, args.bam_analysis_folder, args.number_of_simulations, args.simulation_output_folder, args.motif_file, args.filter_genome_data, args.get_SNPs, args.process_bams, args.simulate_ptc_snps, args.motif_complement, args.overwrite_intersect, args.use_old_sims, args.out_of_frame, args.simulate_ptcs_with_monomorphic, args.generate_monomorphic_indices, args.ignore_determine_snp_type, args.ignore_psi_calculation, args.ptc_location_analysis

    start = time.time()

    # create any necessary output diretories
    directory_splits = out_prefix.split('/')
    directory_paths = "/".join(directory_splits[:-1])
    gen.create_output_directories(directory_paths)
    gen.create_directory('temp_data/')

    CDS_fasta = "{0}_CDS.fasta".format(out_prefix)
    CDS_bed = "{0}_CDS.bed".format(out_prefix)
    exon_bed = "{0}_exons.bed".format(out_prefix)
    filtered_exon_bed = "{0}_filtered_exons.bed".format(out_prefix)
    exon_junctions_file = "{0}_exon_junctions.bed".format(out_prefix)
    coding_exon_bed = "{0}_coding_exons.bed".format(out_prefix)

    if filter_genome_data:
        #extract and filter CDS coordinates and sequences
        print("Extracting and filtering CDSs...")
        bo.extract_cds(gtf, CDS_bed, CDS_fasta, genome_fasta, all_checks = True, uniquify = True, clean_chrom_only = True, full_chr_name = True)
        gen.get_time(start)

        #group the CDS sequences into families based on sequence similarity
        print("Grouping sequences into families...")
        names = gen.read_fasta(CDS_fasta)[0]
        gen.find_families_ensembl("../source_data/GRCh37_ensembl_protein_families.txt", names, "{0}_families.txt".format(out_prefix))
        gen.get_time(start)

        print("Extracting and filtering exons...")
        #extract exon coordinates
        bo.extract_exons(gtf, exon_bed)
        #only leave exons from transcripts that passed quality control in the extract_cds step above.
        #also only leave a single gene per family
        bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed, families_file = "{0}_families.txt".format(out_prefix))
        gen.get_time(start)

        #extract exon-exon junction coordinates
        print("Extracting exon-exon junctions...")
        bo.extract_exon_junctions(exon_bed, exon_junctions_file, window_of_interest = 2)
        gen.get_time(start)

        #make another exons bed that only contains fully coding exons.
        #This is because in the final analysis, we should only consider fully protein-coding exons.
        #However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
        #be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
        print("Filtering out overlapping, non-coding and partially coding, as well as terminal exons...")
        bo.check_coding(filtered_exon_bed, CDS_bed, coding_exon_bed, remove_overlapping = True)
        gen.get_time(start)

    SNP_file = "{0}_SNP_file.txt".format(out_prefix)
    if out_of_frame:
        out_prefix = out_prefix + "_out_of_frame"
    PTC_file = "{0}_ptc_file.txt".format(out_prefix)
    syn_nonsyn_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)
    CDS_interval_file = "{0}_intervals{1}".format(os.path.splitext(CDS_fasta)[0], os.path.splitext(CDS_fasta)[1])
    #check which individuals were included in Geuvadis
    full_sample_names = os.listdir(bams_folder)
    full_sample_names = [i for i in full_sample_names if i[-4:] == ".bam" and "proc" not in i]
    sample_names = [(i.split("."))[0] for i in full_sample_names]
    sample_names = [i for i in sample_names if len(i) > 0]
    print('{0} samples included in Geuvadis...'.format(len(sample_names)))
    #for some reason, 17 of the samples from Geuvadis don't appear in the 1000genomes vcf
    #I'm gonna have to get to the bottom of this at some point
    #but at the moment I'm just gonna filter them out


    with open("../source_data/samples_in_vcf.txt") as file:
        samples_in_vcf = file.readlines()
    samples_in_vcf = [i.rstrip("\n") for i in samples_in_vcf]
    sample_names = [i for i in sample_names if i in samples_in_vcf]
    print('{0} samples also in vcf...'.format(len(sample_names)))
    sample_file = "{0}_sample_file.txt".format(out_prefix)


    # create a fasta containing all sequences for exons with snp
    coding_exons_fasta = "{0}_coding_exons.fasta".format(out_prefix)
    bo.fasta_from_intervals(coding_exon_bed, coding_exons_fasta, genome_fasta, names=True)

    if get_SNPs:
        #get SNPs for the sample
        intersect_file = "{0}_SNP_CDS_intersect.bed".format(out_prefix)
        print("Getting SNP data...")
        so.get_snps_in_cds(coding_exon_bed, CDS_bed, vcf_folder, panel_file, sample_names, sample_file, intersect_file, out_prefix)
        print("Calculating SNP positions...")
        so.get_snp_positions(sample_file, SNP_file, CDS_bed, intersect_file, out_prefix)
        gen.get_time(start)

    if ignore_determine_snp_type:
        pass
    else:
        print("Determining SNP type...")
        so.get_snp_change_status(SNP_file, CDS_fasta, PTC_file, syn_nonsyn_file, out_of_frame = out_of_frame, ref_check = True, headers=True)
        gen.get_time(start)

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    print("Filtering exon-exon junctions to only leave those that flank exons with a PTC variant...")
    PTC_exon_junctions_file = "{0}_filtered_exon_junctions.bed".format(out_prefix)
    bo.filter_exon_junctions(exon_junctions_file, PTC_file, PTC_exon_junctions_file)

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = ["{0}/{1}".format(bams_folder, i) for i in full_sample_names if (i.split("."))[0] in sample_names]

    #in parallel, do the processing on individual .bam files
    if bam_analysis_folder == "None":
        bam_analysis_folder = "{0}__analysis_bam_analysis".format(out_prefix)
    gen.create_directory(bam_analysis_folder)
    if process_bams:
        print("Processing RNA-seq data...")
        exon_junctions_bam_output_folder = "{0}__analysis_exon_junction_bams".format(out_prefix)
        if out_of_frame:
            splits = exon_junctions_bam_output_folder.split('/')
            splits[-1] = splits[-1].replace('_out_of_frame', '')
            exon_junctions_bam_output_folder = "/".join(splits)
        gen.create_directory(exon_junctions_bam_output_folder)
        #we have to do it like this because you can't pass flags into run_in_parallel
        keyword_dict = {"overwrite_intersect": overwrite_intersect}
        processes = gen.run_in_parallel(bam_files, ["foo", exon_junctions_file, PTC_exon_junctions_file, bam_analysis_folder, PTC_file, syn_nonsyn_file, out_prefix, exon_junctions_bam_output_folder, keyword_dict], nao.process_bam_per_individual, workers = 36)
        for process in processes:
            process.get()
        gen.get_time(start)

    #if required, filter PTCs to only leave ones that overlap motifs from a specified set
    motif_filtering = False
    if motif_file != "None":
        print("Filtering SNPs based on whether or not they overlap a motif from the specified set...")
        motif_suffix = ((motif_file.split("/"))[-1]).split(".")[0]
        if motif_complement:
            out_prefix = "{0}_{1}_complement".format(out_prefix, motif_suffix)
        else:
            out_prefix = "{0}_{1}".format(out_prefix, motif_suffix)
        filtered_ptc = "{0}_ptc_file.txt".format(out_prefix)
        so.filter_motif_SNPs(CDS_fasta, PTC_file, motif_file, filtered_ptc, complement = motif_complement)
        PTC_file = filtered_ptc

    final_file = "{0}__analysis_final_output.txt".format(out_prefix)
    if ignore_psi_calculation:
        pass
    else:
        print("Calculating PSI...")
        bmo.compare_PSI(PTC_file, bam_analysis_folder, final_file)

    #run the simulation that swaps ptcs for nonsynonymous snps
    if simulate_ptc_snps:
        if simulate_ptc_snps and not number_of_simulations:
            print("Please specify the number of simulations")
            raise Exception
        nao.ptc_snp_simulation(out_prefix, simulation_output_folder, PTC_file, syn_nonsyn_file, exon_junctions_file, bam_files, number_of_simulations, use_old_sims = use_old_sims)

    # run the simulation that picks monomorphic sites
    if simulate_ptcs_with_monomorphic:
        if simulate_ptcs_with_monomorphic and not number_of_simulations:
            print("Please specify the number of simulations")
            raise Exception

        coding_exon_fasta = "{0}_coding_exons.fasta".format(out_prefix)
        if not os.path.exists(coding_exon_fasta):
            print('Coding exon fasta is required...')
            raise Exception
        nao.ptc_monomorphic_simulation(out_prefix, simulation_output_folder, sample_file, genome_fasta, PTC_file, syn_nonsyn_file, coding_exon_bed, coding_exon_fasta, exon_junctions_file, bam_files, number_of_simulations, generate_indices = generate_monomorphic_indices, use_old_sims = use_old_sims)

    # get the locations of the ptcs
    if ptc_location_analysis:
        print("PTC locations analysis...")
        snp_relative_exon_position_file = "{0}_SNP_relative_exon_position.bed".format(out_prefix)
        ptc_location_analysis_output_file ="{0}_ptc_location_analysis.csv".format(out_prefix)
        coding_exon_fasta = "{0}_coding_exons.fasta".format(out_prefix)
        if not os.path.exists(coding_exon_fasta) or not os.path.exists(snp_relative_exon_position_file) or not os.path.exists(PTC_file):
            print("Please run --filter_genome_data and --get_SNPs first...")
            raise Exception
        # need to work out where and what the analysis outputs need to do
        so.ptc_locations(PTC_file, snp_relative_exon_position_file, ptc_location_analysis_output_file)


if __name__ == "__main__":
    main()
