'''
Author: Rosina Savisaar and Liam Abrahams
Operations use in NAS_analysis.
'''

import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import random
import SNP_ops as so
import time
import numpy as np


def get_non_mutation_indices(simulation_output_folder, sample_file, coding_exon_bed, out_prefix, genome_fasta, nt_indices_files):
    '''
    Get the indices of all the positions where no SNPs/other variants are reported.
    '''

    # to do:
    # 1. need to convert bed to not have chr name
    # for each vcf file:
    # 1. subtract bed to get all locations where there isnt a mutation using intersect_bed with subtract=True
    # bmo.intersect_bed(file1, file2, output_file, subtract=True)
    # 2. using that, extract all sequence pieces using fasta_from_intervals, need to decide if names=true or names=false
    # fasta_from_intervals(intersect_bed_file, fasta_file, genome_fasta, force_strand = True, names = False)
    # 3. from this fasta file, create a file that contains for each exon, the indices at which each nt resides
    # split this so each chr has its own file?
    # bo.extract_nt_indices(fasta_file, output_file)

    # set up the new file to contain the regions without a mutation
    coding_exon_bed_out = "{0}/format_{1}".format(simulation_output_folder, coding_exon_bed.split('/')[-1])
    # set up the new file to contain the regions without a mutation
    coding_exon_bed_out_subtract = "{0}/subtract_{1}".format(simulation_output_folder, coding_exon_bed.split('/')[-1])
    # set up the new file to contain the regions without a mutation
    coding_exon_bed_out_subtract_format = "{0}/subtract_format_{1}".format(simulation_output_folder, coding_exon_bed.split('/')[-1])
    # file to contain the fasta output of regions containing no mutation
    fasta_no_mutations = "{0}/exon_regions_no_mutations.fasta".format(simulation_output_folder)
    # change the names in the bed file to correspond to the mutation vcf
    bo.change_bed_names(coding_exon_bed, coding_exon_bed_out, full_names=True, header=False)
    # intersect the bed, leaving only regions that contain no mutations
    bmo.intersect_bed(coding_exon_bed, sample_file, write_both = False, output_file = coding_exon_bed_out_subtract, no_dups = False, subtract=True, intersect=True)
    # file = gen.read_many_fields(coding_exon_bed_out_subtract, "\t")
    bo.change_bed_names(coding_exon_bed_out_subtract, coding_exon_bed_out_subtract_format, full_names=False, header=False)
    # generate fasta of all the regions
    bo.fasta_from_intervals(coding_exon_bed_out_subtract, fasta_no_mutations, genome_fasta, force_strand = True, names = True)
    # # extract the indices of each location a mutation doesnt occur
    bo.extract_nt_indices(fasta_no_mutations, nt_indices_files)

def run_ptc_monomorphic_simulation_instance(simulations, out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, syn_nonsyn_file, exon_junctions_file, bam_files, nt_indices_files, coding_exon_fasta, parallel = False, use_old_sims = False):
    '''
    Run the ptc simulations for the required number.
    '''

    #iterate over simulations
    counter = 0
    for simulation_number in simulations:

        counter = gen.update_counter(counter, 10, "SIMULATION ")

        #setup a folder to contain the individual simulation inside the simulations output
        simulation_instance_folder = "{0}/ptc_monomorphic_simulation_run_{1}".format(simulation_output_folder, simulation_number)
        if not use_old_sims:
            gen.create_strict_directory(simulation_instance_folder)
        else:
            gen.create_directory(simulation_instance_folder)

        # copy ptc file to directory
        real_ptcs_for_sim_file = "{0}/{1}".format(simulation_output_folder, ptc_file.split('/')[-1])
        gen.copy_file(ptc_file, real_ptcs_for_sim_file)
        ptc_file = real_ptcs_for_sim_file

        #get list of exons
        exon_list = bo.get_fasta_exon_intervals(coding_exon_fasta)

        #generate pseudo ptc snps
        pseudo_monomorphic_ptc_file = "{0}/pseudo_monomorphic_ptc_file_{1}.txt".format(simulation_instance_folder, simulation_number)
        if (not use_old_sims) or (not(os.path.isfile(pseudo_monomorphic_ptc_file))):
            so.generate_pseudo_monomorphic_ptcs(ptc_file, nt_indices_files, exon_list, pseudo_monomorphic_ptc_file)

        #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step when generating pseudo ptcs
        pseudo_monomorphic_ptc_exon_junctions_file = "{0}/filtered_exon_junctions_{1}.bed".format(simulation_instance_folder, simulation_number)
        if (not use_old_sims) or (not(os.path.isfile(pseudo_monomorphic_ptc_exon_junctions_file))):
            bo.filter_exon_junctions(exon_junctions_file, pseudo_monomorphic_ptc_file, pseudo_monomorphic_ptc_exon_junctions_file)

        exon_junctions_bam_output_folder = "{0}__analysis_exon_junction_bams".format(out_prefix)
        gen.create_directory(exon_junctions_bam_output_folder)
        #run the bam analysis for each
        #(don't parallelize if you're doing the simulations in parallel)
        kw_dict = {"ptc_snp_simulation": True, "simulation_instance_folder": simulation_instance_folder, "simulation_number": simulation_number}
        if parallel:
            process_bam_per_individual(bam_files, exon_junctions_file, pseudo_monomorphic_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_monomorphic_ptc_file, syn_nonsyn_file, out_prefix, exon_junctions_bam_output_folder, kw_dict)
        else:
            processes = gen.run_in_parallel(bam_files, ["foo", exon_junctions_file, pseudo_monomorphic_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_monomorphic_ptc_file, syn_nonsyn_file, out_prefix, exon_junctions_bam_output_folder, kw_dict], process_bam_per_individual, workers = 36)
            for process in processes:
                process.get()

        #process final psi for simulation
        final_file = "{0}/final_output_simulation_{1}.txt".format(simulation_bam_analysis_output_folder, simulation_number)
        bmo.compare_PSI(pseudo_monomorphic_ptc_file, simulation_bam_analysis_output_folder, final_file, sim_number = simulation_number)

def ptc_monomorphic_simulation(out_prefix, simulation_output_folder, sample_file, genome_fasta, ptc_file, syn_nonsyn_file, coding_exon_bed, coding_exon_fasta, exon_junctions_file, bam_files, required_simulations, generate_indices = False, use_old_sims = False):
    '''
    Set up the PTC simulations and then run.
    if use_old_sims is True, don't pick new simulant SNPs from monomorphic sites.
    '''

    print("Running simulation picking monomorphic sites that have the same ancestral allele as a PTC snp...")

    #setup up simulation output folder
    if simulation_output_folder == "None":
        simulation_output_folder = "{0}_simulate_ptc_monomorphic_sites".format(out_prefix)
    if not use_old_sims and generate_indices:
        #if the simulation folder we are specifying already exists, delete and start again
        gen.create_strict_directory(simulation_output_folder)
    else:
        gen.create_directory(simulation_output_folder)

    #setup up simulation bam analysis output folder
    simulation_bam_analysis_output_folder = "{0}_simulate_ptc_monomorphic_sites_bam_analysis".format(out_prefix)

    # create the filepaths to hold to positions of non mutated sites
    nt_indices_files = {
        "A": "{0}/nt_indices_no_mutations_A.fasta".format(simulation_output_folder),
        "C": "{0}/nt_indices_no_mutations_C.fasta".format(simulation_output_folder),
        "G": "{0}/nt_indices_no_mutations_G.fasta".format(simulation_output_folder),
        "T": "{0}/nt_indices_no_mutations_T.fasta".format(simulation_output_folder),
    }

    if generate_indices and not use_old_sims:
        get_non_mutation_indices(simulation_output_folder, sample_file, coding_exon_bed, out_prefix, genome_fasta, nt_indices_files)

    if not use_old_sims:
        #if the simulation folder we are specifying already exists, delete and start again
        gen.create_strict_directory(simulation_bam_analysis_output_folder)
    else:
        gen.create_directory(simulation_bam_analysis_output_folder)

    # #create a list of simulations to iterate over
    simulations = list(range(1, required_simulations+1))
    # #if you're only doing one simulation, don't parallelize the simulations
    # #parallelize the processing of bams like for true data
    if required_simulations > 1:
        processes = gen.run_in_parallel(simulations, ["foo", out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, syn_nonsyn_file, exon_junctions_file, bam_files, nt_indices_files, coding_exon_fasta, True, use_old_sims], run_ptc_monomorphic_simulation_instance, workers = 36)
        for process in processes:
            process.get()
    else:
        run_ptc_monomorphic_simulation_instance([1], out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, syn_nonsyn_file, exon_junctions_file, bam_files, nt_indices_files, coding_exon_fasta, False, use_old_sims)

def ptc_snp_simulation(out_prefix, simulation_output_folder, ptc_file, syn_nonsyn_file, exon_junctions_file, bam_files, required_simulations, exon_junctions_bam_output_folder, use_old_sims = False):
    '''
    Set up the PTC simulations and then run.
    if use_old_sims is True, don't pick new simulant SNPs.
    '''

    #setup up simulation output folder
    if simulation_output_folder == "None":
        simulation_output_folder = "{0}_simulate_ptc_snps".format(out_prefix)
    if not use_old_sims:
        #if the simulation folder we are specifying already exists, delete and start again
        gen.create_strict_directory(simulation_output_folder)
    else:
        gen.create_directory(simulation_output_folder)



    #setup up simulation bam analysis output folder
    simulation_bam_analysis_output_folder = "{0}__analysis_simulation_ptc_snps_bam_analysis".format(out_prefix)
    if not use_old_sims:
        #if the simulation folder we are specifying already exists, delete and start again
        gen.create_strict_directory(simulation_bam_analysis_output_folder)
    else:
        gen.create_directory(simulation_bam_analysis_output_folder)

    #get all nonsynonymous snps and put them in the simulation output folder
    nonsynonymous_snps_file = "{0}/nonsynonymous_snps.txt".format(simulation_output_folder)
    so.filter_by_snp_type(syn_nonsyn_file, nonsynonymous_snps_file, "non")

    #create a list of simulations to iterate over
    simulations = list(range(1, required_simulations+1))
    #if you're only doing one simulation, don't parallelize the simulations
    #parallelize the processing of bams like for true data
    if required_simulations > 1:
        processes = gen.run_in_parallel(simulations, ["foo", out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, exon_junctions_bam_output_folder, True, use_old_sims], run_ptc_simulation_instance)
        for process in processes:
            process.get()
    else:
        run_ptc_simulation_instance([1], out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, exon_junctions_bam_output_folder, False, use_old_sims)

def process_bam_per_individual(bam_files, global_exon_junctions_file, PTC_exon_junctions_file, out_folder, PTC_file, syn_nonsyn_file, out_prefix, exon_junctions_bam_output_folder, kw_dict):
    '''
    Do all of the processing on an individual bam, from filtering out low quality data to mapping reads to
    exon-exon junctions.
    For each exon, return information on how many reads fall at different exon-exon junctions.
    '''

    #parse keyword_dict
    #it's done like this to make it easier to parallelize this process
    if "ptc_snp_simulation" in kw_dict:
        ptc_snp_simulation = kw_dict["ptc_snp_simulation"]
    else:
        ptc_snp_simulation = False
    if "simulation_instance_folder" in kw_dict:
        simulation_instance_folder = kw_dict["simulation_instance_folder"]
    else:
        simulation_instance_folder = None
    if "simulation_number" in kw_dict:
        simulation_number = kw_dict["simulation_number"]
    else:
        simulation_number = None
    if "overwrite_intersect" in kw_dict:
        overwrite_intersect = kw_dict["overwrite_intersect"]
    else:
        overwrite_intersect = False
    if "phase" in kw_dict:
        phase = kw_dict["phase"]
    else:
        phase = False

    bam_file_number = len(bam_files)
    for pos, bam_file in enumerate(bam_files):

        #Process:
        # 1. get the number of reads in bam
        # 2. Filter out reads that don't overlap exon-exon junctions
        # 3. Filter out reads that don't overlap exon-exon junctions flanking PTC-containing exons
        # 4. Filter bams by quality
        # This gives us a set of "good" quality reads.
        # 5. scale down total read number proportionally to how many reads were lost in the quality filtering
        # 6. Count reads either skipping or including each exon

        print("{0}/{1}: {2}".format(pos, bam_file_number, bam_file))
        sample_name = (bam_file.split("/")[-1]).split(".")[0]
        if ptc_snp_simulation:
            output_file = "{0}/{1}_simulation_{2}.txt".format(out_folder, sample_name, simulation_number)
        else:
            output_file = "{0}/{1}.txt".format(out_folder, sample_name)

        #folder that will contain all of the intermediate steps in the processing of the bam file
        if ptc_snp_simulation:
            proc_folder = "{0}/bam_proc_files".format(simulation_instance_folder)

        else:
            proc_folder = "{0}__analysis_bam_proc_files".format(out_prefix)

        gen.create_output_directories(proc_folder)

        bam_file_parts = os.path.split(bam_file)
        mapq_filtered_bam = "{0}/{1}_filtered_mapq.bam".format(proc_folder, bam_file_parts[1])
        mapq_flag_filtered_bam = "{0}_flag.bam".format(mapq_filtered_bam[:-4])
        mapq_flag_xt_filtered_bam = "{0}_xt.bam".format(mapq_flag_filtered_bam[:-4])
        mapq_flag_xt_nm_filtered_bam = "{0}_nm.bam".format(mapq_flag_xt_filtered_bam[:-4])


        if not os.path.isfile(output_file):


            #1: We get a count of the total reads in the sample which can be used for normalisation
            #I'm initializing it to None for safety. That way, if the process fails,
            #it won't just silently go with whatever the value was at the end of the previous loop.
            #also, writing it down cause this bit takes forever, don't want to do it again every time.
            read_count_file_name = "{0}/read_count_sample_name.txt".format(exon_junctions_bam_output_folder)
            read_count = None
            if os.path.isfile(read_count_file_name):
                with open(read_count_file_name) as file:
                    read_count = int("".join(file))
            else:
                read_count = int(gen.run_process(["samtools", "view", "-c", bam_file]))
                with open(read_count_file_name, "w") as file:
                    file.write(str(read_count))

            #2: intersect the bam with all exon-exon junctions
            #only has to be done once for each bam
            #also removing "_out_of_frame" from out_prefix if it is present
            global_out_prefix = out_prefix
            if "out_of_frame" in global_out_prefix:
                global_out_prefix = global_out_prefix[:6]
            global_intersect_bam = "{0}/{1}_exon_junctions.bam".format(exon_junctions_bam_output_folder, bam_file_parts[1][:-4])
            if not os.path.isfile(global_intersect_bam) or overwrite_intersect:
                #intersect the filtered bam and the global exon junctions file
                # print(global_intersect_bam)
                bmo.intersect_bed(bam_file, global_exon_junctions_file, output_file = global_intersect_bam, intersect_bam = True)

            #3: filter to relevant exon-exon junctions
            ##Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
            #this uses intersect bed, with the intersect bam parameter
            intersect_bam = "{0}/{1}_exon_junction_bam_intersect.bam".format(proc_folder, bam_file_parts[1][:-4])

            #intersect the filtered bam and the ptc exon junctions file
            bmo.intersect_bed(global_intersect_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)

            #count how many reads there are in the sample after filtering to relevant exon-exon junctions but before quality filtering
            read_count_junctions_no_filter = int(gen.run_process(["samtools", "view", "-c", intersect_bam]))
            #4. filter .bam alignments by quality.
            #takes both upper and lower bam thresholds
            #outputs bam file with "_quality_filter_{lower_lim}_{upper_lim}" appended
            # need to do this twice and merge, so we use both intervals used by Geuvadis
            #set the mapq filter parameters here
            mapq_intervals = [[251, 255], [175, 181]]
            mapq_filter_filelist = []

            for mapq_interval in mapq_intervals:
                lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
                mapq_filter_file = "{0}/{1}_mapq_filter_{2}_{3}.bam".format(proc_folder, bam_file_parts[1][:-4], lower_threshold, upper_threshold)
                mapq_filter_filelist.append(mapq_filter_file)
                ##run the mapq filter
                bmo.bam_quality_filter(intersect_bam, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)

            ##merge files in filelist
            bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)

            ##filter by flags: get all mapped reads
            #Leaves: mapped unpaired and paired reads
            bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)

            ##filter bam by xt tag XT=U
            bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")

            ##filter bam by nm tag NM<=6
            bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)

            #5. scale down the initial count of reads in the sample by the proportion lost during quality filtering
            read_count_junctions_filter = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))
            prop_kept = np.divide(read_count_junctions_filter, read_count_junctions_no_filter)
            read_count = prop_kept * read_count

            #convert to sam format and phase reads
            intersect_sam = "{0}_phased.sam".format(mapq_flag_xt_nm_filtered_bam[:-4])
            if phase:
                temp_snp_file = "temp_data/snps{0}.txt".format(random.random())
                so.merge_and_header(PTC_file, syn_nonsyn_file, temp_snp_file)
                bmo.phase_bams(temp_snp_file, mapq_flag_xt_nm_filtered_bam, sample_name, intersect_sam)
                gen.remove_file(temp_snp_file)
            else:
                gen.run_process(["samtools", "view", mapq_flag_xt_nm_filtered_bam], file_for_output = intersect_sam)

            #6. count the number of reads supporting either the skipping or the inclusion of each exon
            junctions = bmo.read_exon_junctions(PTC_exon_junctions_file)
            bmo.count_junction_reads(intersect_sam, junctions, output_file, read_count)



def run_ptc_simulation_instance(simulations, out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, exon_junctions_bam_output_folder, parallel = False, use_old_sims = False):
    '''
    Run the ptc simulations for the required number.
    '''

    #iterate over simulations
    counter = 0
    for simulation_number in simulations:

        counter = gen.update_counter(counter, 10, "SIMULATION ")

        #setup a folder to contain the individual simulation inside the simulations output
        simulation_instance_folder = "{0}/ptc_simulation_run_{1}".format(simulation_output_folder, simulation_number)
        if not use_old_sims:
            gen.create_strict_directory(simulation_instance_folder)
        else:
            gen.create_directory(simulation_instance_folder)

        #generate pseudo ptc snps
        #also need to remove these snps from the file they started in so create a new remaining snps file
        #we can tweak these if we start running out of snps
        pseudo_ptc_file = "{0}/pseudo_ptc_file_{1}.txt".format(simulation_instance_folder, simulation_number)
        remaining_snps_file = "{0}/remaining_snps_file_{1}.txt".format(simulation_instance_folder, simulation_number)
        if (not use_old_sims) or (not(os.path.isfile(pseudo_ptc_file))):
            so.generate_pseudo_ptc_snps(ptc_file, nonsynonymous_snps_file, pseudo_ptc_file, remaining_snps_file, group_by_gene=False, without_replacement=True, match_allele_frequency=True, match_allele_frequency_window=0.05)

        #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step when generating pseudo ptcs
        pseudo_ptc_exon_junctions_file = "{0}/filtered_exon_junctions_{1}.bed".format(simulation_instance_folder, simulation_number)
        if (not use_old_sims) or (not(os.path.isfile(pseudo_ptc_file))):
            bo.filter_exon_junctions(exon_junctions_file, pseudo_ptc_file, pseudo_ptc_exon_junctions_file)

        # # print('parallel', parallel)
        # #run the bam analysis for each
        # #(don't parallelize if you're doing the simulations in parallel)
        # kw_dict = {"ptc_snp_simulation": True, "simulation_instance_folder": simulation_instance_folder, "simulation_number": simulation_number}
        # if parallel:
        #     process_bam_per_individual(bam_files, exon_junctions_file, pseudo_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_ptc_file, remaining_snps_file, out_prefix, exon_junctions_bam_output_folder, kw_dict)
        # else:
        #     # (bam_files, global_exon_junctions_file, PTC_exon_junctions_file, out_folder, PTC_file, syn_nonsyn_file, out_prefix, exon_junctions_bam_output_folder, kw_dict)
        #     processes = gen.run_in_parallel(bam_files, ["foo", exon_junctions_file, pseudo_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_ptc_file, remaining_snps_file, out_prefix, exon_junctions_bam_output_folder, kw_dict], process_bam_per_individual)
        #     for process in processes:
        #         process.get()
        #
        # #process final psi for simulation
        # final_file = "{0}/final_output_simulation_{1}.txt".format(simulation_bam_analysis_output_folder, simulation_number)
        # bmo.compare_PSI(pseudo_ptc_file, simulation_bam_analysis_output_folder, final_file, sim_number = simulation_number)
