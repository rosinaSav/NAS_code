import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import random
import SNP_ops as so
import time


def run_ptc_simulation_instance(simulations, out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, parallel = False):
    '''
    Run the ptc simulations for the required number.
    '''

    #iterate over simulations
    counter = 0
    for simulation_number in simulations:

        counter = gen.update_counter(counter, 10, "SIMULATION ")

        #setup a folder to contain the individual simulation inside the simulations output
        simulation_instance_folder = "{0}/ptc_simulation_run_{1}".format(simulation_output_folder, simulation_number)
        gen.create_strict_directory(simulation_instance_folder)

        #generate pseudo ptc snps
        #also need to remove these snps from the file they started in so create a new remaining snps file
        #we can tweak these if we start running out of snps
        pseudo_ptc_file = "{0}/pseudo_ptc_file_{1}.txt".format(simulation_instance_folder, simulation_number)
        remaining_snps_file = "{0}/remaining_snps_file_{1}.txt".format(simulation_instance_folder, simulation_number)
        so.generate_pseudo_ptc_snps(ptc_file, nonsynonymous_snps_file, pseudo_ptc_file, remaining_snps_file, group_by_gene=False, without_replacement=True, match_allele_frequency=True, match_allele_frequency_window=0.05)

        #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step when generating pseudo ptcs
        pseudo_ptc_exon_junctions_file = "{0}/filtered_exon_junctions_{1}.bed".format(simulation_instance_folder, simulation_number)
        bo.filter_exon_junctions(exon_junctions_file, pseudo_ptc_file, pseudo_ptc_exon_junctions_file)

        #run the bam analysis for each
        #(don't parallelize if you're doing the simulations in parallel)
        kw_dict = {"ptc_snp_simulation": True, "simulation_instance_folder": simulation_instance_folder, "simulation_number": simulation_number}
        if parallel:
            process_bam_per_individual(bam_files, exon_junctions_file, pseudo_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_ptc_file, remaining_snps_file, out_prefix, kw_dict)
        else:
            processes = gen.run_in_parallel(bam_files, ["foo", exon_junctions_file, pseudo_ptc_exon_junctions_file, simulation_bam_analysis_output_folder, pseudo_ptc_file, remaining_snps_file, out_prefix, kw_dict], process_bam_per_individual)
            for process in processes:
                process.get()

        #process final psi for simulation
        final_file = "{0}/final_output_simulation_{1}.txt".format(simulation_bam_analysis_output_folder, simulation_number)
        bmo.compare_PSI(pseudo_ptc_file, simulation_bam_analysis_output_folder, final_file)


def ptc_snp_simulation(out_prefix, simulation_output_folder, ptc_file, syn_nonsyn_file, exon_junctions_file, bam_files, required_simulations):
    '''
    Set up the PTC simulations and then run.
    '''

    #setup up simulation output folder
    if simulation_output_folder == "None":
        simulation_output_folder = "{0}_simulate_ptc_snps".format(out_prefix)
    #if the simulation folder we are specifying already exists, delete and start again
    gen.create_strict_directory(simulation_output_folder)

    #setup up simulation bam analysis output folder
    simulation_bam_analysis_output_folder = "{0}_simulate_ptc_snps_bam_analysis".format(out_prefix)
    #if the simulation folder we are specifying already exists, delete and start again
    gen.create_strict_directory(simulation_bam_analysis_output_folder)

    #get all nonsynonymous snps and put them in the simulation output folder
    nonsynonymous_snps_file = "{0}/nonsynonymous_snps.txt".format(simulation_output_folder)
    so.filter_by_snp_type(syn_nonsyn_file, nonsynonymous_snps_file, "non")

    #create a list of simulations to iterate over
    simulations = list(range(1, required_simulations+1))
    #if you're only doing one simulation, don't parallelize the simulations
    #parallelize the processing of bams like for true data
    if required_simulations > 1:
        processes = gen.run_in_parallel(simulations, ["foo", out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, True], run_ptc_simulation_instance)
        for process in processes:
            process.get()
    else:
        run_ptc_simulation_instance([1], out_prefix, simulation_output_folder, simulation_bam_analysis_output_folder, ptc_file, nonsynonymous_snps_file, exon_junctions_file, bam_files, False)

def process_bam_per_individual(bam_files, global_exon_junctions_file, PTC_exon_junctions_file, out_folder, PTC_file, syn_nonsyn_file, out_prefix, kw_dict):
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

    for bam_file in bam_files:

        #Process:
        # 1. get the number of reads in bam
        # 2. Filter out reads that don't overlap exon-exon junctions
        # 3. Filter out reads that don't overlap exon-exon junctions flanking PTC-containing exons
        # 4. Filter bams by quality
        # This gives us a set of "good" quality reads.
        # 5. scale down total read number proportionally to how many reads were lost in the quality filtering 
        # 6. Count reads either skipping or including each exon

        print(bam_file)
        sample_name = (bam_file.split("/")[-1]).split(".")[0]
        if ptc_snp_simulation:
            output_file = "{0}/{1}_simulation_{2}.txt".format(out_folder, sample_name, simulation_number)
        else:
            output_file = "{0}/{1}.txt".format(out_folder, sample_name)

        if ptc_snp_simulation:
            proc_folder = "{0}/bam_proc_files".format(simulation_instance_folder)

        else:
            proc_folder = "{0}_bam_proc_files".format(out_prefix)
        gen.create_directory(proc_folder)
        
        bam_file_parts = os.path.split(bam_file)
        mapq_filtered_bam = "{0}/{1}_filtered_mapq.bam".format(proc_folder, bam_file_parts[1])
        mapq_flag_filtered_bam = "{0}_flag.bam".format(mapq_filtered_bam[:-4])
        mapq_flag_xt_filtered_bam = "{0}_xt.bam".format(mapq_flag_filtered_bam[:-4])
        mapq_flag_xt_nm_filtered_bam = "{0}_nm.bam".format(mapq_flag_xt_filtered_bam[:-4])

        #1: We get a count of the total reads in the sample which can be used for normalisation
        read_count = int(gen.run_process(["samtools", "view", "-c", bam_file]))

        #2: intersect the bam with all exon-exon junctions
        #only has to be done once for each bam
        global_intersect_bam = "{0}_{1}_exon_junctions.bam".format(out_prefix, bam_file_parts[1][:-4])
        if (not os.path.isfile(global_intersect_bam)) or overwrite_intersect:
            print("Filtering .bam to all exon-exon junctions...")
            start = time.time()
            #intersect the filtered bam and the ptc exon junctions file
            bmo.intersect_bed(bam_file, global_exon_junctions_file, output_file = global_intersect_bam, intersect_bam = True)
            gen.get_time(start)

        global_start = time.time()

        #3: filter to relevant exon-exon junctions
        print("Filtering .bam to relevant exon-exon junctions...")
        ##Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
        #this uses intersect bed, with the intersect bam parameter
        start = time.time()

        intersect_bam = "{0}/{1}_exon_junction_bam_intersect.bam".format(proc_folder, bam_file_parts[1][:-4])
        #intersect the filtered bam and the ptc exon junctions file
        bmo.intersect_bed(global_intersect_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)
        gen.get_time(start)

        #1: count how many reads there are in the sample after filtering to relevant exon-exon junctions but before quality filtering
        read_count_junctions_no_filter = int(gen.run_process(["samtools", "view", "-c", intersect_bam]))


        #4. filter .bam alignments by quality.
        #takes both upper and lower bam thresholds
        #outputs bam file with "_quality_filter_{lower_lim}_{upper_lim}" appended
        # need to do this twice and merge, so we use both intervals used by Geuvadis


        #set the mapq filter parameters here
        mapq_intervals = [[251, 255], [175, 181]]
        mapq_filter_filelist = []

        print("Filtering by MAPQ...")

        for mapq_interval in mapq_intervals:
            start = time.time()
            print(mapq_interval)
            lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
            mapq_filter_file = "{0}/{1}_mapq_filter_{2}_{3}.bam".format(proc_folder, bam_file_parts[1][:-4], lower_threshold, upper_threshold)
            mapq_filter_filelist.append(mapq_filter_file)
            ##run the mapq filter
            bmo.bam_quality_filter(intersect_bam, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)
            gen.get_time(start)

        ##merge files in filelist
        print("Merging bams...")
        start = time.time()
        bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)
        gen.get_time(start)

        print("Filtering by .bam flags...")

        ##filter by flags: get all mapped reads
        #Leaves: mapped unpaired and paired reads
        start = time.time()
        bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)
        gen.get_time(start)

        ##filter bam by xt tag XT=U
        print("xt filter...")
        start = time.time()
        bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")
        gen.get_time(start)

        ##filter bam by nm tag NM<=6
        print("nm filter...")
        start = time.time()
        bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)
        gen.get_time(start)

        #5. scale down the initial count of reads in the sample by the proportion lost during quality filtering
        read_count_junctions_filter = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))
        prop_kept = read_count_junctions_filter / read_count_junctions_no_filter
        read_count = prop_kept * read_count
                                                                                                                  
        print("Phasing reads...")
        #convert to sam format and phase reads
        start = time.time()
        intersect_sam = "{0}_phased.sam".format(mapq_flag_xt_nm_filtered_bam[:-4])
        temp_snp_file = "temp_data/snps{0}.txt".format(random.random())
        so.merge_and_header(PTC_file, syn_nonsyn_file, temp_snp_file)
        bmo.phase_bams(temp_snp_file, mapq_flag_xt_nm_filtered_bam, sample_name, intersect_sam)
        gen.remove_file(temp_snp_file)
        gen.get_time(start)

        print("Counting reads at junctions...")
        start = time.time()
        #6. count the number of reads supporting either the skipping or the inclusion of each exon
        junctions = bmo.read_exon_junctions(PTC_exon_junctions_file)
        bmo.count_junction_reads(intersect_sam, junctions, output_file, read_count)
        gen.get_time(start)

        print("TOTAL:")
        gen.get_time(global_start)

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix", "bam_analysis_folder", "number_of_simulations", "simulation_output_folder", "motif_file", "filter_genome_data", "get_SNPs", "process_bams", "simulate_ptc_snps", "motif_complement", "overwrite_intersect"], flags = [10, 11, 12, 13, 14, 15], ints = [7])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix, bam_analysis_folder, number_of_simulations, simulation_output_folder, motif_file, filter_genome_data, get_SNPs, process_bams, simulate_ptc_snps, motif_complement, overwrite_intersect = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix, args.bam_analysis_folder, args.number_of_simulations, args.simulation_output_folder, args.motif_file, args.filter_genome_data, args.get_SNPs, args.process_bams, args.simulate_ptc_snps, args.motif_complement, args.overwrite_intersect

    start = time.time()

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
    PTC_file = "{0}_ptc_file.txt".format(out_prefix)
    syn_nonsyn_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)
    CDS_interval_file = "{0}_intervals{1}".format(os.path.splitext(CDS_fasta)[0], os.path.splitext(CDS_fasta)[1])
    #check which individuals were included in Geuvadis
    full_sample_names = os.listdir(bams_folder)
    full_sample_names = [i for i in full_sample_names if i[-4:] == ".bam" and "proc" not in i]
    sample_names = [(i.split("."))[0] for i in full_sample_names]
    sample_names = [i for i in sample_names if len(i) > 0]
    #for some reason, 17 of the samples from Geuvadis don't appear in the 1000genomes vcf
    #I'm gonna have to get to the bottom of this at some point
    #but at the moment I'm just gonna filter them out
    with open("../source_data/samples_in_vcf.txt") as file:
        samples_in_vcf = file.readlines()
    samples_in_vcf = [i.rstrip("\n") for i in samples_in_vcf]
    sample_names = [i for i in sample_names if i in samples_in_vcf]
    sample_file = "{0}_sample_file.txt".format(out_prefix)

    if get_SNPs:
        #get SNPs for the sample
        print("Getting SNP data...")
        so.get_snps_in_cds(coding_exon_bed, CDS_bed, vcf_folder, panel_file, sample_names, sample_file, SNP_file, out_prefix)
        gen.get_time(start)
        print("Determining SNP type...")
        so.get_snp_change_status(SNP_file, CDS_fasta, PTC_file, syn_nonsyn_file)
        gen.get_time(start)

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    print("Filtering exon-exon junctions to only leave those that flank exons with a PTC variant...")
    PTC_exon_junctions_file = "{0}_filtered_exon_junctions.bed".format(out_prefix)
    bo.filter_exon_junctions(exon_junctions_file, PTC_file, PTC_exon_junctions_file)

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = ["{0}/{1}".format(bams_folder, i) for i in full_sample_names if (i.split("."))[0] in sample_names]

    #in parallel, do the processing on individual .bam files
    if bam_analysis_folder == "None":
        bam_analysis_folder = "{0}_bam_analysis".format(out_prefix)
    gen.create_directory(bam_analysis_folder)
    if process_bams:
        print("Processing RNA-seq data...")
        #we have to do it like this because you can't pass flags into run_in_parallel
        keyword_dict = {"overwrite_intersect": overwrite_intersect}
        processes = gen.run_in_parallel(bam_files, ["foo", exon_junctions_file, PTC_exon_junctions_file, bam_analysis_folder, PTC_file, syn_nonsyn_file, out_prefix, keyword_dict], process_bam_per_individual)
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

    print("Calculating PSI...")
    final_file = "{0}_final_output.txt".format(out_prefix)
    bmo.compare_PSI(PTC_file, bam_analysis_folder, final_file)

    #run the simulation that swaps ptcs for nonsynonymous snps
    if simulate_ptc_snps:
        if simulate_ptc_snps and not number_of_simulations:
            print("Please specify the number of simulations")
            raise Exception
        ptc_snp_simulation(out_prefix, simulation_output_folder, PTC_file, syn_nonsyn_file, exon_junctions_file, bam_files, number_of_simulations)


if __name__ == "__main__":
    main()
