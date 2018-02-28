import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import random
import SNP_ops as so
import time

def process_bam_per_individual(bam_files, PTC_exon_junctions_file, out_folder, PTC_file, syn_nonsyn_file, filter_bams):
    #add other arguments

    for bam_file in bam_files:

        #Process:
        # 1. Filter bams
        # This gives us a set of "good" quality reads.
        # extra1: Get a count of all reads
        # 2. Map to exon junctions
        # 3. Count reads either skipping or including each exon

        #####1. filter .bam alignments by quality.
        # upper and lower bam thresholds required
        # need to do this twice and merge, so we use both intervals used by Geuvadis

        bam_file_parts = os.path.split(bam_file)
        mapq_filtered_bam = "{0}/{1}_filtered_mapq_proc.bam".format(bam_file_parts[0], bam_file_parts[1])
        mapq_flag_filtered_bam = "{0}_flag_proc.bam".format(mapq_filtered_bam[:-4])
        mapq_flag_xt_filtered_bam = "{0}_xt_proc.bam".format(mapq_flag_filtered_bam[:-4])
        mapq_flag_xt_nm_filtered_bam = "{0}_nm_proc.bam".format(mapq_flag_xt_filtered_bam[:-4])

        if filter_bams:
            #set the mapq filter parameters here
            mapq_intervals = [[251, 255], [175, 181]]
            mapq_filter_filelist = []

            print("Filtering by MAPQ and merging bams...")

            for mapq_interval in mapq_intervals:
                lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
                mapq_filter_file = "{0}_mapq_filter_{1}_{2}_proc.bam".format(bam_file[:-4], lower_threshold, upper_threshold)
                mapq_filter_filelist.append(mapq_filter_file)
                ##run the mapq filter
                bmo.bam_quality_filter(bam_file, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)

            ##merge files in filelist
            bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)

            print("Filtering by .bam flags...")

            ##filter by flags: get all mapped reads
            #Leaves: mapped unpaired and paired reads
            bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)

            ##filter bam by xt tag XT=U
            bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")

            ##filter bam by nm tag NM<=6
            bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)

        #extra1: We can then get a count of the total reads possible which can be used for normalisation
        read_count = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))

        print("Filtering .bam to relevant sequence intervals...")
        ##3. Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
	#this uses intersect bed, with the intersect bam parameter
        intersect_bam = "{0}_exon_junction_filtered_bam_intersect_proc.bam".format(mapq_flag_xt_nm_filtered_bam[:-4])
        bmo.intersect_bed(mapq_flag_xt_nm_filtered_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)

        print("Phasing reads...")
        #convert to sam format and phase reads
        intersect_sam = "{0}_phased.sam".format(intersect_bam[:-4])
        temp_snp_file = "temp_data/snps{0}.txt".format(random.random())
        so.merge_and_header(PTC_file, syn_nonsyn_file, temp_snp_file)
        sample_name = (bam_file.split("/")[-1]).split(".")[0]
        bmo.phase_bams(temp_snp_file, intersect_bam, sample_name, intersect_sam)
        gen.remove_file(temp_snp_file)
        print(intersect_sam)

        print("Counting reads at junctions...")
        #4. count the number of reads supporting either the skipping or the inclusion of each exon
        junctions = bmo.read_exon_junctions(PTC_exon_junctions_file)
        bmo.count_junction_reads(intersect_sam, junctions, "{0}/{1}.txt".format(out_folder, sample_name), read_count)

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix", "bam_analysis_folder", "filter_genome_data", "get_SNPs", "filter_bams"], flags = [7, 8, 9])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix, bam_analysis_folder, filter_genome_data, get_SNPs, filter_bams = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix, args.bam_analysis_folder, args.filter_genome_data, args.get_SNPs, args.filter_bams

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
    print("Processing RNA-seq data...")
    processes = gen.run_in_parallel(bam_files, ["foo", PTC_exon_junctions_file, bam_analysis_folder, PTC_file, syn_nonsyn_file, filter_bams], process_bam_per_individual)
    for process in processes:
        process.get()

    gen.get_time(start)

    print("Calculating PSI...")
    final_file = "{0}_final_output.txt".format(out_prefix)
    bmo.compare_PSI(PTC_file, bam_analysis_folder, final_file)       


if __name__ == "__main__":
    main()
