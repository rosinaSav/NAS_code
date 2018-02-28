import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import SNP_ops as so
import time


##simulation

#have 2 files, 1 containing all PTCs, 1 containing all other SNPs
#1get a file containing all nonsyn snp
#2generate a psuedo PTC file containing a choice of non synonymous snps with ancestral allele and derived allele matched (need to think about matching allele frequency)
#3map new psudo PTC to exon junctions
#4convert to bam
#5haplotype stuff
#6count reads and output to output folder
#7then run compare psi




def process_bam_per_individual(bam_files, PTC_exon_junctions_file, out_folder, PTC_file, syn_nonsyn_file):
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

        #set the mapq filter parameters here
        mapq_intervals = [[251, 255], [175, 181]]
        mapq_filter_filelist = []

        print("Filtering by MAPQ and merging bams...")

        for mapq_interval in mapq_intervals:
            lower_threshold, upper_threshold = mapq_interval[0], mapq_interval[1]
            mapq_filter_file = "{0}_mapq_filter_{1}_{2}_proc.bam".format(bam_file[-4:], lower_threshold, upper_threshold)
            mapq_filter_filelist.append(mapq_filter_file)
            ##run the mapq filter
            bmo.bam_quality_filter(bam_file, mapq_filter_file, quality_greater_than_equal_to=lower_threshold, quality_less_than_equal_to=upper_threshold)

        ##merge files in filelist
        bam_file_parts = os.path.split(bam_file)
        #setup the filename for the mapq filtered bam file
        mapq_filtered_bam = "{0}/{1}_filtered_mapq_proc.bam".format(bam_file_parts[0], bam_file_parts[1])
        #bmo.merge_bams(mapq_filter_filelist, mapq_filtered_bam)

        print("Filtering by .bam flags...")

        ##filter by flags: get all mapped reads
        #Leaves: mapped unpaired and paired reads
        mapq_flag_filtered_bam = "{0}_flag_proc.bam".format(mapq_filtered_bam[:-4])
        bmo.bam_flag_filter(mapq_filtered_bam, mapq_flag_filtered_bam, get_mapped_reads=True)

        ##filter bam by xt tag XT=U
        mapq_flag_xt_filtered_bam = "{0}_xt_proc.bam".format(mapq_flag_filtered_bam[:-4])
        bmo.bam_xt_filter(mapq_flag_filtered_bam, mapq_flag_xt_filtered_bam, xt_filter="U")

        ##filter bam by nm tag NM<=6
        mapq_flag_xt_nm_filtered_bam = "{0}_nm_proc.bam".format(mapq_flag_xt_filtered_bam[:-4])
        bmo.bam_nm_filter(mapq_flag_xt_filtered_bam, mapq_flag_xt_nm_filtered_bam, nm_less_equal_to=6)

        #extra1: We can then get a count of the total reads possible which can be used for normalisation
        read_count = int(gen.run_process(["samtools", "view", "-c", mapq_flag_xt_nm_filtered_bam]))

        print("Filtering .bam to relevant sequence intervals...")
        ##3. Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
		#this uses intersect bed, with the intersect bam parameter
        intersect_bam = "{0}_exon_junction_filtered_bam_intersect_proc.bam".format(mapq_flag_xt_nm_filtered_bam[-4:])
        bmo.intersect_bed(mapq_flag_xt_nm_filtered_bam, PTC_exon_junctions_file, output_file = intersect_bam, intersect_bam = True)

        #convert to sam format
        intersect_sam = "{0}.sam".format(intersect_bam[-4:])
        gen.run_process(["samtools", "view", intersect_bam], file_for_output = intersect_sam)

        '''
        **********
        MISSING
        **********
        '''
        #when counting junction reads, if the individual is heterozygous, get the big tabix samples file
        #which you have previously intersected with the .bam
        #for all the SNPs in the sample, check which allele the read has
        #assign read to one of the two haplotypes based on that

        print("Counting reads at junctions...")
        #4. count the number of reads supporting either the skipping or the inclusion of each exon
        junctions = bmo.read_exon_junctions(PTC_exon_junctions_file)
        bmo.count_junction_reads(intersect_sam, junctions, "{0}/{1}.txt".format(out_folder, bam_file.split(".")[0]), read_count)

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix"])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix

    start = time.time()

    #extract and filter CDS coordinates and sequences
    CDS_fasta = "{0}_CDS.fasta".format(out_prefix)
    CDS_bed = "{0}_CDS.bed".format(out_prefix)
    print("Extracting and filtering CDSs...")
    #bo.extract_cds(gtf, CDS_bed, CDS_fasta, genome_fasta, all_checks = True, uniquify = True, clean_chrom_only = True, full_chr_name = True)
    gen.get_time(start)

    #group the CDS sequences into families based on sequence similarity
    print("Grouping sequences into families...")
    names = gen.read_fasta(CDS_fasta)[0]
    #gen.find_families_ensembl("../source_data/GRCh37_ensembl_protein_families.txt", names, "{0}_families.txt".format(out_prefix))
    gen.get_time(start)

    print("Extracting and filtering exons...")
    #extract exon coordinates
    exon_bed = "{0}_exons.bed".format(out_prefix)
    #bo.extract_exons(gtf, exon_bed)
    #only leave exons from transcripts that passed quality control in the extract_cds step above.
    #also only leave a single gene per family
    filtered_exon_bed = "{0}_filtered_exons.bed".format(out_prefix)
    #bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed, families_file = "{0}_families.txt".format(out_prefix))
    gen.get_time(start)

    #extract exon-exon junction coordinates
    print("Extracting exon-exon junctions...")
    exon_junctions_file = "{0}_exon_junctions.bed".format(out_prefix)
    #bo.extract_exon_junctions(exon_bed, exon_junctions_file, window_of_interest = 2)
    gen.get_time(start)

    #make another exons bed that only contains fully coding exons.
    #This is because in the final analysis, we should only consider fully protein-coding exons.
    #However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
    #be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
    print("Filtering out overlapping, non-coding and partially coding, as well as terminal exons...")
    coding_exon_bed = "{0}_coding_exons.bed".format(out_prefix)
    #bo.check_coding(filtered_exon_bed, CDS_bed, coding_exon_bed, remove_overlapping = True)
    gen.get_time(start)

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
    SNP_file = "{0}_SNP_file.txt".format(out_prefix)
    PTC_file = "{0}_ptc_file.txt".format(out_prefix)
    syn_nonsyn_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)
    #get SNPs for the sample
    print("Getting SNP data...")
    CDS_interval_file = "{0}_intervals{1}".format(os.path.splitext(CDS_fasta)[0], os.path.splitext(CDS_fasta)[1])
    #so.get_snps_in_cds(coding_exon_bed, CDS_bed, vcf_folder, panel_file, sample_names, sample_file, SNP_file, out_prefix)
    gen.get_time(start)
    print("Determining SNP type...")
##    so.get_snp_change_status(SNP_file, CDS_fasta, PTC_file, syn_nonsyn_file)
    gen.get_time(start)

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    print("Filtering exon-exon junctions to only leave those that flank exons with a PTC variant...")
    PTC_exon_junctions_file = "{0}_filtered_exon_junctions.bed".format(out_prefix)
##    bo.filter_exon_junctions(exon_junctions_file, PTC_file, PTC_exon_junctions_file)

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = ["{0}/{1}".format(bams_folder, i) for i in full_sample_names if (i.split("."))[0] in sample_names]

    #in parallel, do the processing on individual .bam files
    bam_analysis_folder = "{0}_bam_analysis".format(out_prefix)
    gen.create_directory(bam_analysis_folder)
    print("Processing RNA-seq data...")
    processes = gen.run_in_parallel(bam_files, ["foo", PTC_exon_junctions_file, bam_analysis_folder, PTC_file, syn_nonsyn_file], process_bam_per_individual)
    for process in processes:
        process.get()

    print("Calculating PSI...")
    final_file = "{0}_final_output.txt".format(out_prefix)
    bmo.compare_PSI(PTC_file, bam_analysis_folder, final_file)


if __name__ == "__main__":
    main()
