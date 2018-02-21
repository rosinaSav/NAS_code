import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import SNP_ops as so
import time

def process_bam_per_individual(bam_files, PTC_exon_junctions_file, out_folder):
    #add other arguments

    for bam_file in bam_files:

		##Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
		#this uses intersect bed, with the intersect bam paramter
		bmo.intersect_bed(PTC_exon_junctions_file, bam_file, output_file="{0}_exon_junction_intersect.bed".format(bam_file[-4:]), intersect_bam=True)

        '''
        **********
        MISSING: filter remaining .bam alignments by quality.
        **********
        '''

        '''
        **********
        MISSING: Find cases where both of the reads in a pair have been retained in step 1, and randomly remove one of them.
        **********
        '''

                #count the number of reads supporting either the skipping or the inclusion of each exon
                junctions = bmo.read_exon_junctions(PTC_exon_junctions_file)
                bmo.count_junction_reads(clean_sam, junctions, "{0}/{1}.txt".format(out_folder, bam_file.split(".")[0])

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix"])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix

    start = time.time()

    #extract and filter CDS coordinates and sequences
    CDS_fasta = "{0}_CDS.fasta".format(out_prefix)
    CDS_bed = "{0}_CDS.bed".format(out_prefix)
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
    exon_bed = "{0}_exons.bed".format(out_prefix)
    bo.extract_exons(gtf, exon_bed)
    #only leave exons from transcripts that passed quality control in the extract_cds step above.
    #also only leave a single gene per family
    filtered_exon_bed = "{0}_filtered_exons.bed".format(out_prefix)
    bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed, families_file = "{0}_families.txt".format(out_prefix))
    gen.get_time(start)

    #extract exon-exon junction coordinates
    print("Extracting exon-exon junctions...")
    exon_junctions_file = "{0}_exon_junctions.bed".format(out_prefix)
    bo.extract_exon_junctions(exon_bed, exon_junctions_file, window_of_interest = 2)
    gen.get_time(start)

    #make another exons bed that only contains fully coding exons.
    #This is because in the final analysis, we should only consider fully protein-coding exons.
    #However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
    #be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
    print("Filtering out overlapping, non-coding and partially coding, as well as terminal exons...")
    coding_exon_bed = "{0}_coding_exons.bed".format(out_prefix)
    bo.check_coding(filtered_exon_bed, CDS_bed, coding_exon_bed, remove_overlapping = True)
    gen.get_time(start)

    #check which individuals were included in Geuvadis
    sample_names = os.listdir(bams_folder)
    sample_names = [(i.split("."))[0] for i in sample_names]
    sample_names = [i for i in sample_names if len(i) > 0]
    #for some reason, 17 of the samples from Geuvadis don't appear in the 1000genomes vcf
    #I'm gonna have to get to the bottom of this at some point
    #but at the moment I'm just gonna filter them out
    with open("../source_data/samples_in_vcf.txt") as file:
        samples_in_vcf = file.readlines()
    samples_in_vcf = [i.rstrip("\n") for i in samples_in_vcf]
    sample_names = [i for i in sample_names if i in samples_in_vcf]
    print(len(sample_names))
    sample_names = []
    sample_file = "{0}_sample_file.txt".format(out_prefix)
    SNP_file = "{0}_SNP_file.txt".format(out_prefix)
    PTC_file = "{0}_ptc_file.txt".format(out_prefix)
    syn_nonsyn_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)
    #get SNPs for the sample
    print("Getting SNP data...")
    CDS_interval_file = "{0}_intervals{1}".format(os.path.splitext(CDS_fasta)[0], os.path.splitext(CDS_fasta)[1])
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
    bam_files = ["{0}/{1}".format(bam_folder, i) for i in sample_names if i[-4:] == ".bam"]

    #in parallel, do the processing on individual .bam files
    bam_analysis_folder = "{0}_bam_analysis".format(out_prefix)
    gen.create_directory(bam_analysis_folder)
    processes = gen.run_in_parallel(bam_files, ["foo", PTC_exon_junctions_file, bam_analysis_folder], process_bam_per_individual)
    print("Processing RNA-seq data...")
    for process in processes:
        process.get()

    '''
    **********
    MISSING: for each of the exon-exon junctions that have been retained, determine median PSI (or count/total sample size) for individuals that have/don't have the PSI.
    Would be good to also record sample size for either group.
    Average results across paralogous families.
    Perform some sort of a paired comparison between PTC+ and PTC-.
    **********
    '''

if __name__ == "__main__":
    main()
