import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os
import SNP_ops as so

def process_bam_per_individual(bam_files, other_arguments):
    #add other arguments

    for bam_file in bam_files:

        '''
        **********
        MISSING: code to filter the bam file to only leave split alignments (i.e. ones that overlap an exon-exon junction).
        Information contained in the cigar.
        **********
        '''

        '''
        **********
        MISSING: filter remaining .bam alignments by quality.
        **********
        '''

        '''
        **********
        MISSING: our reads are paired-end, so we need to split the intersect step in two.
        1) Intersect junctions and .bam, and write down the overlapping .bam alignments, without counting.
        Intermediate stage: Find cases where both of the reads in a pair have been retained in step 1, and randomly remove one of them.
        This is because you can theoretically imagine a situation where one read overlaps with one junction and
        the other read maps with another junction. I think it's better to only count one of them because the two reads from the same
        fragment are not independent data points. Do you agree?
        2) An intersect like the one below where you just count the number of overlapping alignments.
        NB! The RNA-seq data we're using is not stranded so contrary to what I had written before, we should NOT match strand in the intersect.
        **********
        '''

        #count how many .bam alignments overlap each exon-exon junction
        bmo.intersect_bed(PTC_exon_junctions_file, bam_file, overlap = 1, output_file = "{0}_junction_hit_count.bed".format(bam_file[:-4]), force_strand = False, no_dups = False, hit_count = True, use_bedops = False)

        '''
        **********
        MISSING: for each fully protein-coding exon,
        1) estimate PSI in this individual (no. of alignments that support exon inclusion/(no. of alignments that support exon inclusion + number of alignments that support exon skipping))
        2) calculate (no. of alignments that support exon inclusion)/(total no. of alignments in sample)
        **********
        '''

def main():

    description = "Check whether PTCs are associated with greater rates of exon skipping."
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "bams_folder", "vcf_folder", "panel_file", "out_prefix"])
    gtf, genome_fasta, bams_folder, vcf_folder, panel_file, out_prefix = args.gtf, args.genome_fasta, args.bams_folder, args.vcf_folder, args.panel_file, args.out_prefix

    start = time.time()

    #extract and filter CDS coordinates and sequences
    CDS_fasta = "{0}_CDS.fasta".format(out_prefix)
    CDS_bed = "{0}_CDS.bed".format(out_prefix)
    print("Extracting and filtering CDSs...")
    #bo.extract_cds(gtf, CDS_bed, CDS_fasta, genome_fasta, all_checks = True, uniquify = True)
    gen.get_time(start)

    #match transcript identifiers with gene names so that you coud also generate
    #a families output file that is meaningful to a human
    descriptions_file = "{0}_descriptions.txt".format(out_prefix)
    names = gen.read_fasta(CDS_fasta)[0]
    print("Matching transcript identifiers with gene names...")
    #bo.get_descriptions(names, gtf, descriptions_file)
    gen.get_time.start()

    #group the CDS sequences into families based on sequence similarity
    print("Grouping sequences into families...")
    names = gen.read_fasta(CDS_fasta)[0]
##    gen.find_families_ensembl("../source_data/GRCh37_ensembl_protein_families.txt", names, "{0}_families.txt".format(out_prefix))
    gen.get_time(start)

    print("Extracting and filtering exons...")
    #extract exon coordinates
    exon_bed = "{0}_exons.bed".format(out_prefix)
##    bo.extract_exons(gtf, exon_bed)
    #only leave exons from transcripts that passed quality control in the extract_cds step above.
    #also only leave a single gene per family
    filtered_exon_bed = "{0}_filtered_exons.bed".format(out_prefix)
##    bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed, families_file = "{0}_families.txt".format(out_prefix))
    gen.get_time(start)

    #extract exon-exon junction coordinates
    print("Extracting exon-exon junctions...")
    exon_junctions_file = "{0}_exon_junctions.bed".format(out_prefix)
##    bo.extract_exon_junctions(exon_bed, exon_junctions_file, window_of_interest = 2)
    gen.get_time(start)

    #make another exons bed that only contains fully coding exons.
    #This is because in the final analysis, we should only consider fully protein-coding exons.
    #However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
    #be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
    print("Filtering out overlapping, non-coding and partially coding, as well as terminal exons...")
    coding_exon_bed = "{0}_coding_exons.bed".format(out_prefix)
##    bo.check_coding(filtered_exon_bed, CDS_bed, coding_exon_bed, remove_overlapping = True)
    gen.get_time(start)
        
    #check which individuals were included in Geuvadis
    sample_names = os.listdir(bams_folder)
    sample_names = [(i.split("."))[0] for i in sample_names]
    sample_names = [i for i in sample_names if len(i) > 0]
    sample_file = "{0}_sample_file.txt".format(out_prefix)
    SNP_file = "{0}_SNP_file.txt".format(out_prefix)
    PTC_file = "{0}_ptc_file.txt".format(out_prefix)
    syn_nonsyn_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)   
    #get SNPs for the sample
    print("Getting SNP data...")
    so.get_snps_in_cds(coding_exon_bed, CDS_fasta, vcf_folder, panel_file, sample_names, sample_file, SNP_file, out_prefix)
    gen.get_time(start)
    print("Determining SNP type...")
    so.get_snp_change_status(SNP_file, CDS_fasta, PTC_file, syn_nonsyn_file)
    gen.get_time(start)

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    print("Filtering exon-exon junctions to only leave those that flank exons with a PTC variant...")
    bo.filter_exon_junctions(exon_junctions_file, PTC_file, PTC_exon_junctions_file)

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = os.listdir(bam_folder)
    bam_files = ["{0}/{1}".format(bam_folder, i) for i in bam_files if i[-4:] == ".bam"]

    #in parallel, do the processing on individual .bam files
    #fill in 'other_arguments' later
    processes = gen.run_in_parallel(bam_files, ["foo", other_arguments], process_bam_per_individual)
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
