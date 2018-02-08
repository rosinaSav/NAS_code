import bed_ops as bo
import bam_ops as bmo
import generic as gen
import os

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
    args = gen.parse_arguments(description, ["gtf", "genome_fasta", "out_prefix"])
    gtf, genome_fasta, out_prefix = args.gtf, args.genome_fasta, args.out_prefix

    #extract and filter CDS coordinates and sequences
    #NB! We need to add in a step during quality control to only leave one transcript isoform per gene!
    #and also remove overlapping genes
    CDS_fasta = "{0}_CDS.fasta"
    bo.extract_cds(gtf, CDS_fasta, genome_fasta, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)
    print("Extracted and filtered CDSs.")

    #group the CDS sequences into families based on sequence similarity
    print("Grouping sequences into families...")
    gen.find_families(CDS_fasta, out_prefix)

    #extract exon coordinates
    bo.extract_exons(gtf, exon_bed)

    #only leave exons from transcripts that passed quality control in the extract_cds step above.
    bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed)
    print("Extracted and filtered exons.")


    #extract exon-exon junction coordinates
    bo.extract_exon_junctions(filtered_exon_bed, exon_junctions_file, window_of_interest = 2)
    print("Extracted exon-exon junctions.")

    #make another exons bed that only contains fully coding exons.
    #This is because in the final analysis we should only consider fully protein-coding exons.
    #However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
    #be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
    bo.check_coding(filtered_exon_bed, CDS_bed, filtered_coding_exon_bed)
    print("Filtered out non-coding and partially coding, as well as terminal exons.")
        
    '''
    **********
    MISSING: out of the exons retained in the previous step, check which ones contain a PTC in some individuals and not others.
    NB! Only consider individuals included in Geuvadis.
    **********
    '''

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    bo.filter_exon_junctions(exon_junctions_file, PTC_exons_file, PTC_exon_junctions_file)
    print("Filtered exon-exon junctions to only leave those that flank exons with a PTC variant.")

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = os.listdir(bam_folder)
    bam_files = ["{0}/{1}".format(bam_folder, i) for i in bam_files if i[-4:] == ".bam"]

    #in parallel, do the processing on individual .bam files
    #fill in 'other_arguments' later
    processes = gen.run_in_parallel(bam_files, ["foo", other_arguments], process_bam_per_individual)
    for process in processes:
        process.get()

    print("Processed RNA-seq data.")

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
