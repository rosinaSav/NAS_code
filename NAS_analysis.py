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
        NB! I believe the sequencing was paired-end so at some point, we need to figure out what that
        means for our data and what to do with that information.
        **********
        '''

        '''
        **********
        MISSING: filter remaining .bam alignments by quality.
        **********
        '''
        
        #count how many .bam alignments overlap each exon-exon junction
        bmo.intersect_bed(PTC_exon_junctions_file, bam_file, overlap = 1, output_file = "{0}_junction_hit_count.bed".format(bam_file[:-4]), force_strand = True, no_dups = False, hit_count = True, use_bedops = False)

        '''
        **********
        MISSING: for each fully protein-coding exon,
        1) estimate PSI in this individual (no. of alignments that support exon inclusion/(no. of alignments that support exon inclusion + number of alignments that support exon skipping))
        2) calculate (no. of alignments that support exon inclusion)/(total no. of alignments in sample)
        **********
        '''

def main():

    '''
    **********
    MISSING: Parse input arguments.
    **********
    '''

    #extract and filter CDS coordinates and sequences
    #NB! We need to add in a step during quality control to only leave one transcript isoform per gene!
    #and also remove overlapping genes
    bo.extract_cds(gtf, CDS_fasta, genome_fasta, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)

    '''
    **********
    MISSING: code that will group sequences into paralogous families.
    **********
    '''

    #extract exon coordinates
    bo.extract_exons(gtf, exon_bed)

    #only leave exons from transcripts that passed quality control in the extract_cds step above.
    bo.filter_bed_from_fasta(exon_bed, CDS_fasta, filtered_exon_bed)

    #extract exon-exon junction coordinates
    bo.extract_exon_junctions(filtered_exon_bed, exon_junctions_file, window_of_interest = 2)

    '''
    **********
    MISSING: make another exons bed that only contains fully coding exons (easy to do by intersecting with CDS file).
    This is because in the final analysis we should only consider fully protein-coding exons.
    However, for getting the exon junctions we need the full exons file because fully protein-coding exons might
    be flanked by exons that are not. This is why we couldn't do this filtering step earlier.
    **********
    '''

    '''
    **********
    MISSING: out of the exons retained in the previous step, check which ones contain a PTC in some individuals and not others.
    NB! Only consider individuals included in Geuvadis.
    **********
    '''

    #filter the exon junctions file to only leave those junctions that flank exons retained in the previous step.
    bo.filter_exon_junctions(exon_junctions_file, PTC_exons_file, PTC_exon_junctions_file)

    #make a list of all the .bam files and modify them to have the full path rather than just the file name
    bam_files = os.listdir(bam_folder)
    bam_files = ["{0}/{1}".format(bam_folder, i) for i in bam_files if i[-4:] == ".bam"]

    #in parallel, do the processing on individual .bam files
    #fill in 'other_arguments' later
    processes = gen.run_in_parallel(bam_files, ["foo", other_arguments], process_bam_per_individual)
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
