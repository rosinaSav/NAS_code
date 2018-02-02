import bed_ops as bo
import csv
import generic as gen
import random
import re

def convert2bed(input_file_name, output_file_name, group_flags = None):
    '''
    Converts an input file (sam, bam, gtf, gff...) to a bed file using bedops.
    Set 'group_flags' to an integer if you want to group all the fields from a certain field onwards.
    For instance, if you set group_flags to 5, all of the fields from the 5th onward will be turned
    into a comma-separated string and stored as one field in the.bed file.
    Note that you cannot group all the fields in a row (i.e. you can't set it to 0.)
    '''
    extension = gen.get_extension(input_file_name, 3, ["sam", "bam", "gtf", "gff"])
    bed_data = gen.run_process(["convert2bed", "--input={0}".format(extension)], file_for_input = input_file_name, file_for_output = output_file_name)
    if group_flags:
        temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
        group_flags(output_file_name, temp_file_name, group_flags)
        gen.run_process(["mv", temp_file_name, output_file_name])
        print("Grouped flags.")      
    print("Converted data from {0} to bed.".format(extension))

def extract_exons(gtf, bed):
    '''Given a GTF file, extract exon coordinates and write them to .bed.
    EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf", "../source_data/Homo_sapiens.GRCh37.87_exons.bed")'''
    #extract exons from GTF
    exons = gen.run_process(["grep", "\texon\t", gtf])
    #filter down to only protein-coding ones
    exons = gen.run_process(["grep", "gene_biotype \"protein_coding\""], input_to_pipe = exons)
    #split lines
    exons = [i.split("\t") for i in exons.split("\n")]
    #format as .bed. Switch to base 0.
    exons = [["chr{0}".format(i[0]), int(i[3]) - 1, i[4], i[8], ".", i[6]] for i in exons if len(i) >= 8]
    #pre-compile regex
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
    #extract transcript IDs and exon numbers
    for pos, exon in enumerate(exons):
        to_parse = exon[3]
        trans = re.search(trans_regex, to_parse).group(0)
        exon_no = re.search(exon_no_regex, to_parse).group(0)
        exons[pos][3] = "{0}.{1}".format(trans, exon_no)
    #write to bed
    with open(bed, "w") as file:
        for exon in exons:
            file.write("{0}\n".format("\t".join([str(i) for i in exon])))

def group_flags(input_bed, output_bed, flag_start):
    '''Takes an input bed file and converts all the fields from the flag_start'th
    onwards into a single field, with the elements separated by commas.'''
    with open(input_bed) as input_file, open(output_bed, "w") as output_file:
        reader = csv.reader(input_file, delimiter = "\t")
        writer = csv.writer(output_file, delimiter = "\t")
        for i in reader:
            flags = i[flag_start:]
            flags = ",".join(flags)
            new_row = i[0:flag_start]
            new_row.append(flags)
            writer.writerow(new_row)

def intersect_bed(bed_file1, bed_file2, use_bedops = False, overlap = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, no_name_check = False, no_dups = True, chrom = None, intersect = False, hit_count = False):
    '''Use bedtools/bedops to intersect coordinates from two bed files.
    Return those lines in bed file 1 that overlap with intervals in bed file 2.
    OPTIONS
    output_file: write output to this file
    use_bedops: use bedops rather than bedtools. Certain options are only valid with one of the two, see below.
    overlap: minimum overlap required as a fraction of the intervals in bed file 1 (EX: 0.8 means that the
    overlap has to be at least 80% of the intervals in bed file 1).
    write_both: if True, return not only the interval from bed file 1 but, tagged onto the end, also the
    interval from bed file 2 that it overlaps (only
    valid when using bedtools).
    sort: sort bed files before taking the intersection
    force_strand: check that the feature and the bed interval are on the same strand (only valid with bedtools)
    no_name_check: if set to False, checks whether the chromosome names are the same in the too bed files (only valid with bedtools)
    no_dups: if True, only returns each interval once. If set to false, intervals in bed file 1 that overlap several intervals in
    bed file 2 will be returned several times (as many times as there are overlaps with different elements in bed file 2)
    chrom: limit search to a specific chromosome (only valid with bedops, can help in terms of efficiency)
    intersect: rather than returning the entire interval, only return the part of the interval that overlaps an interval in bed file 2.
    hit_count: for each element in bed file 1, return the number of elements it overlaps in bed file 2 (only valid with bedtools)'''
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    if use_bedops:
        bedtools_output = run_bedops(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, no_dups = no_dups)
    else:
        bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count)
    #move it to a permanent location only if you want to keep it
    if output_file:
        gen.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = bedtools_output.splitlines()
        bedtools_output = [i.split("\t") for i in bedtools_output]
    gen.remove_file(temp_file_name)
    return(bedtools_output)

def run_bedops(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, output_file = None, intersect = False, hit_number = None, no_dups = False):
    '''
    See intersect_bed for details.
    '''
    if intersect:
        command = "--intersect"
    else:
        command = "--element-of"
    if sort:
        sort_bed(A_file, A_file)
        sort_bed(B_file, B_file)
    bedops_args = ["bedops", "--chrom", "foo", command, "1", A_file, B_file]
    if overlap:
        bedops_args[4] = overlap
    if chrom:
        bedops_args[2] = chrom
        if intersect:
            del bedops_args[4]
    else:
        del bedops_args[1:3]
        if intersect:
            del bedops_args[2]
    if force_strand:
        print("Bedops can't search by strand! Either use bedtools or separate input data by strand!")
        raise Exception
    if write_both:
        print("Bedops can't write both features!")
        raise Exception
    if hit_number:
        print("Bedops hasn't been set up to count the number of overlapping elements. Use bedtools!")
    if no_dups:
        print("Bedops doesn't print duplicates by default!")
    print(bedops_args)
    bedops_output = gen.run_process(bedops_args, file_for_output = output_file)
    return(bedops_output)

def run_bedtools(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False):
    '''
    See intersect_bed for details.
    '''
    if write_both:
        write_option = "-wo"
    elif hit_number:
        write_option = "-c"
    else:
        write_option = "-wa"
    if sort:
        sort_bed(A_file, A_file)
        sort_bed(B_file, B_file)
    bedtools_args = ["intersectBed", "-a", A_file,"-b", B_file, write_option]
    if overlap:
        bedtools_args.extend(["-f", str(overlap)])
    if force_strand:
        bedtools_args.append("-s")
    if no_name_check:
        bedtools_args.append("-nonamecheck")
    if no_dups:
        bedtools_args.append("-u")
    if chrom:
        print("Bedtools cannot be restricted to a single chromosome. Use bedops!")
        raise Exception
    if intersect:
        print("Bedtools has not been set up to take intersections. Either implement it or use bedops!")
    if hit_number and no_dups:
        print("When counting hits, each interval in the first bed file is only reported once by default. Set no_dups to False!")
        raise(Exception)
    bedtools_output = gen.run_process(bedtools_args, file_for_output = output_file)
    return(bedtools_output)

def sort_bed(input_file_name, output_file_name):
    '''
    Sort a bed file.
    '''
    #This is done via a temp file because that way you can specify the same file as input and output file and thus
    #overwrite the unsorted file with the sorted one.
    temp_file_name = "temp_data/temp_sorted_bed{0}.bed".format(random.random())
    gen.run_process(["sort-bed", input_file_name], file_for_output = temp_file_name)
    gen.run_process(["mv", temp_file_name, output_file_name])
    gen.remove_file(temp_file_name)

def write_hits_at_junctions_per_sample(ftp_site, target_directory, exon_junctions_file, subset = None):
    '''
    For each .bam file at the ftp site, intersect it with the exon junction intervals and
    make a file with the number of overlapping reads for each junction interval.
    subset: only retrieve this many .bam files (useful for testing)
    '''
    #create target directory, if it doesn't exist
    gen.make_dir(target_directory)
    #split the ftp_site address into host and the path
    ftp_site = ftp_site.split("/")
    host = ftp_site[0]
    ftp_directory = "/".join(ftp_site[1:])
    user = "anonymous"
    password = "rs949@bath.ac.uk"
    #connect to FTP server
    #ftp = gen.ftp_connect(host, user, password, directory = ftp_directory)
    #get list of all .bam files
    #all_files = ftp.nlst()
    all_files = ["NA12399.1.M_120209_4.bam"]
    all_files = [i for i in all_files if i[-4:] == ".bam"]
    if subset:
        all_files = all_files[:subset]
    #loop over .bam files
    for pos, bam_file in enumerate(all_files):
        print("{0}/{1}".format(pos, len(all_files)))
        #retrieve current file
        #ftp = gen.ftp_retrieve(ftp, host, user, password, ftp_directory, bam_file, destination = target_directory)
        #note that overlap = 1
        local_bam_file = "{0}/{1}".format(target_directory, bam_file)
        #gen.run_process(["curl", "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/processed/NA12399.1.M_120209_4.bam"], file_for_output = local_bam_file)
        print("Retrieved file {0}.".format(bam_file))
        #local_bed_file = "temp_data/test.bed"
        #convert2bed(local_bam_file, local_bed_file)
        #test = gen.run_process(["head", "-5", local_bed_file])
        #print(test)
        intersect_bed(exon_junctions_file, local_bam_file, overlap = 1, output_file = "{0}/{1}_junction_hit_count.bed".format(target_directory, bam_file[:-4]), force_strand = False, no_dups = False, hit_count = False, use_bedops = True)
        print("Intersected with exon-exon junctions.\n")
        #gen.remove_file(local_bam_file)
        raise(Exception)
    #end the connection to the FTP server
    #make sure the connection is live before you do or you might crash
    ftp = gen.ftp_check(ftp, host, user, password, ftp_directory)
    ftp.quit()

    
