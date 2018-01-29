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
                             force_strand = False, no_name_check = False, no_dups = True, chrom = None, bed_input = False, intersect = False, hit_count = False):
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
        bedtools_output = run_bedops(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, output_file = temp_file_name, intersect = intersect, hit_number = hit_count)
    else:
        bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_number)
    #move it to a permanent location only if you want to keep it
    if output_file:
        gen.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = bedtools_output.splitlines()
        bedtools_output = [i.split("\t") for i in bedtools_output]
    gen.remove_file(temp_file_name)
    return(bedtools_output)

def run_bedops(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, output_file = None, intersect = False, hit_number = hit_number):
    '''
    See intersect_bed for details.
    '''
    if intersect:
        command = "--intersect"
    else:
        command = "--element-of"
    if sort:
        sorting_temp_file_name = "temp_data/temp_sort_{1}.bed".format(random.random())
        sorting_temp_file_name2 = "temp_data/temp_sort2_{1}.bed".format(random.random())
        gen.run_process(["sort-bed", A_file], file_for_output = sorting_temp_file_name)
        gen.run_process(["sort-bed", B_file], file_for_output = sorting_temp_file_name2)
        bedops_args = ["bedops", "--chrom", "foo", command, "1", sorting_temp_file_name, sorting_temp_file_name2]
    else:
        bedops_args = ["bedops", "--chrom", "foo", command, "1", A_file, B_file]
    if overlap:
        bedops_args[4] = overlap
    if chrom:
        bedops_args[2] = chrom
    else:
        del bedops_args[1:3]
    if intersect:
        del bedops_args[4]
    if force_strand:
        print("Bedops can't search by strand! Either use bedtools or separate input data by strand!")
        raise Exception
    if write_both:
        print("Bedops can't write both features!")
        raise Exception
    if hit_number:
        print("Bedops hasn't been set up to count the number of overlapping elements. Use bedtools!")
    bedops_output = gen.run_process(bedops_args, file_for_output = output_file)
    if sort:
        gen.remove_file(sorting_temp_file_name)
        gen.remove_file(sorting_temp_file_name2)
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
        sorting_temp_file_name = "temp_data/temp_sort_{1}.bed".format(random.random())
        sorting_temp_file_name2 = "temp_data/temp_sort2_{1}.bed".format(random.random())
        gen.run_process(["sort-bed", A_file], file_for_output = sorting_temp_file_name)
        gen.run_process(["sort-bed", B_file], file_for_output = sorting_temp_file_name2)
        bedtools_args = ["intersectBed", "-a", sorting_temp_file_name,"-b", sorting_temp-file_name2, write_option]
    else:
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
    bedtools_output = gen.run_process(bedtools_args, file_for_output = output_file)
    return(bedtools_output)
