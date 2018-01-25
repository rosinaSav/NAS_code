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
