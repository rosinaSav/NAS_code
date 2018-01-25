import csv
import generic as gen
import random

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
