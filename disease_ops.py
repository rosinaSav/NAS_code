'''
Authors: Rosina Savisaar and Liam Abrahams.
Contains more specific functions for operations related the processing of disease data.
'''

import bed_ops as beo
import bam_ops as bao
import generic as gen
import os
import random
import collections

def check_line(line, outlist):
    '''
    Check the existence of line features and send to list for output
    '''

    exist_list = ["Chromosome", "Start_position", "dbSNP_RS", "Strand", "Reference_Allele", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Matched_Norm_Sample_Barcode", "Transcript_Strand", "NCBI_Build", "t_ref_count", "t_alt_count"]
    match_list = [["Start_position", "End_position"]]
    not_empty_list = ["Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Matched_Norm_Sample_Barcode", "t_ref_count", "t_alt_count"]

    strands = ["+", "-", "1", "-1"]

    # check that the required entries exist
    for entry in exist_list:
        if entry not in line:
            return False, None
    # check not empty
    for entry in not_empty_list:
        if line[entry] in ["-", ".", ""]:
            return False, None
    # check the strands exist
    if line["Strand"] not in strands or line["Transcript_Strand"] not in strands:
        return False, None
    # check that the items in a pair exist and they match
    for pair in match_list:
        if pair[0] not in line or pair[1] not in line or line[pair[0]] != line[pair[1]]:
            return False, None
    # check the genome build
    if line["NCBI_Build"] != "37":
        return False, None

    # write the entry
    entry_out = []
    for entry in outlist:
        if entry in line and line[entry] not in ["", " "]:
            if entry in ["Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"] and line["Transcript_Strand"] == "-":
                # *****
                # After lots of checking, for some reason the alleles
                # are listed as on the + strand in the files
                # but are actually on the correct strands. I dont think we need
                # to convert here

                # if on the - strand, need to read on the - strand
                base = line[entry]
                # base = gen.reverse_complement(base)
                entry_out.append(base)
            else:
                entry_out.append(line[entry])
        else:
            entry_out.append('.')

    return True, entry_out

def concatenate_files(filelist, output_file):
    gen.remove_file(output_file)
    args = ["cat"]
    for file in filelist:
        args.append(file)
    gen.run_process(args, file_for_output = output_file)

def refactor_files(dir, output_dir, filename_prefix, full_mutation_file, limit=None, subset=None, subset_no=None, clean_directory=None):
    '''
    Refactor the files so they resemble something like
    a vcf file that we can use.
    dir: the directory containing the mutation files
    output_dir: the directory to output the processed files
    filename_prefix: prefix for the processed files
    '''

    if clean_directory:
        # clean the directory
        [gen.remove_file("{0}/{1}".format(output_dir, file)) for file in os.listdir(output_dir)]

    counter = 0
    passed = 0
    useable_count = 0
    temp_file_list = []
    outlist = [
        "Chromosome", "Start_position", "dbSNP_RS", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification", "Variant_Type",
        "Refseq_mRNA_Id", "Refseq_prot_Id", "Strand", "Transcript_Strand", "t_ref_count", "t_alt_count"
    ]

    current_temp_file = "temp_data/{0}_{1}.txt${2}".format(filename_prefix, len(temp_file_list)+1, random.random())
    temp_file_list.append(current_temp_file)
    outfile = open(current_temp_file, "w")

    with open("{0}/{1}_readme.txt".format(output_dir, filename_prefix), "w") as readme:
        readme.write("{0}\n".format("\t".join(outlist)))

    mutations_filelist = os.listdir(dir)
    if subset:
        mutations_filelist = mutations_filelist[:int(subset_no)]

    # remove any extras
    mutations_filelist = [file for file in mutations_filelist if '.DS_Store' not in file]


    print("{0} mutation files...".format(len(mutations_filelist)))

    for i, file in enumerate(mutations_filelist):
        filepath = "{0}/{1}".format(dir, file)
        mutations = gen.read_many_fields(filepath, "\t")

        file_count = 0
        # read the header line
        headers = mutations[0]
        # for each of the mutations, write to a dictionary to query
        for mut in mutations[1:]:
            line_dict = {}
            for j, entry in enumerate(mut):
                line_dict[headers[j]] = entry

            # check that the required information is present
            passed_checks, entry_out = check_line(line_dict, outlist)
            if passed_checks:
                file_count += 1
                passed += 1
                # write to file, and reset counter if we want to split the files up
                outfile.write("{0}\n".format("\t".join(entry_out)))
                if limit:
                    counter = gen.update_reset_count(counter, limit)
                    if counter == 0:
                        outfile.close()
                        gen.run_process(["cp", current_temp_file, "{0}/{1}".format(output_dir, current_temp_file.split("$")[0].split('/')[-1])])
                        current_temp_file = "temp_data/{0}_{1}.txt${2}".format(filename_prefix, len(temp_file_list)+1, random.random())
                        temp_file_list.append(current_temp_file)
                        outfile = open(current_temp_file, "w")
        if file_count > 0:
            useable_count += 1

    print("{0} files with useable mutations...".format(useable_count))

    # close the last output file
    outfile.close()
    gen.run_process(["cp", current_temp_file, "{0}/{1}".format(output_dir, current_temp_file.split("$")[0].split('/')[-1])])
    print("{0} sub files generated...".format(len(temp_file_list)))

    temp_file = "temp_data/{0}.txt".format(random.random())
    temp_file_list = [temp_file] + temp_file_list
    with open(temp_file, "w") as temp_out:
        temp_out.write("##fileformat=VCFv4.1\n")
        temp_out.write("#{0}\n".format("\t".join(outlist)))

    # get one resulting file
    concatenate_files(temp_file_list, full_mutation_file)
    # clean the temp directory if any files remain
    [gen.remove_file(file) for file in temp_file_list]

    # sort file by first, then second column
    temp_sorted_file = "temp_data/{0}.vcf".format(random.random())
    gen.run_process(["sort", "-k1,1", "-k2,2n", full_mutation_file], file_for_output = temp_sorted_file)
    gen.run_process(["cp", temp_sorted_file, full_mutation_file])
    gen.remove_file(temp_sorted_file)

    zip = "{0}.gz".format(full_mutation_file)
    gen.run_process(["bgzip", "-c", full_mutation_file], file_for_output = zip)

    print("{0} mutations remain...".format(passed))


def junction_raw_counts_to_bed(file, dir, output_file, output_sample_file):
    '''
    Write the exon junction raw counts to a bed file format.
    files: list of files to process
    dir: directory containing files with raw counts
    output_dir: directory to contain the processed files
    '''

    filepath = "{0}/{1}".format(dir, file)
    # read file
    file_lines = gen.read_many_fields(filepath, "\t")
    # get a list of sample barcodes and the lines
    samples = file_lines[0][1:]
    lines = file_lines[2:]

    # write the samples to a file
    with open(output_sample_file, "w") as sample_outfile:
        sample_outfile.write("{0}\n".format("\t".join(samples)))

    with open(output_file, "w") as outfile:
        # write each of the lines to a bed like format
        for line in lines:
            info = line[0].split(':')
            chr = info[0]
            start = info[1]
            strand = info[2][0]
            stop = info[3]
            line_items = [chr, start, stop, ".", strand]
            line_items.extend(line[1:])
            outfile.write("{0}\n".format("\t".join(line_items)))

def raw_counts_to_samples(intersect_file, sample_file, output_dir):
    '''
    Take a list of samples and write the read counts for each sample
    to its own file
    intersect_file: file containing intersect of exon junctions and raw counts
    sample_file: file containing sample_names
    output_dir: directory of the output folder
    '''

    # read in the sample names
    samples = gen.read_many_fields(sample_file, "\t")[0]

    # create an output directory to hold the samples
    gen.create_output_directories(output_dir)

    # read in the raw counts
    raw_counts = gen.read_many_fields(intersect_file, "\t")

    # for each sample
    for i, sample in enumerate(samples):
        sample_file = "{0}/{1}.bed".format(output_dir, sample)

        with open(sample_file, "w") as outfile:
            # for each of the raw counts, create list with raw count and add to file
            for entry in raw_counts:
                outlist = entry[:11]
                outlist.append(entry[11 + i])
                outfile.write("{0}\n".format("\t".join(outlist)))


def get_exon_locations(file):

    entries = gen.read_many_fields(file, "\t")
    # create dictionaries
    location_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))

    for entry in entries:
        chr = entry[0]
        t = entry[3].split('.')[0]
        e = int(entry[3].split('.')[1])
        start = int(entry[1])
        stop = int(entry[2])
        # create all sites in a transcript
        for i in list(range(start, stop)):
            location_list[chr][i] = [t, e]

    return location_list


def get_full_raw_counts(file):

    raw_counts = {}
    counts = gen.read_many_fields(file, "\t")
    for sample in counts:
        raw_counts[sample[0]] = int(sample[1])
    return(raw_counts)


def process_files(files, dir, output_dir, output_junction_suffix, exon_junctions_file, coding_exons_file, exon_counts_file):
    '''
    Wrapper for processing each raw counts file
    '''

    raw_counts = get_full_raw_counts(exon_counts_file)

    location_list = get_exon_locations(coding_exons_file)

    for file in files:
        print("Processing {0}...".format(file))
        sample = file.split('.')[0]
        file_prefix = "{0}/{1}".format(output_dir, file[:-4])
        # convert the file to a bed like format
        file_bed_format = "{0}.bed".format(file_prefix)
        file_sample_list = "{0}_sample_list.bed".format(file_prefix)
        junction_raw_counts_to_bed(file, dir, file_bed_format, file_sample_list)
        # # intersect with exon junctions
        exon_junction_intersect = "{0}_exon_junction_intersect.bed".format(file_prefix)
        bao.intersect_bed(exon_junctions_file, file_bed_format, output_file=exon_junction_intersect, write_both=True, no_dups = False)

        sample_dir = "{0}/{1}".format(output_dir, sample)
        raw_counts_to_samples(exon_junction_intersect, file_sample_list, sample_dir)

        sample_junctions_output_dir = "{0}/{1}_{2}".format(output_dir, sample, output_junction_suffix)
        count_junction_reads(sample_dir, sample_junctions_output_dir, raw_counts, location_list)

def process_counts(dir, output_dir, output_junction_suffix, exon_junctions_file, results_prefix, exon_counts_file):
    '''
    Wrapper for processing reads
    dir: directory containing files with raw counts
    output_dir: directory to contain the processed files
    '''
    gen.create_output_directories(output_dir)
    # get a list of all the files
    filelist = [file for file in os.listdir(dir) if file != ".DS_Store"]
    # run the processing
    coding_exons_file = "{0}_coding_exons.bed".format(results_prefix)

    # gen.run_in_parallel(filelist, ["foo", dir, output_dir, output_junction_suffix, exon_junctions_file, coding_exons_file], process_files)
    process_files(filelist, dir, output_dir, output_junction_suffix, exon_junctions_file, coding_exons_file, exon_counts_file)


def count_junction_reads(sample_dir, output_dir, full_counts, location_list):
    '''
    Given a sample file and a dictionary of exon-exon junctions, count how many reads overlap each junction.
    For each exon, count how many reads support its skipping and how many support its inclusion.
    Multiply the former count by 2.
    '''

    sex_chr = ["chrX", "chrY"]

    gen.create_output_directories(output_dir)

    out_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: [0,0]))

    samples = os.listdir(sample_dir)

    for sample in samples:
        sample_name = sample[:-4]
        filepath = "{0}/{1}".format(sample_dir, sample)
        lines = gen.read_many_fields(filepath, "\t")

        for line in lines:
            chrom = line[6]
            start = int(line[7])
            stop = int(line[8])
            raw_count = int(line[-1])

            if chrom not in sex_chr:
                if len(location_list[chrom][start]) and len(location_list[chrom][stop]):
                    transcript = location_list[chrom][start][0]
                    exon1 = location_list[chrom][start][1]
                    exon2 = location_list[chrom][stop][1]
                    if abs(exon1-exon2) == 1:
                        out_dict[transcript][exon1][0] += raw_count
                        out_dict[transcript][exon2][0] += raw_count
                    else:
                        # add 1 to the start of the skipped exons to exlucde the first exon
                        # and get all inbetween
                        skipped_exons = list(range(min(exon1, exon2)+1, max(exon1, exon2)))
                        for j in skipped_exons:
                            out_dict[transcript][j][1] += raw_count

        outfile = "{0}/{1}".format(output_dir, sample)
        with open(outfile, "w") as file:
            file.write("exon\tskippedx2\tincluded\ttotal_reads\n")
            for transcript in sorted(out_dict):
                for exon in sorted(out_dict[transcript]):
                    file.write("{0}.{1}\t{2}|0|0\t{3}|0|0\t{4}\n".format(transcript, exon, out_dict[transcript][exon][1]*2, out_dict[transcript][exon][0], full_counts[sample_name]))




def check_ptcs(ptc_file, processed_dir, processed_suffix):

    # create a dictionary of the samples and their corresponding files
    processed_dirs = [dir for dir in os.listdir(processed_dir) if processed_suffix in dir]
    processed_filelist = {}
    for dir in processed_dirs:
        dir = "{0}/{1}".format(processed_dir, dir)
        filelist = os.listdir(dir)
        for file in filelist:
            if file != ".DS_Store":
                processed_filelist[file.strip('.bed')[:15]] = "{0}/{1}".format(dir, file)


    ptcs = gen.read_many_fields(ptc_file, "\t")

    ptcs_with_samples = []

    for ptc in ptcs:
        t_sample = ptc[15][:15]

        if t_sample in processed_filelist:
            ptcs_with_samples.append(ptc)

    for i, ptc in enumerate(ptcs_with_samples):

        if i:
            transcript = ptc[3].split('.')[0]
            t_sample = ptc[15][:15]
            n_sample = ptc[18][:15]

            t_sample_file = processed_filelist[t_sample]
            t_samples = gen.read_many_fields(t_sample_file, "\t")


<<<<<<< HEAD
            for sample in t_samples:
                sample_transcript = sample[0].split('.')[0]
                if transcript == sample_transcript:
                    incl_count = int(sample[2].split('|')[0])
                    skip_count = int(sample[1].split('|')[0])
=======
def process_raw_counts(input_dir, output_file):
    '''
    Process the raw counts for all exons for each sample.
    input_dir: the directory containing the raw count files.
    '''
>>>>>>> 5b2ee73730f001f8748a2ab82b4ef5904724f0ef

                    if incl_count != 0 and skip_count != 0:
                        print(ptc)
                        print(sample)

<<<<<<< HEAD
                    # for n_sample in n_samples:
                    #     if n_sample[0] == exon:
                    #         print(sample, n_sample)
=======
    temp_filelist = []

    for i, file in enumerate(filelist):
        print('Processing {0}/{1}: {2}...'.format(i, len(filelist), file))
        samples_raw_counts = collections.defaultdict(lambda: 0)

        filepath = "{0}/{1}".format(input_dir, file)
        with open(filepath) as infile:
            lines = infile.readlines()

            samples_list = lines[0].split('\t')[1:][::3]
            info_list = lines[1].split('\t')[1:][::3]

            for line in lines[2:]:
                raw_counts = line.split('\t')[1:][::3]
                for j, count in enumerate(raw_counts):
                    samples_raw_counts[samples_list[j]] += int(count)

        temp_file = "temp_data/{0}".format(random.random())
        with open(temp_file, "w") as temp:
            for sample in samples_raw_counts:
                temp.write("{0}\t{1}\n".format(sample, samples_raw_counts[sample]))
        temp_filelist.append(temp_file)

    # now write all raw counts to one file
    all_samples = collections.defaultdict(lambda: 0)
    with open(output_file, "w") as outfile:
        for file in temp_filelist:
            counts = gen.read_many_fields(file, "\t")
            for sample in counts:
                all_samples[sample[0]] += int(sample[1])

        for sample in all_samples:
            outfile.write("{0}\t{1}\n".format(sample, all_samples[sample]))

    for file in temp_filelist:
        gen.remove_file(file)
>>>>>>> 5b2ee73730f001f8748a2ab82b4ef5904724f0ef
