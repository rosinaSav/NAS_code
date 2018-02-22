import bed_ops as bo
import csv
import generic as gen
import os
import random
import re
import time
import copy

def bam_flag_filter(input_bam, output_bam, get_paired_reads=None, get_unpaired_reads=None, get_mapped_reads=None, get_unmapped_reads=None, get_mate_mapped_reads=None):
    '''
    Filters bam reads by flag.
    get_paired_reads: get all reads that are tagged as paired
    get_unpaired_reads: get all reads that are not tagged as paired
    get_mapped_reads: get all reads that are not tagged unmapped
    get_unmapped_reads: get all reads that are tagged unmapped
    '''
    samtools_args = ["samtools", "view", "-h"]

    include_flags = []
    exclude_flags = []

    if get_paired_reads:
        #add 1 to include for paired reads, code -f 1
        include_flags.append(1)
    if get_unpaired_reads:
        #add 1 to exclude for unpaired reads, code -F 1
        exclude_flags.append(1)
    if get_mapped_reads:
        #get reads not unmapped, i.e. get mapped reads, code -F 4
        exclude_flags.append(4)
    if get_unmapped_reads:
        #get reads unmapped, code -f 4
        include_flags.append(4)

    include_flags = sum(include_flags)
    exclude_flags = sum(exclude_flags)


    if include_flags > 0:
        samtools_args.extend(["-f", include_flags])
    if exclude_flags > 0:
        samtools_args.extend(["-F", exclude_flags])

    samtools_args.append(input_bam)
    gen.run_process(samtools_args, file_for_output=output_bam)


def bam_quality_filter(input_bam, output_bam, quality_greater_than_equal_to=None, quality_less_than_equal_to=None):
    '''
    Filters bam reads by quality.
    quality_less_than_equal_to: the lower threshold for quality control
    quality_greater_than_equal_to: the upper threshold for quality control
    '''

    samtools_args = ["samtools", "view", "-h"]
    #if neither thresholds are specified
    if not quality_greater_than_equal_to and not quality_less_than_equal_to:
        print("You must specify one threshold to filter reads by.")
        raise Exception
    #if both thresholds are specified
    if quality_greater_than_equal_to and quality_less_than_equal_to:
        #create temp file
        gen.create_directory("temp_data/")
        temp_file = "temp_data/{0}.{1}.bam".format(os.path.split(output_bam)[1][:-4], random.random())
        #first get everything below the upper threshold
        #need to account for the fact samtools removes everything below threshold
        #so when inversing need to add 1 to total
        args = samtools_args.copy()
        upper_limit = quality_less_than_equal_to + 1
        args.extend(["-q", upper_limit, input_bam, "-U", temp_file])
        gen.run_process(args)
        #second get everything above the lower threshold
        args = samtools_args.copy()
        args.extend(["-bq", quality_greater_than_equal_to, temp_file])
        gen.run_process(args, file_for_output=output_bam)
        # #cleanup files
        gen.remove_file(temp_file)
    #if only the lower threshold is specified
    elif quality_greater_than_equal_to and not quality_less_than_equal_to:
        samtools_args.extend(["-bq", quality_greater_than_equal_to, input_bam])
        gen.run_process(samtools_args, file_for_output=output_bam)
    #if only the upper threshold is specified
    elif quality_less_than_equal_to and not quality_greater_than_equal_to:
        #need to account for the fact samtools removes everything below threshold
        #so when inversing need to add 1 to total
        upper_limit = quality_less_than_equal_to + 1
        samtools_args.extend(["-q", upper_limit, input_bam, "-U", output_bam])
        gen.run_process(samtools_args)

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

def intersect_bed(bed_file1, bed_file2, use_bedops = False, overlap = False, overlap_rec = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, no_name_check = False, no_dups = True, chrom = None, intersect = False, hit_count = False, bed_path = None, intersect_bam=None):
    '''Use bedtools/bedops to intersect coordinates from two bed files.
    Return those lines in bed file 1 that overlap with intervals in bed file 2.
    OPTIONS
    output_file: write output to this file
    use_bedops: use bedops rather than bedtools. Certain options are only valid with one of the two, see below.
    overlap: minimum overlap required as a fraction of the intervals in bed file 1 (EX: 0.8 means that the
    overlap has to be at least 80% of the intervals in bed file 1).
    overlap_rec: require that the overlap as a fraction of the intervals in file 2 be at least as high as
    the threshold indicated in -f.
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
    hit_count: for each element in bed file 1, return the number of elements it overlaps in bed file 2 (only valid with bedtools)
    intersect_bam: intersect a bam file with a bed file. Requires bam file to be called first.'''
    gen.create_directory("temp_data/")
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    if use_bedops:
        bedtools_output = run_bedops(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, no_dups = no_dups, intersect_bam = intersect_bam, overlap_rec = overlap_rec)
    else:
        bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, bed_path = bed_path, intersect_bam = intersect_bam, overlap_rec = overlap_rec)
    #move it to a permanent location only if you want to keep it
    if output_file:
        gen.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = bedtools_output.splitlines()
        bedtools_output = [i.split("\t") for i in bedtools_output]
    gen.remove_file(temp_file_name)
    return(bedtools_output)

def retrieve_bams(ftp_site, local_directory, remote_directory, password_file, subset = None):
    '''
    For each .bam file at the ftp site, downsload it, transfer it to
    a remote server and delete it.
    ftp_site: the remote site that contains the files
    local_directory: the local directory where you want to temporarily store the files
    remote_directory: path to directory on remote server where you want to transfer the files
    password_file: path to file that contains Watson password
    subset: only retrieve this many .bam files (useful for testing)
    '''
    #create local directory, if it doesn't exist
    gen.create_directory(local_directory)
    #split the ftp_site address into host and the path
    ftp_site = ftp_site.split("/")
    host = ftp_site[0]
    ftp_directory = "/".join(ftp_site[1:])
    user = "anonymous"
    password = "rs949@bath.ac.uk"
    #connect to FTP server
    ftp = gen.ftp_connect(host, user, password, directory = ftp_directory)
    #get list of all .bam files
    all_files = ftp.nlst()
    all_files = [i for i in all_files if i[-4:] == ".bam"]
    print(len(all_files))
    ftp = gen.ftp_check(ftp, host, user, password, ftp_directory)
    ftp.quit()
    #get password for Watson
    with open(password_file) as file:
        expect_password = "".join(file)
        expect_password = expect_password.rstrip("\n")
    #I will use expect to run scp from the script
    #the way this works is you write an expect script
    #and then use the expect programme to run it
    #this is the string that will be in the script
    #each time, you replace "foo" with the name of the file you want to transfer
    expect_string = "#!/usr/bin/expect\nset timeout -1\nspawn rsync {0}/foo {1}\nexpect \"rs949@bssv-watson's password:\"\nsend \"{2}\\n\";\nexpect eof\nexit".format(local_directory, remote_directory, expect_password)
    if subset:
        all_files = all_files[:subset]
    #retrieve and transfer .bams in parallel
    processes = gen.run_in_parallel(all_files, ["foo", local_directory, host, user, password, ftp_directory, expect_string], retrieve_bams_core, workers = 6)
    for process in processes:
        process.get()

def retrieve_bams_core(all_files, local_directory, host, user, password, ftp_directory, expect_string):
    '''
    Core function parallelized in retrieve_bams above.
    '''
    #connect to FTP server
    ftp = gen.ftp_connect(host, user, password, directory = ftp_directory)
    #loop over .bam files
    for pos, bam_file in enumerate(all_files):
        expect_file = "temp_data/expect_file{0}.txt".format(random.random())
        start_time = time.time()
        print("{0}/{1}".format(pos, len(all_files)))
        local_bam_file = "{0}/{1}".format(local_directory, bam_file)
        #retrieve current file
        if not os.path.isfile(local_bam_file):
            ftp = gen.ftp_retrieve(ftp, host, user, password, ftp_directory, bam_file, destination = local_directory)
        #transfer file to Watson
        current_expect_string = str.replace(expect_string, "foo", bam_file)
        with open(expect_file, "w") as e_file:
            e_file.write(current_expect_string)
        gen.run_process(["expect", expect_file])
        print("Transferred to Watson.")
        gen.remove_file(expect_file)
        gen.remove_file(local_bam_file)
        print("Time spent: {0} minutes.\n".format(round((time.time() - start_time)/60), 3))
    ftp = gen.ftp_check(ftp, host, user, password, ftp_directory)
    ftp.quit()

def run_bedops(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, output_file = None, intersect = False, hit_number = None, no_dups = False, overlap_rec = None, intersect_bam = None):
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
        raise Exception
    if no_dups:
        print("Bedops doesn't print duplicates by default!")
    if overlap_rec:
        print("Bedops hasn't been set up to filter by overlap in second file!")
    if intersect_bam:
        print("Use bedtools to intersect bam and bed!")
        raise Exception
    bedops_output = gen.run_process(bedops_args, file_for_output = output_file)
    return(bedops_output)

def run_bedtools(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False, bed_path = None, overlap_rec = None, intersect_bam = None):
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
    if intersect:
        del bedtools_args[-1]
    if overlap:
        bedtools_args.extend(["-f", str(overlap)])
    if overlap_rec:
        bedtools_args.append("-r")
    if force_strand:
        bedtools_args.append("-s")
    if no_name_check:
        bedtools_args.append("-nonamecheck")
    if no_dups:
        bedtools_args.append("-u")
    if chrom:
        print("Bedtools cannot be restricted to a single chromosome. Use bedops!")
        raise Exception
    if hit_number and no_dups:
        print("When counting hits, each interval in the first bed file is only reported once by default. Set no_dups to False!")
        raise(Exception)
    if bed_path:
        bedtools_args[0] = "{0}{1}".format(bed_path, bedtools_args[0])
    if intersect_bam:
        if A_file[-4:] != ".bam":
            print("Bam file must be called first")
            raise Exception
        if B_file[-4:] != ".bed":
            print("Bed file must be called second")
            raise Exception
        bedtools_args = ["intersectBed", "-wa", "-abam", A_file, "-b", B_file]
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
