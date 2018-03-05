import bed_ops as bo
import csv
import generic as gen
import numpy as np
import os
import random
import re
import time
import copy

def bam_flag_filter(input_bam, output_bam, get_paired_reads=None, get_unpaired_reads=None, get_mapped_reads=None, get_unmapped_reads=None, get_proper_paired_reads=None, get_improper_paired_reads=None, get_mate_mapped_reads=None, get_mate_unmapped_reads=None):
    '''
    Filters bam reads by flag.
    get_paired_reads: get all reads from paired-end (or multiple-segment) sequencing technology
    get_unpaired_reads: get all reads that not from paired-end (or multiple-segment) sequencing technology
    get_mapped_reads: get all reads where the segment is mapped
    get_unmapped_reads: get all reads where the segment is unmapped
    get_proper_paired_reads: get all reads that where each segment is properly aligned according to the aligner
    get_improper_paired_reads: get all reads that where each segment is not properly aligned according to the aligner
    get_mate_mapped_reads: get all reads where the next segment in the template is mapped
    get_mate_unmapped_reads: get all reads where the next segment in the template is unmapped
    '''
    samtools_args = ["samtools", "view", "-bh"]

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
    if get_proper_paired_reads:
        if not get_paired_reads:
            print("get_paired_reads must be set!")
            raise Exception
        #get all reads where each segment is properly aligned, code -f 2
        include_flags.append(2)
    if get_improper_paired_reads:
        if not get_paired_reads:
            print("get_paired_reads must be set!")
            raise Exception
        #get all reads where each segment is not properly aligned, code -F 2
        exclude_flags.append(2)
    if get_mate_mapped_reads:
        if not get_paired_reads:
            print("get_paired_reads must be set!")
            raise Exception
        #get all reads where the next segment in the template is not unmapped, i.e. mapped, code -F 8
        exclude_flags.append(8)
    if get_mate_unmapped_reads:
        if not get_paired_reads:
            print("get_paired_reads must be set!")
            raise Exception
        #get all reads where the next segment in the template is unmapped, code -f 8
        include_flags.append(8)

    #get the final flag code for filter
    include_flags = sum(include_flags)
    exclude_flags = sum(exclude_flags)
    #include arguments
    if include_flags > 0:
        samtools_args.extend(["-f", include_flags])
    if exclude_flags > 0:
        samtools_args.extend(["-F", exclude_flags])
    #run samtools with flag arguments
    samtools_args.append(input_bam)
    gen.run_process(samtools_args, file_for_output=output_bam)

def bam_nm_filter(input_bam, output, nm_less_equal_to=None):
    '''
    Filters bam reads by NM value.
    nm_less_equal_to: the NM value you wish to filter by.
    '''
    if not nm_less_equal_to:
        print("Please provide NM filter value.")
        raise Exception

    #create output file
    if output[-4:] == ".bam":
        output_file = "{0}.sam".format(output[:-4])
    else:
        output_file = output
    sam_output = gen.run_process(["samtools", "view", "-h", input_bam])
    #create grep args and include header fields if they exist
    grep_args = ["^@"]
    #for each nm less than equal to threshold, create grep arg
    for i in range(nm_less_equal_to+1):
        grep_args.append("\|\tNM:i:{0}\t".format(i))
    grep_args = "".join(grep_args)
    gen.run_process(["grep", grep_args], input_to_pipe=sam_output, file_for_output = output_file)

    #if wanting to create bam, create bam and delete sam
    if output != output_file:
        samtools_args = ["samtools", "view", "-bh", output_file]
        gen.run_process(samtools_args, file_for_output=output)
        gen.remove_file(output_file)

def bam_xt_filter(input_bam, output, xt_filter=None):
    '''
    Filter a bam/sam file by XT tag.
    '''
    if not xt_filter:
        print("Please specify XT filter.")
        raise Exception
    #create output file
    if output[-4:] == ".bam":
        output_file = "{0}.sam".format(output[:-4])
    else:
        output_file = output


    sam_output = gen.run_process(["samtools", "view", "-h", input_bam])
    grep_args = []
    #get header lines
    grep_args.append("^@")
    #get XT values with xt_filter
    grep_args.append("\|\tXT:A:{0}\t".format(xt_filter))
    grep_args = "".join(grep_args)
    gen.run_process(["grep", grep_args], input_to_pipe=sam_output, file_for_output = output_file)

    #if wanting to create bam, create bam and delete sam
    if output != output_file:
        samtools_args = ["samtools", "view", "-bh", output_file]
        gen.run_process(samtools_args, file_for_output=output)
        gen.remove_file(output_file)

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

def compare_PSI(SNP_file, bam_folder, out_file, round_norm_count = None):
    '''
    Given PTC-generating SNPs, as well as read counts at exon-exon junctions, compare exon skipping rates
    within samples that do or do not have a PTC within a given exon.
    If round_norm_count is specified, the normalized read counts will be rounded to
    that number of digits after the decimal point.
    '''
    SNPs = gen.read_many_fields(SNP_file, "\t")
    samples = SNPs[0][15:-1]
    #note that if two SNPs appear in the same exon, the one that appears later
    #will overwrite the one that appears first so only one of the SNPs will be analyzed
    SNPs = {i[3]: i[15:-1] for i in SNPs[1:]}
    results = {i: {"PSI_w_PTC": [], "PSI_het_PTC": [], "PSI_no_PTC": [], "norm_count_w_PTC": [],
                   "norm_count_het_PTC": [], "norm_count_no_PTC": [], "ptc_count": 0,
                   "sample_count": 0} for i in SNPs}
    for pos, sample in enumerate(samples):
        file_name = "{0}/{1}.txt".format(bam_folder, sample)
        #in case you're analyzing simulation output
        if not os.path.isfile(file_name):
            file_name = "{0}/{1}_simulation_{1}.txt".format(bam_folder, sample)
        #maybe that bam simply hasn't been processed yet, skip it if that's the case
        if os.path.isfile(file_name):
            with open(file_name) as file:
                for line in file:
                    line = line.split("\t")
                    exon = line[0]
                    #this also filters out the header line
                    if exon in SNPs:
                        skipped = np.sum([int(i) for i in line[1].split("|")])
                        included = np.sum([int(i) for i in line[2].split("|")])
                        total = int(line[3])
                        genotype = SNPs[exon][pos]
                        if skipped > 0 or included > 0:
                            results[exon]["sample_count"] = results[exon]["sample_count"] + 1
                            #if this sample contains a PTC
                            if "1" in genotype:
                                #if it's heterozygous:
                                if "0" in genotype:
                                    results[exon]["ptc_count"] = results[exon]["ptc_count"] + 0.5                               
                                    results[exon]["PSI_het_PTC"].append(included/(skipped + included))
                                    results[exon]["norm_count_het_PTC"].append(skipped/total)
                                #if it's homozygous for PTC
                                else:
                                    results[exon]["ptc_count"] = results[exon]["ptc_count"] + 1                               
                                    results[exon]["PSI_w_PTC"].append(included/(skipped + included))
                                    results[exon]["norm_count_w_PTC"].append(skipped/total)
                            #if it's homozygous for lack of PTC
                            else:
                                results[exon]["PSI_no_PTC"].append(included/(skipped + included))
                                results[exon]["norm_count_no_PTC"].append(skipped/total)
    header = "exon\tptc_count\tsample_count\tPSI_w_PTC\tPSI_het_PTC\tPSI_no_PTC\tnorm_count_w_PTC\tnorm_count_het_PTC\tnorm_count_no_PTC\n"
    header_split = header.split("\t")
    header_split[-1] = header_split[-1].rstrip("\n")
    with open(out_file, "w") as file:
        file.write(header)
        for exon in sorted(results):
            if results[exon]["sample_count"] > 0:
                file.write("{0}\t".format(exon))
                #:-1 cause you don't want a \t at the end of the line
                for info in header_split[1:-1]:
                    to_write = results[exon][info]
                    if type(to_write) == list:
                        to_write = np.mean(to_write)
                        if round_norm_count:
                            to_write = round(to_write, round_norm_count)
                    file.write("{0}\t".format(to_write))
                to_write = results[exon][header_split[-1]]
                to_write = np.mean(to_write)
                if round_norm_count:
                    to_write = round(to_write, round_norm_count)
                file.write("{0}\n".format(to_write))

def compare_PSI_haplotypes(SNP_file, bam_folder, out_file):
    '''
    Given PTC-generating SNPs, as well as read counts at exon-exon junctions, compare exon skipping rates
    within samples that do or do not have a PTC within a given exon. Consider the two haplotypes in a sample as separate data points.
    Discard information from unphased reads.
    '''
    SNPs = gen.read_many_fields(SNP_file, "\t")
    samples = SNPs[0][15:-1]
    #note that if two SNPs appear in the same exon, the one that appears later
    #will overwrite the one that appears first so only one of the SNPs will be analyzed
    SNPs = {i[3]: i[15:-1] for i in SNPs[1:]}
    results = {i: {"PSI_w_PTC": [], "PSI_no_PTC": [], "norm_count_w_PTC": [],
                   "norm_count_no_PTC": [], "ptc_count": 0,
                   "sample_count": 0} for i in SNPs}
    for pos, sample in enumerate(samples):
        with open("{0}/{1}.txt".format(bam_folder, sample)) as file:
            for line in file:
                line = line.split("\t")
                exon = line[0]
                #this also filters out the header line
                if exon in SNPs:
                    skipped = [int(i) for i in line[1].split("|")]
                    included = [int(i) for i in line[2].split("|")]
                    total = int(line[3])
                    genotype = SNPs[exon][pos].split("|")
                    #we're ignoring data from unphased reads (third subfield in the skipped and included fields)
                    for haplotype in range(2):
                        #if any phased reads mapped to this exon-exon junction
                        if skipped[haplotype] > 0 or included[haplotype] > 0:
                            results[exon]["sample_count"] = results[exon]["sample_count"] + 1
                            #if this sample contains a PTC
                            if genotype[haplotype] == "1":
                                results[exon]["ptc_count"] = results[exon]["ptc_count"] + 1                               
                                results[exon]["PSI_w_PTC"].append(included[haplotype]/(skipped[haplotype] + included[haplotype]))
                                results[exon]["norm_count_w_PTC"].append(skipped[haplotype]/total)
                            else:
                                results[exon]["PSI_no_PTC"].append(included[haplotype]/(skipped[haplotype] + included[haplotype]))
                                results[exon]["norm_count_no_PTC"].append(skipped[haplotype]/total)
    header = "exon\tptc_count\tsample_count\tPSI_w_PTC\tPSI_no_PTC\tnorm_count_w_PTC\tnorm_count_no_PTC\n"
    header_split = header.split("\t")
    header_split[-1] = header_split[-1].rstrip("\n")
    with open(out_file, "w") as file:
        file.write(header)
        for exon in sorted(results):
            if results[exon]["sample_count"] > 0:
                file.write("{0}\t".format(exon))
                #:-1 cause you don't want a \t at the end of the line
                for info in header_split[1:-1]:
                    to_write = results[exon][info]
                    if type(to_write) == list:
                        to_write = round(np.mean(to_write), 3)
                    file.write("{0}\t".format(to_write))
                to_write = results[exon][header_split[-1]]
                to_write = round(np.mean(to_write), 3)
                file.write("{0}\n".format(to_write))


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

def count_junction_reads(sam, junctions, outfile, read_count):
    '''
    Given a sam file and a dictionary of exon-exon junctions, count how many reads overlap each junction.
    For each exon, count how many reads support its skipping and how many support its inclusion.
    Multiply the former count by 2.
    '''
    out_dict = {}
    haplotypes = ["0", "1", "N"]
    with open(sam) as file:
        for line in file:
            line = line.split("\t")
            sam_start = int(line[3])
            cigar = line[5]
            chrom = line[2]
            #0, 1 or N
            haplotype = (line[0].split("|"))[1]
            if chrom in junctions:
                #get intron position in chromosome coordinates based on the alignment cigar
                putative_junctions = map_intron_from_cigar(cigar, sam_start)
                if putative_junctions:
                    for junction in putative_junctions:
                        #if the 3'end of the exon is in the junctions dict
                        if junction[0] in junctions[chrom]:
                            #if the 5' end of the exon is in the junctions dict
                            if junction[1] in junctions[chrom][junction[0]]:
                                current_dict = junctions[chrom][junction[0]][junction[1]]
                                for pos, exon in enumerate(current_dict["exon"]):
                                    if exon not in out_dict:
                                        #we have three counters instead of one because we're counting the two copies separately
                                        #and there's the third category for reads that couldn't be phased
                                        out_dict[exon] = {"skip": [0, 0, 0], "incl": [0, 0, 0]}
                                    haplotype_pos = haplotypes.index(haplotype)
                                    out_dict[exon][current_dict["type"][pos]][haplotype_pos] = out_dict[exon][current_dict["type"][pos]][haplotype_pos] + 1
            else:
                print("Chromosome {0} not found!".format(chrom))
    with open(outfile, "w") as file:
        file.write("exon\tskippedx2\tincluded\ttotal_reads\n")
        for exon in sorted(out_dict):
            file.write("{0}\t{1}\t{2}\t{3}\n".format(exon, "|".join([str(i * 2) for i in out_dict[exon]["skip"]]), "|".join([str(i) for i in out_dict[exon]["incl"]]), read_count))

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
                             force_strand = False, no_name_check = False, no_dups = True, chrom = None, intersect = False, hit_count = False, bed_path = None, intersect_bam=None,
                  write_zero = False, write_bed = False):
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
    intersect_bam: intersect a bam file with a bed file. Requires bam file to be called first
    write_zero: like write_both but also write A intervals that don't overlap with any B intervals,
    write_bed: when intersecting a bam file, write output as bed.'''
    gen.create_directory("temp_data/")
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    if use_bedops:
        bedtools_output = run_bedops(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, no_dups = no_dups, intersect_bam = intersect_bam, overlap_rec = overlap_rec)
    else:
        bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, bed_path = bed_path, intersect_bam = intersect_bam, write_zero = write_zero, overlap_rec = overlap_rec, write_bed = write_bed)
    #move it to a permanent location only if you want to keep it
    if output_file:
        gen.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = gen.read_many_fields(temp_file_name, "\t")
    gen.remove_file(temp_file_name)
    return(bedtools_output)

def map_from_cigar(cigar, sam_start, focal_pos):
    '''
    Given the cigar and the start coordinate (base 1) of an exon-exon read in a bam file,
    as well as the (base 1) coordinate of a focal position, determine where the focal
    position falls in the read sequence (in base 0). 
    '''
    consume_query = ["M", "I", "S"]
    consume_ref = ["M", "D", "N"]
    split_cigar = re.findall("\d+\w", cigar)
    ref_count = sam_start
    query_count = 0
    for part in split_cigar:
        if ref_count == focal_pos:
            return(query_count)
        letter = part[-1]
        number = int(part[:-1])
        for occurrence in range(number):
            if ref_count == focal_pos:
                return(query_count)
            if letter in consume_query:
                query_count = query_count + 1
            if letter in consume_ref:
                ref_count = ref_count + 1

def map_intron_from_cigar(cigar, sam_start):
    '''
    Given the cigar and the start coordinate (base 1) of an exon-exon read in a bam file,
    determine the (base 0) coordinates of both the 3'-most nucleotide in the 5' exon and the 5'-most
    nucleotide in the 3' exon. Can handle reads that map multiple exon-exon junctions.
    '''
    result = []
    #check if one gap in alignment
    gaps = [i for i in re.finditer("\d+N", cigar)]
    if len(gaps) < 1:
        return None
    for gap in gaps:
        #intron length
        length = int(gap.group(0)[:-1])
        #where the intron starts in the alignment
        before_gap = cigar[:gap.start()]
        pairs = re.findall("\d+(?=[MDN])", before_gap)
        rel_start = 0
        if pairs:
            pairs = [int(i) for i in pairs]
            rel_start = np.sum(pairs)
        #where the intron starts on the chromosome
        abs_start = (sam_start - 1) + rel_start - 1
        #where the intron ends (actually the base after that)
        abs_end = abs_start + length + 1
        result.append([abs_start, abs_end])
    return(result)

def merge_bams(bam_list, output_file):
    '''
    Merge a list of bam files to defined output file.
    '''
    #setup args, add -r to attach filename rg tag
    args = ["samtools", "merge", "-r"]
    if os.path.exists(output_file):
        args.append("-f")
    args.append(output_file)
    #loop through each input file and add to argument list
    for file in bam_list:
        args.append(file)
    gen.run_process(args)

def phase_bams(snps, bam, sample_name, out_sam):
    '''
    Given a bed file of SNPs and a sam file of reads, annotate the reads according to haplotype.
    '''
    temp_snps = "temp_data/{0}.bed".format(random.random())
    gen.run_process(["cut", "-f7-", snps], file_for_output = temp_snps)
    intersecting_reads = intersect_bed(bam, temp_snps, write_zero = True, intersect_bam = True, write_bed = True)
    sam = (gen.run_process(["samtools", "view", bam])).split("\n")
    sam = [i.split("\t") for i in sam]
    #store cigar and sequence info cause that is lost in the bam to bed conversion
    #you can't just put down the name of the read cause those aren't unique
    sam = {" ".join([i[0], i[3]]): [i[5], i[9]] for i in sam if len(i) > 1 and i[0][0] != "#"}
    #get order of samples in SNP file
    samples = ((((gen.run_process(["head", "-2", temp_snps])).split("\n"))[1]).split("\t"))[9:-1]
    sample_pos = samples.index(sample_name)
    
    out_dict = {}
    for read in intersecting_reads:
        #because the bam to bed conversion appends read pair member identifier
        read_name = read[3].split("/")[0]
        read_code = "|".join([read_name, read[1], read[0], read[5]])
        if read_code not in out_dict:
            out_dict[read_code] = []
        #read overlaps with a SNP
        if read[14] != ".":
            genotypes = read[21:-2]
            genotype = genotypes[sample_pos]
            #a site where this individual is homozygous won't tell us anything
            if "0" in genotype and "1" in genotype:
                sam_read_name = " ".join([read_name, str(int(read[1]) + 1)])
                cigar, seq = sam[sam_read_name]
##                print("\n")
##                print(read)
##                print(cigar)
##                print(seq)
##                print(map_from_cigar(cigar, int(read[1]) + 1, int(read[13])))
##                print(seq[map_from_cigar(cigar, int(read[1]) + 1, int(read[13]))])
                current_base = seq[map_from_cigar(cigar, int(read[1]) + 1, int(read[13]))]
                reference = read[15]
                variant = read[16]
                genotype = genotype.split("|")
                #don't count it if the base in the read is neither the reported reference nor the reported variant base
                if current_base == reference:
                    haplotype = genotype.index("0")
                    out_dict[read_code].append(haplotype)
                elif current_base == variant:
                    haplotype = genotype.index("1")
                    out_dict[read_code].append(haplotype)
    with open(out_sam, "w") as file:
        for read in sorted(out_dict):
            count_0 = out_dict[read].count(0)
            count_1 = out_dict[read].count(1)
            if count_0 > count_1:
                haplotype = 0
            elif count_1 > count_0:
                haplotype = 1
            else:
                haplotype = "N"
            read = read.split("|")
            #I'm using this strange format to maintain compatibility with count_junction_reads
            #(that is also why I'm taking the read start position to base 1)
            file.write("{0}|{1}\t.\t{2}\t{3}\t.\t{4}\n".format(read[0], haplotype, read[2], int(read[1]) + 1, sam[" ".join([read[0], str(int(read[1]) + 1)])][0]))
    gen.remove_file(temp_snps)

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

def read_exon_junctions(junctions_file):
    '''
    Read the exon junctions in a bed file into a dictionary that tells you how the junctions are matched
    and what information they contain on skipping events.
    '''
    out_dict = {}
    work_dict = {}
    ends_list = []
    junctions = gen.read_many_fields(junctions_file, "\t")
    #loop over intervals, only consider 3' exon ends on first go
    for junction in junctions:
        if len(junction) > 1:
            chrom = junction[0]
            if chrom not in out_dict:
                out_dict[chrom] = {}
            strand = junction[5]
            start_id = "5"
            end_id = "3"
            #we're gonna ignore strand in all of the analysis to come
            #so we need antisense exon-exon junctions to look like sense exon-exon junctions
            if strand == "-":
                start_id = "3"
                end_id = "5"
            #if it's a 3' exon end, just create the record in the output dictionary
            if junction[3][-1] == end_id:
                start = int(junction[1])
                if start not in out_dict[chrom]:
                    out_dict[chrom][start] = {}
                work_dict[junction[3]] = start
            #otherwise put it aside
            elif junction[3][-1] == start_id:
                ends_list.append(junction)

    #now go over the 5' exon ends and match them up with the 3' exon ends
    for junction in ends_list:
        name = junction[3].split(".")
        trans = name[0]
        exon = int(name[1])
        strand = junction[5]
        end_id = "3"
        mult_factor = 1
        if strand == "-":
            end_id = "5"
            mult_factor = -1
        #the corresponding 3' exon end can either be that of the previous exon or of that of
        #the exon before it
        #for exons on the antisense strand, you will be looking at exons to come rather than preceding exons,
        #hence the use of mult_factor
        #note that the append might lead to an exon number of -1
        #but expected_starts like that won't appear in work_dict so they'll be skipped
        expected_starts = ["{0}.{1}.{2}".format(trans, exon + (-1 * mult_factor), end_id)]
        expected_starts.append("{0}.{1}.{2}".format(trans, exon + (-2 * mult_factor), end_id))
        end = int(junction[1])
        chrom = junction[0]
        for expected_start in expected_starts:
            if expected_start in work_dict:
                start = work_dict[expected_start]
                current_exon = int(expected_start.split(".")[1])
                #if it's the previous exon, reads overlapping that junction support the splicing in of
                #both this and the current exon
                if abs(exon - current_exon) == 1:
                    exons = ["{0}.{1}".format(trans, current_exon), "{0}.{1}".format(trans, exon)]
                    types = ["incl", "incl"]
                #if it's the 3'end of n-2th exon, then reads overlapping that junction support the skipping
                #of n-1th exon
                elif abs(exon - current_exon) == 2:
                    exons = ["{0}.{1}".format(trans, current_exon + (1 * mult_factor))]
                    types = ["skip"]
                out_dict[chrom][start][end] = {"exon": exons, "type": types}
    return(out_dict)

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

def run_bedtools(A_file, B_file, force_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False, bed_path = None, overlap_rec = None, intersect_bam = None, write_zero = None, write_bed = False):
    '''
    See intersect_bed for details.
    '''
    if write_both:
        write_option = "-wo"
    elif hit_number:
        write_option = "-c"
    elif write_zero:
        write_option = "-wao"
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
        bedtools_args = ["intersectBed", write_option, "-abam", A_file, "-b", B_file]
        if write_bed:
            bedtools_args.append("-bed")
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
