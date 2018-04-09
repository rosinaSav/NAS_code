'''
Authors: Rosina Savisaar and Liam Abrahams.
Contains functions for operations on sequence coordinates.
Some functions wrap bedtools or bedops, others are all Python.
'''

import generic as gen
import bam_ops as bmo
import re
import collections
import copy
import numpy as np
import os
from pathlib import Path
import random
import shutil


def change_bed_names(input_bed, output_bed, full_names, header):
    '''
    Convert between full chromosome names and just chromosome numbers.
    As is, output_bed needs to be different to input_bed.
    full_names: Set to true if the input_bed has the full chr name, e.g. chr1
    header: Set to true if the bed file has a header
    '''
    lines = gen.read_many_fields(input_bed, "\t")

    # Set whether the bed file has a header or not
    if header:
        start_index = 1
    else:
        start_index = 0

    with open(output_bed, "w") as outfile:
        if header:
            outfile.write("{0}\n".format("\t".join(lines[0])))
        for line in lines[start_index:]:
            if full_names:
                line[0] = line[0].strip("chr")
            else:
                line[0] = "chr{0}".format(line[0])
            outfile.write("{0}\n".format("\t".join(line)))




def check_coding(exons_file, CDSs_file, outfile, remove_overlapping = False):
        '''
        Given a bed file of exon coordinates and a bed file of CDS coordinates,
        writes a new bed file that only contains those exon coordinates form the former file that
        1) are fully coding
        2) are internal
        NB! Assumes that all the coordinates are from non-overlapping transcripts.
        If this is not the case, set remove_overlaps to True and it'll remove overlapping
        intervals.
        '''
        if remove_overlapping:
                bmo.sort_bed(exons_file, exons_file)
                remove_overlaps(exons_file, exons_file)
        #filter out anything that isn't fully coding
        #you have to write_both because you want to make sure that they
        #haven't been kept because of an overlap to a transcript that doesn't appear in the exons file
        temp_file = "temp_data/temp{0}.txt".format(random.random())
        bmo.intersect_bed(exons_file, CDSs_file, overlap = 1, overlap_rec = True, output_file = temp_file, force_strand = True, write_both = True, no_dups = False, no_name_check = False)
        #filter out terminal exons
        #in theory, there shouldn't be any left after the previous step
        #in practice, there may be unannotated UTRs, so it looks like we have a fully coding terminal exon,
        #whereas in reality, the exon is only partially coding
        temp_file2 = "temp_data/temp{0}.txt".format(random.random())
        with open(temp_file2, "w") as o_file:
                #figure out the rank of the last exon for each transcript
                filt_exons = gen.read_many_fields(exons_file, "\t")
                filt_exons = [i for i in filt_exons if len(i) > 3]
                names = [i[3].split(".") for i in filt_exons]
                names = gen.list_to_dict(names, 0, 1, as_list = True)
                names = {i: max([int(j) for j in names[i]]) for i in names}
                coding_exons = gen.read_many_fields(temp_file, "\t")
                for exon in coding_exons:
                        overlap_name = exon[9].split(".")
                        if overlap_name[0] in names:
                                name = exon[3].split(".")
                                if name[-1] != "1":
                                        last_exon = names[name[0]]
                                        if int(name[-1]) != last_exon:
                                                exon = [str(i) for i in exon[:6]]
                                                o_file.write("\t".join(exon))
                                                o_file.write("\n")
        bmo.sort_bed(temp_file2, temp_file2)
        gen.run_process(["mergeBed", "-i", temp_file2, "-c", "4,5,6", "-o", "distinct,distinct,distinct"], file_for_output = outfile)
        gen.remove_file(temp_file)
        gen.remove_file(temp_file2)

def check_sequence_quality(names, seqs, check_acgt=None, check_stop=None, check_start=None, check_length=None, check_inframe_stop=None, all_checks=None):
        '''
        Check sequence quality of fasta entries.
        Options:
        check_acgt: check if sequence only containts A,C,G,T nucleotides
        check_stop: check if the sequence has a recognised stop codon
        check_start: check if the sequence has a correctly defined ATG start codon
        check_length: check if the lenth of the sequence is a multiple of 3
        check_inframe_stop: check if the sequence contains an in frame stop codon
        '''

        stop_codons = ['TAA', 'TAG', 'TGA']
        passed_names = []
        passed_seqs = []

        print("Will filter sequences for quality. Starting out with {0} sequences.".format(len(seqs)))
        for i, seq in enumerate(seqs):
                actg_pass, stop_pass, start_pass, length_pass, inframe_stop_pass = True, True, True, True, True
                #check whether sequence only containts ACTG
                if check_acgt or all_checks:
                        actg_regex = re.compile('[^ACTG]')
                        non_actg = re.subn(actg_regex, '0', seq)[1]
                        if non_actg != 0:
                                actg_pass = False
                #check whether sequence contains a standard stop codon
                if check_stop or all_checks:
                        if seq[-3:] not in stop_codons:
                                stop_pass = False
                #check whether sequence containts correct start codon
                if check_start or all_checks:
                        if seq[:3] not in ['ATG']:
                                start_pass = False
                #check whether the sequence is a multiple of 3
                if check_length or all_checks:
                        if len(seq) % 3 != 0:
                                length_pass = False
                #check whether the sequence contains an in frame stop
                if check_inframe_stop or all_checks:
                        in_frame_stop_regex = re.compile('.{3}')
                        codons = re.findall(in_frame_stop_regex, seq[:-3])
                        if len(np.intersect1d(stop_codons, codons)) > 0:
                                inframe_stop_pass = False
                #if the sequence passes all checks, return
                if [actg_pass, stop_pass, start_pass, length_pass, inframe_stop_pass] == [True, True, True, True, True]:
                        passed_names.append(names[i])
                        passed_seqs.append(seq)
        print("After filtering, {0} sequences remain".format(len(passed_seqs)))
        return(passed_names, passed_seqs)

def extract_cds(gtf, bed_output, output_fasta, genome_fasta, full_chr_name=None, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None, all_checks=None, uniquify = False,clean_chrom_only=False):
    '''
    Given a .gtf file, extract the coding sequences to a fasta file.
    clean_chrom_only: do not include CDSs from contigs that haven't been assigned to a particular chromosome.
    EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf", "../output_data/Homo_sapiens.GRCh37.87_cds.fasta", "../source_data/Genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
    '''
    # extract required cds features to bed file
    bed_output_not_filtered = "{0}_not_filtered{1}".format(os.path.splitext(bed_output)[0], os.path.splitext(bed_output)[1])
    extract_features(gtf, bed_output_not_filtered, ['CDS', 'stop_codon'], full_chr_name, clean_chrom_only = clean_chrom_only)

    # extract features to fasta
    # args: bed file of features, fasta output file for cds, genome fasta, file that contains the parts that make up the cdss
    fasta_intervals_file_not_filtered = "{0}_intervals_not_filtered{1}".format(os.path.splitext(output_fasta)[0], os.path.splitext(output_fasta)[1])
    extract_cds_from_bed(bed_output_not_filtered, output_fasta, genome_fasta, fasta_intervals_file_not_filtered, check_acgt, check_start, check_length, check_stop, check_inframe_stop, all_checks, uniquify)

    # filter the previous bed file to only include those that passed the filters
    filter_bed_from_fasta(bed_output_not_filtered, output_fasta, bed_output)

    # filter the fasta intervals to only keep those that passed filtering
    fasta_intervals_file = "{0}_intervals{1}".format(os.path.splitext(output_fasta)[0], os.path.splitext(output_fasta)[1])
    filter_fasta_intervals_from_fasta(fasta_intervals_file_not_filtered, output_fasta, fasta_intervals_file)



def extract_cds_from_bed(bed_file, output_fasta, genome_fasta, fasta_interval_file, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None, all_checks=None, uniquify = False):
        '''
        Extract the CDS to fasta file
        Ex.: extract_cds('../feature_file.bed', '../output_file_fasta.fasta', '../source_data/genome_fasta_file.fa')
        '''
        #create dictionaries to hold cds parts
        cds_list = collections.defaultdict(lambda: collections.defaultdict())
        stop_list = {}
        concat_list = collections.defaultdict(lambda: collections.UserList())
        #create fasta file with extracted parts
        fasta_from_intervals(bed_file, fasta_interval_file, genome_fasta, names = True)
        #read the fasta interval file
        entries = gen.read_fasta(fasta_interval_file)
        # get the entry names and seqs
        sample_names = entries[0]
        seqs = entries[1]
        #set up the regex to get entry meta needed
        entry_regex = re.compile("(\w+)\.(\d+)(\..*)*")
        #iterate through the samples
        for i, sample in enumerate(sample_names):
                entry_meta = re.search(entry_regex, sample)
                #set the sample name: sample(.exon)
                sample_name = entry_meta.group(1)
                #if stop, set sample stop or send sample name to dict, with each part and its seq
                if seqs[i] in ['TAA', 'TAG', 'TGA']:
                        stop_list[sample_name] = seqs[i]
                else:
                        cds_list[sample_name][entry_meta.group(2)] = seqs[i]
        #get sorted list of seq parts
        for sample in sorted(cds_list):
                for part in sorted(cds_list[sample]):
                        concat_list[sample].append(cds_list[sample][part])
                #append the stop codon if it exists
                if sample in stop_list:
                        concat_list[sample].append(stop_list[sample])
        #concatenate and write to output
        names = []
        seqs = []
        for sample in sorted(concat_list):
                names.append(sample)
                seqs.append("".join(concat_list[sample]))
        #perform sequence quality control checks
        if check_acgt or check_stop or check_start or check_length or check_inframe_stop or all_checks:
                names, seqs = check_sequence_quality(names, seqs, check_acgt, check_stop, check_start, check_length, check_inframe_stop, all_checks)
        if uniquify:
                #leave only one transcript per gene
                gene_to_trans = link_genes_and_transcripts(bed_file)
                names, seqs = uniquify_trans(names, seqs, gene_to_trans)
                print("After leaving only one transcript per gene, {0} sequences remain.".format(len(seqs)))
        #write to output fasta file
        gen.write_to_fasta(names, seqs, output_fasta)

def extract_exons(gtf, bed):
        '''
        Given a GTF file, extract exon coordinates and write them to .bed.
        EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf",
        "../source_data/Homo_sapiens.GRCh37.87_exons.bed")
        '''
        #extract exons from GTF
        exons = gen.run_process(["grep", "\texon\t", gtf])
        #filter down to only protein-coding ones
        exons = gen.run_process(["grep", "transcript_biotype \"protein_coding\""], input_to_pipe = exons)
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

def extract_exon_junctions(exons, bed, window_of_interest=None):
        '''
        Given the file of extracted exons (generated using extract_exons), extract the coordinates of the junctions and write to .bed
        Set window_of_interest to a number of nucletides that you wish to examine across the junction
        EX.: extract_exon_junctions("../source_data/Homo_sapiens.GRCh37.87_exons.bed", "../source_data/Homo_sapiens.GRCh37.87_exon_junctions.bed", 30)
        '''

        #set up default dict to store info
        exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))
        #precompile regex to extract transcript id and exon id
        trans_exon_regex = re.compile(r"(?<=ENST)([0-9]*)\.([0-9]*)")
        #iterate over all exons and sort
        with open(exons, "r") as file:
                #read the line
                readLines = file.readlines()
                #for each exon
                for line in readLines:
                        #split line and return info
                        splits = line.strip('\n').split("\t")
                        chro = splits[0]
                        start = int(splits[1])
                        stop = int(splits[2])
                        trans = splits[3]
                        strand = splits[5]
                        #extract identifiers
                        trans = re.search(trans_exon_regex, trans)
                        trans_id = trans.group(1)
                        exon_id = int(trans.group(2))
                        #add to the dictionary
                        exon_list[chro][strand][trans_id][exon_id] = [start, stop]
        #open the output file
        out_file = open(bed, "w")

        #this is a bit clunky
        # for each chromosome, strand, transcript, exon, see if there is a 'next' exon
        # if there is write to file
        # I changed "chr" to "chrom" to avoid conflict with the in-built type "chr"
        for chrom in sorted(exon_list):
                for strand in sorted(exon_list[chrom]):
                        for trans_id in sorted(exon_list[chrom][strand]):
                                #create blank transcript output so we arent writing to file twice
                                for exon_id in sorted(exon_list[chrom][strand][trans_id]):
                                        if(exon_id+1 in exon_list[chrom][strand][trans_id]):

                                                #get exons for ease
                                                exon1 = exon_list[chrom][strand][trans_id][exon_id]
                                                exon2 = copy.deepcopy(exon_list[chrom][strand][trans_id][exon_id+1])

                                                #if window is defined, extract junction of size defined
                                                if window_of_interest:
                                                        #ensure window is even number
                                                        if window_of_interest % 2 != 0:
                                                                window_of_interest = window_of_interest +1
                                                        #get half the window interval
                                                        window_half = int(window_of_interest/2)

                                                        #if exon1 is bigger than the window interval, redefine window of interest
                                                        if(exon1[1] - exon1[0] > window_half):
                                                                exon1[0] = exon1[1]-window_half
                                                        #if exon2 is bigger than the window interval, redefine window of interest
                                                        if(exon2[1] - exon2[0] > window_half):
                                                                exon2[1] = exon2[0]+window_half

                                                if strand == "+":
                                                        exon1_site, exon2_site = 3,5
                                                elif strand == "-":
                                                        exon1_site, exon2_site = 5,3

                                                #write exon1 window to file
                                                out_file.write('{}\t{}\t{}\tENST{}.{}.{}\t.\t{}\n'.format(chrom,exon1[0],exon1[1],trans_id,exon_id,exon1_site,strand))
                                                #write exon2 window to file
                                                out_file.write('{}\t{}\t{}\tENST{}.{}.{}\t.\t{}\n'.format(chrom,exon2[0],exon2[1],trans_id,exon_id+1,exon2_site,strand))

        #close file
        out_file.close()

def extract_features(gtf_file, out_file, features, full_chr_name=None, clean_chrom_only = False):
        '''
        Given a GTF file, extract exon coordinates for specific features and write to .bed.
        EX.: extract_features("../source_data/Homo_sapiens.GRCh37.87.gtf", "../source_data/Homo_sapiens.GRCh37.87_exons.bed", ['CDS', 'stop_codon'])
        '''
        feature_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.UserList())))))
        if len(features) > 1:
                list_feature = True
        else:
                list_feature = False

        if clean_chrom_only:
                allowed = [str(i) for i in range(1, 23)] + ["X", "Y"]

        lines = gen.read_many_fields(gtf_file, "\t")
        #ensure protein coding and not entry meta
        lines = [line for line in lines if not line[0].startswith("#") and "pseudogene" not in line[-1]]
        #compile regex to find genes, transcripts and exons
        gene_regex = re.compile("(?<=gene_id \")ENSG[0-9]*")
        trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
        exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")

        for line in lines:
                #if the feature identifier is in the list
                if line[2] in features:
                        #get entry meta
                        meta = line[-1]
                        gene = re.search(gene_regex, meta).group(0)
                        trans = re.search(trans_regex, meta).group(0)
                        exon_no = re.search(exon_no_regex, meta).group(0)
                        chr_no = line[0]
                        add = True
                        if clean_chrom_only:
                                if chr_no not in allowed:
                                        add = False
                        if add:
                                feature_list[chr_no][gene][trans][exon_no][line[2]].append([line[3], line[4], line[6]])
        #output features sorted by feature, chr, gene, transcript id
        with open(out_file, 'w') as output:
                for chr_no in sorted(feature_list):
                        for gene in sorted(feature_list[chr_no]):
                                for trans in sorted(feature_list[chr_no][gene]):
                                        for exon in sorted(feature_list[chr_no][gene][trans]):
                                                for feature in sorted(feature_list[chr_no][gene][trans][exon]):
                                                        for item in feature_list[chr_no][gene][trans][exon][feature]:
                                                                if not list_feature:
                                                                        feature = '.'
                                                                if full_chr_name:
                                                                        chr_name = "chr{0}".format(chr_no)
                                                                else:
                                                                        chr_name = str(chr_no)
                                                                #output and convert to base 0
                                                                output.write('\t'.join([chr_name, str(int(item[0])-1), item[1], '{0}.{1}.{2}'.format(trans, exon, gene), feature, item[2]]) + '\n')


def extract_nt_indices(fasta_file, output_files):

    '''
    Extract the indices for each nt given a fasta file
    Output files need to be of format: output_files: "A": "filepath_for_A", "C", "filepath_for_C" etc
    '''

    nts = ["A", "C", "G", "T"]
    names, seqs = gen.read_fasta(fasta_file)
    # pos_regex = re.compile('^(chr\d+):(\d+)-(\d+)(?=\([+-]\))');

    indices = collections.defaultdict(lambda: collections.defaultdict(lambda: []))

    for i, seq in enumerate(seqs):
        id = names[i].strip('>')
        for nt in nts:
            indices[id][nt] = [str(m.start(0)) for m in re.finditer('{0}'.format(nt), seq)]

    outfiles = {}
    for nt in nts:
        outfiles[nt] = open(output_files[nt], "w")


    for id in indices:
        for nt in indices[id]:
            if len(indices[id][nt]) > 0:
                outfiles[nt].write(">{0}\n".format(id))
                outfiles[nt].write("{0}\n".format(",".join(indices[id][nt])))

    for nt in nts:
        outfiles[nt].close()

        # positions = re.search(pos_regex, id)
        # sampleid = positions.group(1)
        # start = int(positions.group(2))
        # nts = list(seq)
        # for j, nt in enumerate(nts):
        #     indicies[sampleid][nt].append(start+j)

    # with open(output_file, "w") as output:
    #     for sampleid in sorted(indicies):
    #         output.write('>{0}\n'.format(sampleid))
    #         line = ""
    #         for nt in sorted(indicies[sampleid]):
    #             # remove any duplicates
    #             indicies[sampleid][nt] = sorted(list(set(indicies[sampleid][nt])))
    #             indicies[sampleid][nt] = [str(x) for x in indicies[sampleid][nt]]
    #             line += "{0}:{1};".format(nt, ",".join(indicies[sampleid][nt]))
    #         output.write("{0}\n".format(line))


def fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = False):
    '''
    Takes a bed file and creates a fasta file with the corresponding sequences.
    If names == False, the fasta record names will be generated from the sequence coordinates.
    If names == True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    '''

    #if the index file exists, check whether the expected features are present
    genome_fasta_index = genome_fasta + '.fai'
    if(os.path.exists(genome_fasta_index)):
        bed_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(bed_file, "\t")])))
        index_chrs = sorted(list(set([entry[0] for entry in gen.read_many_fields(genome_fasta_index, "\t")])))
        if(not set(bed_chrs).issubset(set(index_chrs))):
            gen.remove_file(genome_fasta_index)

    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    print(" ".join(bedtools_args))
    gen.run_process(bedtools_args)
    names, seqs = gen.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    gen.write_to_fasta(names, seqs, fasta_file)

def fasta_from_intervals_temp_file(bed_file, output_fasta, genome_fasta, random_directory=None):
        '''
        Create a temporary file to hold the fasta extractions
        '''
        random_int = np.random.randint(9999999,size=2)
        if random_directory:
                temp_directory_path = './temp_files/temp_fasta_files_{0}'.format(random_int[0])
        else:
                temp_directory_path = './temp_files/temp_fasta_files'
        #create temp directory if doesnt already exist
        gen.create_directory('./temp_files/')
        #delete temp fasta directory and create new
        gen.create_strict_directory(temp_directory_path)
        #set the temporary fasta file path
        temp_fasta_file = '{0}/{1}_{2}{3}'.format(temp_directory_path, os.path.splitext(os.path.basename(output_fasta))[0], random_int[1], os.path.splitext(os.path.basename(output_fasta))[1])
        temp_fasta_file = output_fasta
        fasta_from_intervals(bed_file, temp_fasta_file, genome_fasta, force_strand = True, names = True)
        return(temp_fasta_file, temp_directory_path)

def filter_bed_from_fasta(bed, fasta, out_bed, families_file = None):
        '''
        Given a bed file and a fasta file, filter the bed file to only leave records where the 'name' field appears
        among the names in the fasta file. Write to out_bed.
        If a families_file is given, leave only one (randomly picked) member per family.
        '''
        #add feature in here that enables overwrite of current file
        output_exists = False
        if Path(out_bed).exists():
            output_exists = True
            temp_file_name = "{0}.{1}{2}".format(os.path.splitext(out_bed)[0], random.random(), os.path.splitext(out_bed)[1])
        else:
            temp_file_name = out_bed

        print(temp_file_name)


        fasta_names, fasta_seqs = gen.read_fasta(fasta)

        #read in family information and pick one transcript per family
        if families_file:
                families = gen.read_families(families_file)
                #make sure the families file doesn't contain transcripts that are not in the fasta
                for pos, family in enumerate(families):
                        families[pos] = [i for i in family if i in fasta_names]
                flat_families = gen.flatten(families)
                #first fish out singletons
                fasta_names_new = [i for i in fasta_names if i not in flat_families]
                for family in families:
                        if family:
                                family_lengths = [len(fasta_seqs[fasta_names.index(i)]) for i in family]
                                fasta_names_new.append(family[family_lengths.index(max(family_lengths))])
                fasta_names = fasta_names_new.copy()
        bed_data = gen.read_many_fields(bed, "\t")
        #remove empty lines
        bed_data = [i for i in bed_data if len(i) > 3]
        #I have to ask Liam cause I'm having a hard time understanding this regex
        id_regex = re.compile("^(\w+).*")
        with open(temp_file_name, "w") as file:
            for line in bed_data:
                idn = re.search(id_regex, line[3])
                if idn:
                    #filter bed data
                    if idn.group(1) in fasta_names:
                        file.write("\t".join(line))
                        file.write("\n")
        #remove old file, replace with new
        if(output_exists):
            os.remove(out_bed)
            shutil.move(temp_file_name, out_bed)

def filter_exon_junctions(junctions_file, exons_file, out_file):
        '''
        Given two bed files, one containing exon junction coordinates and one containing exon coordinates,
        filter the former to only leave intervals that either overlap exons in the latter or form part of an exon-exon
        junction with an exon that appears in the exon file.
        '''
        #read in exons file
        exons = gen.read_many_fields(exons_file, "\t")
        #only leave name field, parse, remove empty lines
        exons = [i[3].split(".") for i in exons if len(i) > 1]
        #exons[1:] because the first line is the header
        exons = gen.list_to_dict(exons[1:], 0, 1, as_list = True)
        #open output file
        with open(out_file, "w") as o_file:
                #loop over exon junctions file
                with open(junctions_file) as ej_file:
                        for line in ej_file:
                                #ignore empty lines
                                if len(line) > 1:
                                        name_field = line.split("\t")[3]
                                        name_field = name_field.split(".")
                                        #check if the transcript appears
                                        if name_field[0] in exons:
                                                #check if the exon appears
                                                if name_field[1] in exons[name_field[0]]:
                                                        o_file.write(line)
                                                #check if the upstream/downstream exon appears (depending on whether it's the 3' or 5' part of the junction)
                                                else:
                                                        if name_field[2] == "3":
                                                                if str(int(name_field[1]) + 1) in exons[name_field[0]]:
                                                                        o_file.write(line)
                                                        #elif safer than else
                                                        elif name_field[2] == "5":
                                                                if str(int(name_field[1]) - 1) in exons[name_field[0]]:
                                                                        o_file.write(line)

def filter_fasta_intervals_from_fasta(intervals_fasta, fasta, output):
        '''
        Given a fasta file and a fasta intervals file, filter the intervals file to only leave records where the 'name' field appears
        among the names in the fasta file. Write to fasta.
        '''
        #add feature in here that enables overwrite of current file
        output_exists = False
        if Path(output).exists():
            output_exists = True
            temp_file_name = "{0}.{1}{2}".format(os.path.splitext(output)[0], random.random(), os.path.splitext(output)[1])
        else:
            temp_file_name = output

        #fish out the names in the fasta
        fasta_names = gen.run_process(["grep", ">", fasta])
        fasta_names = fasta_names.split("\n")
        #remove tag and newline from each name
        fasta_names = [(i.lstrip("\>")).rstrip("\n") for i in fasta_names]
        #remove and potential blank entries
        fasta_names = [i for i in fasta_names if len(fasta_names) > 3]

        #read in the interval data
        fasta_interval_names, fasta_interval_seqs = gen.read_fasta(intervals_fasta)
        id_regex = re.compile("^(\w+).*")
        with open(temp_file_name, "w") as file:
            for i, interval in enumerate(fasta_interval_names):
                #search for the sample name
                id = re.search(id_regex, interval)
                if id:
                    trans_id = id.group(1)
                    #if the sample name is in the fasta names, output to file
                    if trans_id in fasta_names:
                        file.write(">{0}\n{1}\n".format(fasta_interval_names[i], fasta_interval_seqs[i]))
        #remove old file, replace with new
        if(output_exists):
            os.remove(output)
            shutil.move(temp_file_name, output)

def get_descriptions(names, gtf, out_file):
        '''
        Given a set of Ensembl transcript identifiers and a GTF file,
        determine the corresponding "gene name" for each transcript identifier.
        '''
        name_regex = re.compile("(?<=gene_name \")[A-z0-9\.\-\/\(\)]*(?=\")")
        trans_regex = re.compile("(?<=transcript_id \")[A-z0-9]*(?=\")")
        transcript_lines = gen.run_process(["grep", "\ttranscript\t", gtf])
        transcript_lines = transcript_lines.split("\n")
        with open(out_file, "w") as file:
                for line in transcript_lines:
                        if len(line) > 1:
                                trans = re.search(trans_regex, line).group(0)
                                if trans in names:
                                        description = re.search(name_regex, line).group(0)
                                        file.write("{0}\t{1}\n".format(trans, description))

def link_genes_and_transcripts(bed):
        '''
        Given a bed file of exons or similar genomic features, make a dictionary with gene IDs as keys and associated transcript IDs as values.
        '''
        bed_data = gen.read_many_fields(bed, "\t")
        bed_data = [i[3].split(".") for i in bed_data if len(i) > 2]
        bed_data = gen.list_to_dict(bed_data, 2, 0, as_list = True, uniquify = True)
        return(bed_data)

def remove_overlaps(in_bed, out_bed):
        '''
        Given a bed file, only leave non-overlapping elements, regardless of the strand of the overlap.
        '''
        #check how many columns there are in the bedfile
        with open(in_bed) as file:
                line = file.readline()
                column_number = line.count("\t") + 1
        #merge overlapping intervals and hav it count how many of the elements from the original file contribute to each
        #interval in the new file
        #note that bedops takes the column numbers in base 1
        if column_number > 3:
                columns = ",".join([str(i) for i in range(4, column_number + 1)] + ["1"])
                operations = ",".join(["distinct" for i in range(4, column_number + 1)] + ["count"])
        else:
                columns = "1"
                operations = "count"
        merge_result = gen.run_process(["bedtools", "merge", "-i", in_bed, "-c", columns, "-o", operations])
        #only leave those intervals that do not result from a merge and delete counts column
        merge_result = merge_result.split("\n")
        with open(out_bed, "w") as file:
                for line in merge_result:
                        if line[-2:] == "\t1":
                                line = line[:-2]
                                file.write(line)
                                file.write("\n")

def uniquify_trans(names, seqs, gene_to_trans):
        '''
        Filter a set of transcript IDs and corresponding sequences to only leave one transcript per gene (the longest).
        '''
        to_keep = []
        for gene in sorted(gene_to_trans):
                lengths = [len(seqs[names.index(i)]) for i in gene_to_trans[gene] if i in names]
                if lengths:
                	to_keep.append(gene_to_trans[gene][lengths.index(max(lengths))])
        seqs = [i for pos, i in enumerate(seqs) if names[pos] in to_keep]
        names = [i for i in names if i in to_keep]
        return(names, seqs)
