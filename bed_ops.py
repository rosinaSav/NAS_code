import generic as gen
import bam_ops as bmo
import re
import collections
import copy
import numpy as np
import os
import shutil
from pathlib import Path
import random

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

	return(passed_names, passed_seqs)

def extract_cds(gtf, bed_output, output_fasta, genome_fasta, full_chr_name=None, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None, all_checks=None):
    '''
    Given a .gtf file, exrtract the coding sequences to a fasta file.
    EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf", "../output_data/Homo_sapiens.GRCh37.87_cds.fasta", "../source_data/Genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
    '''
    #extract required cds features
    extract_features(gtf, bed_output, ['CDS', 'stop_codon'], full_chr_name)
    #extract to fasta
    extract_cds_from_bed(bed_output, output_fasta, genome_fasta, check_acgt, check_start, check_length, check_stop, check_inframe_stop, all_checks)
    #filter the previous files to only include those that passed the filters

    print(bed_output)
    print(output_fasta)
    filter_bed_from_fasta(bed_output, output_fasta, bed_output)
    # fasta_interval_file = "{0}_intervals{1}".format(os.path.splitext(output_fasta)[0], os.path.splitext(output_fasta)[1])
    # filter_fasta_intervals_from_fasta(fasta_interval_file, output_fasta, fasta_interval_file)

def extract_cds_from_bed(bed_file, output_fasta, genome_fasta, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None, all_checks=None):
	'''
	Extract the CDS to fasta file
	Ex.: extract_cds('../feature_file.bed', '../output_file_fasta.fasta', '../source_data/genome_fasta_file.fa')
	'''
	#create dictionaries to hold cds parts
	cds_list = collections.defaultdict(lambda: collections.defaultdict())
	stop_list = {}
	concat_list = collections.defaultdict(lambda: collections.UserList())
	#create fasta file with intervals
	fasta_interval_file = "{0}_intervals{1}".format(os.path.splitext(output_fasta)[0], os.path.splitext(output_fasta)[1])
	fasta_from_intervals(bed_file, fasta_interval_file, genome_fasta, names=True)
	#read the interval fasta file
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
	#write to output fasta file
	gen.write_to_fasta(names, seqs, output_fasta)

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

def extract_features(gtf_file, out_file, features, full_chr_name=None):
	'''
	Given a GTF file, extract exon coordinates for specific features and write to .bed.
	EX.: extract_fetures("../source_data/Homo_sapiens.GRCh37.87.gtf", "../source_data/Homo_sapiens.GRCh37.87_exons.bed", ['CDS', 'stop_codon'])
	'''
	feature_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.UserList()))))
	if len(features) > 1:
		list_feature = True
	else:
		list_feature = False

	lines = gen.read_many_fields(gtf_file, "\t")
	#ensure protein coding and not entry meta
	lines = [line for line in lines if not line[0].startswith("#") and "pseudogene" not in line[-1]]
	#compile regex to fine transcripts and exons
	trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
	exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")

	for line in lines:
		#if the feature identifier is in the list
		if line[2] in features:
			#get entry meta
			meta = line[-1]
			trans = re.search(trans_regex, meta).group(0)
			exon_no = re.search(exon_no_regex, meta).group(0)
			chr_no = line[0]
			feature_list[chr_no][trans][exon_no][line[2]].append([line[3], line[4], line[6]])
	#output features sorted by feature, chr, transcript id
	with open(out_file, 'w') as output:
		for chr_no in sorted(feature_list):
			for trans in sorted(feature_list[chr_no]):
				for exon in sorted(feature_list[chr_no][trans]):
					for feature in sorted(feature_list[chr_no][trans][exon]):
						for item in feature_list[chr_no][trans][exon][feature]:
							if not list_feature:
								feature = '.'
							if full_chr_name:
								chr_name = "chr{0}".format(chr_no)
							else:
								chr_name = str(chr_no)
							#output and convert to base 0
							output.write('\t'.join([chr_name, str(int(item[0])-1), item[1], '{0}.{1}'.format(trans, exon), feature, item[2]]) + '\n')

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
	gen.run_process(bedtools_args)
	names, seqs = gen.read_fasta(fasta_file)
	seqs = [i.upper() for i in seqs]
	gen.write_to_fasta(names, seqs, fasta_file)

def filter_bed_from_fasta(bed, fasta, out_bed):
        '''
        Given a bed file and a fasta file, filter the bed file to only leave records where the 'name' field appears
        among the names in the fasta file. Write to out_bed.
        '''
        #add feature in here that enables overwrite of current file
        output_exists = False
        if Path(out_bed).exists():
            output_exists = True
            temp_file_name = "{0}.{1}{2}".format(os.path.splitext(out_bed)[0], random.random(), os.path.splitext(out_bed)[1])
        else:
            temp_file_name = out_bed

        print(temp_file_name)

        #fish out the names in the fasta
        fasta_names = gen.run_process(["grep", ">", fasta])
        fasta_names = fasta_names.split("\n")
        #remove tag and newline from each name
        fasta_names = [(i.lstrip("\>")).rstrip("\n") for i in fasta_names]
        bed_data = gen.read_many_fields(bed, "\t")
        #remove empty lines
        bed_data = [i for i in bed_data if len(i) > 3]
        id_regex = re.compile("^(\w+).*")
        with open(temp_file_name, "w") as file:
            for line in bed_data:
                id = re.search(id_regex, line[3])
                if id:
                    #filter bed data
                    if id.group(1) in fasta_names:
                        file.write("\t".join(line))
                        file.write("\n")
        #remove old file, replace with new
        if(output_exists):
            os.remove(out_bed)
            shutil.move(temp_file_name, out_bed)

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
        exons = gen.list_to_dict(exons, 0, 1, as_list = True)
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
