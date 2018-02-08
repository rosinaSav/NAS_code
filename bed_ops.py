import generic as gen
import bam_ops as bmo
import re
import collections
import copy
import numpy as np
import os
import random
import shutil

def check_coding(exons_file, CDSs_file, outfile):
        '''
        Given a bed file of exon coordinates and a bed file of CDS coordinates,
        writes a new bed file that only contains those exon coordinates form the former file that
        1) are fully coding
        2) are internal
        NB! Assumes that all the cooridnates are from non-overlapping transcripts.
        '''
        #filter out anything that isn't fully coding
        temp_file = "temp_data/temp{0}.txt".format(random.random())
        bmo.intersect_bed(exons_file, CDSs_file, overlap = 1, output_file = temp_file, force_strand = True, no_dups = True)
        #filter out terminal exons
        #in theory, there shouldn't be any left after the previous step
        #in practice, there may be unannotated UTRs, so it looks like we have a fully coding terminal exon,
        #whereas in reality, the exon is only partially coding
        with open(outfile, "w") as o_file:
                filt_exons = gen.read_many_fields(temp_file, "\t")
                filt_exons = [i for i in filt_exons if len(i) > 3]
                names = [i[3].split(".") for i in filt_exons]
                names = gen.list_to_dict(names, 0, 1, as_list = True)
                names = {i: [int(j) for j in names[i]] for i in names}
                for exon in filt_exons:
                        name = exon[3].split(".")
                        if name[-1] != "1":
                                last_exon = max(names[name[0]])
                                if int(name[-1]) != last_exon:
                                        exon = [str(i) for i in exon]
                                        o_file.write("\t".join(exon))
                                        o_file.write("\n")
        gen.remove_file(temp_file)

def check_sequence_quality(names, seqs, check_acgt=None, check_stop=None, check_start=None, check_length=None, check_inframe_stop=None):
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
		if check_acgt:
			actg_regex = re.compile('[^ACTG]')
			non_actg = re.subn(actg_regex, '0', seq)[1]
			if non_actg != 0:
				actg_pass = False
		#check whether sequence contains a standard stop codon
		if check_stop:
			if seq[-3:] not in stop_codons:
				stop_pass = False
		#check whether sequence containts correct start codon
		if check_start:
			if seq[:3] not in ['ATG']:
				start_pass = False
		#check whether the sequence is a multiple of 3
		if check_length:
			if len(seq) % 3 != 0:
				length_pass = False
		#check whether the sequence contains an in frame stop
		if check_inframe_stop:
			in_frame_stop_regex = re.compile('.{3}')
			codons = re.findall(in_frame_stop_regex, seq[:-3])
			if len(np.intersect1d(stop_codons, codons)) > 0:
				inframe_stop_pass = False
		#if the sequence passes all checks, return
		if [actg_pass, stop_pass, start_pass, length_pass, inframe_stop_pass] == [True, True, True, True, True]:
			passed_names.append(names[i])
			passed_seqs.append(seq)

	return(passed_names, passed_seqs)

def extract_cds(gtf, output_fasta, genome_fasta, random_directory=None, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None):
	'''
	Given a .gtf file, exrtract the coding sequences to a fasta file.
	EX.: extract_cds("../source_data/Homo_sapiens.GRCh37.87.gtf", "../output_data/Homo_sapiens.GRCh37.87_cds.fasta", "../source_data/Genomes/Homo_sapiens.GRCh37.dna.primary_assembly.fa")
	Use random_directory for creating a randomised directory to hold intermediate fasta components
	'''
	#create bed file path
	bed = os.path.splitext(gtf)[0] + ".bed"
	#extract required cds features
	extract_features(gtf, bed, ['CDS', 'stop_codon'])
	#extract to fasta
	extract_cds_from_bed(bed, output_fasta, genome_fasta, random_directory, check_acgt, check_start, check_length, check_stop, check_inframe_stop)

def extract_cds_from_bed(bed_file, output_fasta, genome_fasta, random_directory=None, check_acgt=None, check_start=None, check_length=None, check_stop=None, check_inframe_stop=None):
	'''
	Extract the CDS to fasta file
	Ex.: extract_cds('../feature_file.bed', '../output_file_fasta.fasta', '../source_data/genome_fasta_file.fa')
	'''
	#create dictionaries to hold cds parts
	cds_list = collections.defaultdict(lambda: collections.defaultdict())
	stop_list = {}
	concat_list = collections.defaultdict(lambda: collections.UserList())
	#create temp fasta file with extracted parts
	temp_fasta_file, temp_directory_path = fasta_from_intervals_temp_file(bed_file, output_fasta, genome_fasta, random_directory)
	#read the temp fasta file
	entries = gen.read_fasta(temp_fasta_file)
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
	if check_acgt or check_stop or check_start or check_length or check_inframe_stop:
		names, seqs = check_sequence_quality(names, seqs, check_acgt, check_stop, check_start, check_length, check_inframe_stop)
	#write to output fasta file
	gen.write_to_fasta(names, seqs, output_fasta)
	#remove the temporary directory
	shutil.rmtree(temp_directory_path)

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

def extract_features(gtf_file, out_file, features):
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
							#output and convert to base 0
							output.write('\t'.join([chr_no, str(int(item[0])-1), item[1], '{0}.{1}'.format(trans, exon), feature, item[2]]) + '\n')

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

def filter_bed_from_fasta(bed, fasta, out_bed):
        '''
        Given a bed file and a fasta file, filter the bed file to only leave records where the 'name' field appears
        among the names in the fasta file. Write to out_bed.
        '''
        #fish out the names in the fasta
        fasta_names = gen.run_process(["grep", ">", fasta])
        fasta_names = fasta_names.split("\n")
        #remove tag and newline from each name
        fasta_names = [(i.lstrip("\>")).rstrip("\n") for i in fasta_names]
        bed_data = gen.read_many_fields(bed, "\t")
        #remove empty lines
        bed_data = [i for i in bed_data if len(i) > 3]
        #filter bed data
        bed_data = [i for i in bed_data if i[3] in fasta_names]
        with open(out_bed, "w") as file:
                for line in bed_data:
                        file.write("\t".join(line))
                        file.write("\n")

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

                
