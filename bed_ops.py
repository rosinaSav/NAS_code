import generic as gen
from bam_ops import *
import re
import collections
import copy

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
    Given the file of extacted exons (generated using extract_exons), extact the coordinates of the junctions and write to .bed
    Set window_of_interest to a number of nucletides that you wish to examine across the junction
    EX.: extract_exon_junctions("../source_data/Homo_sapiens.GRCh37.87_exons.bed", "../source_data/Homo_sapiens.GRCh37.87_exon_junctions.bed", 30)
    '''

    #set up default dict to store info
    exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))
    #precompile regex to extact transcript id and exon id
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
    for chr in sorted(exon_list):
        for strand in sorted(exon_list[chr]):
            for trans_id in sorted(exon_list[chr][strand]):
                #create blank transcript output so we arent writing to file twice
                for exon_id in sorted(exon_list[chr][strand][trans_id]):
                    if(exon_id+1 in exon_list[chr][strand][trans_id]):

                        #get exons for ease
                        exon1 = exon_list[chr][strand][trans_id][exon_id]
                        exon2 = copy.deepcopy(exon_list[chr][strand][trans_id][exon_id+1])

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
                            #if exon1 is bigger than the window interval, redefine window of interest
                            if(exon2[1] - exon1[0] > window_half):
                                exon2[1] = exon2[0]+window_half

                        if strand == "+":
                            exon1_site, exon2_site = 3,5
                        elif strand == "-":
                            exon1_site, exon2_site = 5,3

                        #write exon1 window to file
                        out_file.write('{}\t{}\t{}\tENST{}.{}.{}\t.\t{}\n'.format(chr,exon1[0],exon1[1],trans_id,exon_id,exon1_site,strand))
                        #write exon2 window to file
                        out_file.write('{}\t{}\t{}\tENST{}.{}.{}\t.\t{}\n'.format(chr,exon2[0],exon2[1],trans_id,exon_id+1,exon2_site,strand))

    #close file
    out_file.close()

def fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = False):
    '''
    Takes a bed file and creates a fasta file with the corresponding sequences.
    If names == False, the fasta record names will be generated from the sequence coordinates.
    If names == True, the fasta name will correspond to whatever is in the 'name' field of the bed file
    '''
    bedtools_args = ["bedtools", "getfasta", "-s", "-fi", genome_fasta, "-bed", bed_file, "-fo", fasta_file]
    if not force_strand:
        del bedtools_args[2]
    if names:
        bedtools_args.append("-name")
    gen.run_process(bedtools_args)
    names, seqs = gen.read_fasta(fasta_file)
    seqs = [i.upper() for i in seqs]
    gen.write_to_fasta(names, seqs, fasta_file)


def extract_features(bed_file, out_file, features, list_feature=None, full_chromosome=False):

	feature_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.UserList())))

	lines = gen.read_many_fields(bed_file, "\t")
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
			#if you want the full chromosome name
			if full_chromosome:
				chr_regex = re.compile("[XY0-9]*")
				chr_no = re.search(chr_regex, chr_no).group(0)
			#if you want to list the feature extracted
			if list_feature:
				feature = line[2]
			else:
				feature = '.'
			feature_list[line[2]][chr_no][trans].append('{0}\t{1}\t{2}\t{3}\t{4}\t{5}.{6}\n'.format(chr_no,line[3],line[4],feature,line[6],trans,exon_no))

	#output features sorted by feature, chr, transcript id
	with open(out_file, 'w') as output:
		for feature in sorted(feature_list):
			for chr_no in sorted(feature_list[feature]):
				for trans in sorted(feature_list[feature][chr_no]):
					print(feature, chr_no, trans)
					for item in feature_list[feature][trans][chr_no]:
						print(item)
						output.write(item)
