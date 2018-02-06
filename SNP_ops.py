import generic as gen
import numpy as np
import os
import random
import string
import collections
import re

def get_relative_SNPs(CDSs, SNP_file_name, CDS_SNP_file_name, seqs, names, genome, get_new_SNPs = False, parse_SNPs = False, remove_GT = False):
    if get_new_SNPs:
        flat_CDS = [[j[0] for j in i] for i in list(CDSs.values())]
        CDS_bed = "temp_data/temp_{0}.bed".format(random.random())
        write_features_to_bed(flat_CDS, CDS_bed, modify_chr_ids = True)
        tabix(CDS_bed, SNP_file_name, genome)
        os.remove(CDS_bed)
    if parse_SNPs:
        relative_coords_dict = {i: {"positions": [], "to_remove": [], "var_bases": [], "MAFs": [], "allele_counts": []} for i in list(CDSs.keys())}
        counter = 0
        with open(SNP_file_name) as file:
            error_counter = 0
            #loop over bed file with SNPs
            for line in file:
                if line[0] != "#":
                    if counter % 1000 == 0:
                        print(counter)
                    counter = counter + 1
                    bed_record = line.split("\t")
                    info = bed_record[6].split("$")
                    #reported ancestral base
                    ref_base = info[1]
                    #filter out anything that isn't simple SNPs
                    var_base = info[2].split(",")
                    var_base_number = len(var_base)
                    var_base = [i for i in var_base if i in nc._canon_bases_]
                    trans = bed_record[3]
                    #convert bed record to indices relative to the CDS
                    #mismatch_warning_only = True means that if the SNP coordinates don't map onto the CDS,
                    #return an error message but don't crash.
                    current_index = bed_to_CDS_indices(bed_record, CDSs[trans], mismatch_warning_only = True)[0]
                    #this is so you could remove the stuff you filter out for various reasons also from the monomorphic set
                    relative_coords_dict[trans]["to_remove"].append(str(current_index))
                    #the SNP bed files contain information on which transcript they overlap
                    if current_index != ("error"):
                        if len(ref_base) == 1:
                            #filtering out polymorphisms with more than 2 segregating alleles
                            #note that you need to check the number both before and after filtering out non-canonical bases
                            if var_base_number == 1 and len(var_base) == 1:
                            #check if the reported ancestral allele matches the base that is at that position in the CDS
                            #it's expected that sometimes it won't but it should most of the time if everything's worked out properly
                                strand = CDSs[trans][0][0][6]
                                var_base = var_base[0]
                                if strand == "-":
                                    ref_base = nc.rev_comp(ref_base)
                                    var_base = nc.rev_comp(var_base)
                                CDS_base = seqs[names.index(trans)][current_index]
                                if ref_base != CDS_base:
                                    print("PROBLEM!")
                                    print(ref_base)
                                    print(CDS_base)
                                    print(var_base)
                                    print("\n")
                                elif remove_GT and ref_base == "G" and var_base == "T":
                                    pass
                                else:
                                    relative_coords_dict[trans]["positions"].append(str(current_index))
                                    relative_coords_dict[trans]["var_bases"].append(str(var_base))
                                    MAF_MAC = info[4].split(";")
                                    MAC = MAF_MAC[0].split("=")[1]
                                    MAF = MAF_MAC[1].split("=")[1]
                                    relative_coords_dict[trans]["MAFs"].append(str(MAF))
                                    relative_coords_dict[trans]["allele_counts"].append(str(MAC))
                    else:
                        error_counter = error_counter + 1
        with open(CDS_SNP_file_name, "w") as file:
            for trans in relative_coords_dict:
                to_write = zip(relative_coords_dict[trans]["positions"], relative_coords_dict[trans]["var_bases"], relative_coords_dict[trans]["MAFs"], relative_coords_dict[trans]["allele_counts"])
                to_write = "\t".join([trans, "|".join([",".join(i) for i in to_write]), ",".join(relative_coords_dict[trans]["to_remove"])])
                file.write(to_write)
                file.write("\n")
        relative_coords_dict = {}
        print("Number of conversion errors:")
        print(error_counter)

    current_SNPs = rw.read_many_fields(CDS_SNP_file_name, "\t")
    current_SNPs_to_remove = list_to_dict(current_SNPs, 0, 2)
    current_SNPs_to_remove = {i: [int(j) for j in current_SNPs_to_remove[i].split(",") if j != "error"] for i in current_SNPs_to_remove if current_SNPs_to_remove[i]}
    current_SNPs = list_to_dict(current_SNPs, 0, 1)
    current_SNPs = {i: [j.split(",") for j in current_SNPs[i].split("|")] for i in current_SNPs if current_SNPs[i]}
    current_SNPs = {i: {int(j[0]): (j[1], float(j[2])) for j in current_SNPs[i]} for i in current_SNPs}
    return(current_SNPs, current_SNPs_to_remove)

def get_snp_feature_relative_position(intersect_file, out_file):
	'''
	Get the relative position of a snp within the feature it is found. Used as an intermediate step
	before calculating the snp position in the cds using get_snp_cds_relative_position.
	'''
	#read the intersects file
	intersects = gen.read_many_fields(intersect_file, "\t")
	with open(out_file, "w") as output:
		for intersect in intersects:
			#get the features
			feature_start = intersect[1]
			feature = intersect[3]
			snp_start = intersect[7]
			#get the meta about the snp
			snp_meta = intersect[6:]
			#get the position of the snp compared with the feature
			relative_position = int(snp_start) - int(feature_start)
			#replace . field with relative position
			snp_meta[5] = str(relative_position)
			output.write("{0}\n".format("\t".join(snp_meta)))

def get_snp_cds_relative_position(features_relative_position, fasta_parts, out_file):
	'''
	Get the position of the snp within a CDS using the relative positions of snps in the features they are found
	'''
	#read snp positions relative to their feature
	snp_relative_positions = gen.read_many_fields(features_relative_position, "\t")
	#get the fasta entries of the features
	feature_names, feature_seqs =  gen.read_fasta(fasta_parts)
	#set up dict to hold the feature positions relative to the cds
	cds_features_relative_positions = collections.defaultdict(lambda: collections.defaultdict())

	entry_regex = re.compile("(\w+)\.(\d+)(\..*)*")
	for i, feature_name in enumerate(feature_names):
		entry_meta = re.search(entry_regex, feature_name)
		#get the lengths of each feature and those before it, and sum
		rel_pos = sum([len(feature_seqs[x]) for x in range(0,i)])
		cds_features_relative_positions[entry_meta.group(1)][entry_meta.group(2)] = rel_pos

	with open(out_file, "w") as output:
		for snp_relative_position in snp_relative_positions:
			id = snp_relative_position[3]
			id_meta = re.search(entry_regex, id)
			snp_position = int(snp_relative_position[5])
			#get the position of the snp relative to the cds, by summing the feature start position and the snp position
			snp_cds_relative_position = cds_features_relative_positions[id_meta.group(1)][id_meta.group(2)] + snp_position
			snp_relative_position[5] = str(snp_cds_relative_position)
			output.write("{0}\n".format("\t".join(snp_relative_position)))



def tabix(bed_file, output_file, vcf, process_number = None):
    '''
    Given a bed file, use tabix to get overlapping 1000Genomes SNPs.
    bed_file: input bed file
    output_file: output file name
    vcf: vcf summary file
    process_number: number of parallel processes
    '''
    #divide the input bed file into smaller files so that you could parallelize
    #if the number of processes has not been specified, use CPU count/2.
    if not process_number:
        process_number = int(os.cpu_count()/2)
    bed_file_length = gen.line_count(bed_file)
    #if the input bed_file has fewer lines than process_number, use 2 processes
    if bed_file_length <= process_number:
        process_number = 2
    lines_per_file = int(bed_file_length/process_number)
    #split the files
    gen.run_process(["split", "-l", lines_per_file, bed_file, bed_file])
    #split automatically names the files it creates according to the following pattern
    #make all the possible names and then only keep the ones you expect to have given the number of processes
    bed_names = ["{0}a{1}".format(bed_file, i) for i in string.ascii_lowercase]
    if (bed_file_length%process_number) == 0:
        bed_names = bed_names[:process_number]
    else:
        bed_names = bed_names[:(process_number + 1)]
    #run the core tabix function on the bed files
    parallel_tabix = gen.run_in_parallel(bed_names, ["foo", vcf], tabix_core, workers = process_number)
    [i.get() for i in parallel_tabix]
    #this is how the output files have been anmed by the core tabix function
    output_files = ["{0}.out".format(i) for i in bed_names]
    #concatenate the output files
    gen.run_process(["cat {0}??.out".format(bed_file)], file_for_output = output_file, shell = True)
    #clean up temp files
    [os.remove(i) for i in bed_names]
    [os.remove(i) for i in output_files]

def tabix_core(bed_files, vcf):
    '''
    The code that's parallelized in tabix above.
    '''
    #loop over the input bed files
    for curr_bed_file in bed_files:
        curr_output_file = curr_bed_file + ".out"
        with open(curr_bed_file) as file, open(curr_output_file, "w") as file2:
            counter = 0
            #loop over the lines in the current bed file
            for line in file:
                #print out every 100th line number
                counter = gen.update_counter(counter, 100)
                #parse current line
                line = line.split("\t")
                chrom = line[0].lstrip("chr")
                #you need to add 1 because bed files are 0-based, whereas the vcf files are 1-based
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                #get the SNPs within the current bed file coords
                output = gen.run_process(["tabix", vcf, "{0}:{1}-{2}".format(chrom, start, end)])
                if output:
                    #parse output
                    output = output.rstrip("\n")
                    output = output.split("\n")
                    output = [i.split("\t") for i in output]
                    output = [i for i in output if i[2][:2] == "rs"]
                    if output:
                        #format output in .bed format (note the passage to base 0) and write to file
                        output = ["\t".join(["chr{0}".format(i[0]), str(int(i[1]) - 1), i[1], trans, "100", ".", "$".join([i[2], i[3], i[4], i[6], i[7]])]) for i in output]
                        output = "\n".join(output)
                        file2.write(output)
                        file2.write("\n")

def tabix_samples(bed_file, output_file_name, panel_file, vcf_folder, superpop = None, subpop = None, samples = None, downsample_by = None, exclude_xy = False, vcftools_path = None):
    '''
    Extract 1000Genomes SNPs for a subpopulation.
    bed_file: input bed file for the intervals you want
    output_file_name: name of output file
    panel_file: path to panel file
    vcf_folder: directory that contains the per-individual VCF files

    superpop: population code if you want to filter by population. Possible:
    AFR, African
    AMR, Ad Mixed American
    EAS, East Asian
    EUR, European
    SAS, South Asian

    subpop: subpopulation code if you want to filter by subpopulation. Possible:
    CHB 	Han Chinese in Beijing, China
    JPT 	Japanese in Tokyo, Japan
    CHS 	Southern Han Chinese
    CDX 	Chinese Dai in Xishuangbanna, China
    KHV 	Kinh in Ho Chi Minh City, Vietnam
    CEU 	Utah Residents (CEPH) with Northern and Western European Ancestry
    TSI 	Toscani in Italia
    FIN 	Finnish in Finland
    GBR 	British in England and Scotland
    IBS 	Iberian Population in Spain
    YRI 	Yoruba in Ibadan, Nigeria
    LWK 	Luhya in Webuye, Kenya
    GWD 	Gambian in Western Divisions in the Gambia
    MSL 	Mende in Sierra Leone
    ESN 	Esan in Nigeria
    ASW 	Americans of African Ancestry in SW USA
    ACB 	African Caribbeans in Barbados
    MXL 	Mexican Ancestry from Los Angeles USA
    PUR 	Puerto Ricans from Puerto Rico
    CLM 	Colombians from Medellin, Colombia
    PEL 	Peruvians from Lima, Peru
    GIH 	Gujarati Indian from Houston, Texas
    PJL 	Punjabi from Lahore, Pakistan
    BEB 	Bengali from Bangladesh
    STU 	Sri Lankan Tamil from the UK
    ITU 	Indian Telugu from the UK

    samples: list of sample IDs, if you wish to filter in that way. Ex: ["NA21141", "NA21142", "NA21143"]. Also possible to supply a single one: ["NA21141"].

    downsample_by: if you want to randomly pick only a subsample of the individuals. For example, 0.5 would give you half of the individuals.
    exclude_xy: if True, SNPs on sex chromosomes will not be returned
    vcftools_path: path from current location to vcf-tools. Example: "../Software/vcftools"
    '''

    sex_chromosomes = ["Y", "X"]

    if not samples:

        #read and parse panel file
        panel = gen.read_many_fields(panel_file, "\t")
        panel = [i for i in panel if len(i) == 4]

        #filter by (sub)population. The first element of each line is the sample ID so this bit just makes a list of sample IDs.
        if subpop:
            samples = [i[0] for i in panel if i[1] == subpop]
        elif superpop:
            samples = [i[0] for i in panel if i[2] == superpop]
        else:
            samples = [i[0] for i in panel]

    #downsample if needed
    if downsample_by:
        samples = np.random.choice(samples, size = int(len(samples)/downsample_by), replace = False)

    print(len(samples))

    #turn list into comma-separated string
    samples = ",".join(samples)

    if not vcftools_path:
        vcftools_path = ""

    #get a list of all the files in the VCF folder
    vcf_files = os.listdir(vcf_folder)
    #just in case there is nonsense in the directory
    vcf_files = [i for i in vcf_files if "vcf" in i]

    with open(bed_file) as file:
        sample_files = []
        counter = 0
        #loop over lines in bed file
        for line in file:
            #print out every 100th line number
            counter = gen.update_counter(counter, 100)
            #parse line in bed file
            line = line.split("\t")
            chrom = line[0].lstrip("chr")
            if chrom in sex_chromosomes and exclude_xy:
                print("Only autosomes are processed!")
            else:
                #add 1 to start coordinate because bed files are 0-based, whereas the vcf files are 1-based
                start = int(line[1]) + 1
                end = line[2]
                trans = line[3]
                #get the vcf file for the right chromosome, making sure not to get the index file instead
                current_vcf = ["{0}/{1}".format(vcf_folder, i) for i in vcf_files if "chr{0}".format(chrom) in i and ".tbi" not in i]
                #check that you only got a single file
                if len(current_vcf) != 1:
                    print("Ambiguous or missing files in VCF folder!")
                    print(current_vcf)
                    print(chrom)
                    raise Exception
                else:
                    current_vcf = current_vcf[0]
                #generate temporary output file for all SNPs in interval
                temp_output_file = "temp_data/temp_vcf{0}.vcf".format(random.random())
                #get ALL SNPs (that is to say, for all samples) for current interval
                gen.run_process(["tabix", "-h", current_vcf, "{0}:{1}-{2}".format(chrom, start, end)], file_for_output = temp_output_file)
                #uncomment the following line for debug
##                run_process(["cp", temp_output_file, "temp_data/{0}:{1}-{2}_tabix_slice.txt".format(chrom, start, end)])
                #generate temporary output file for SNPs from your seleceted samples
                sample_output_file = "temp_data/temp_sample_tabix{0}.txt".format(random.random())
                sample_files.append(sample_output_file)
                #filter the file you made with all the SNPs to only leave the SNPs that appear in your samples
                gen.run_process(["{0}vcf-subset".format(vcftools_path), "-c", samples, temp_output_file], file_for_output = sample_output_file)
                gen.remove_file(temp_output_file)

    #you want to concatenate the sample files you made (one file per bed interval) but you can't in one go cause there's too many
    #therefore, you take the 10 last files, concatenate those
    #then concatenate the next 10 files (moving from the end of the list towards the beginning) to each-other and to the file you got in the previous step
    #etc.
    #you juggle the two temp concat file names just so you would be overwriting files rather than creating new ones
    concat_files = ["temp_data/temp_concat_file{0}.vcf".format(random.random()), "temp_data/temp_concat_file{0}.vcf".format(random.random())]
    current_sample_files = sample_files[-10:]
    del sample_files[-10:]
    gen.run_process(["{0}vcf-concat".format(vcftools_path)] + current_sample_files, file_for_output = concat_files[0])
    local_counter = 0
    files_left = True
    while files_left:
        local_counter = local_counter + 1
        current_sample_files = sample_files[-10:]
        del sample_files[-10:]
        if len(sample_files) == 0:
            files_left = False
        if local_counter%2 == 0:
            current_concat_file = concat_files[0]
            previous_concat_file = concat_files[1]
        else:
            current_concat_file = concat_files[1]
            previous_concat_file = concat_files[0]
        gen.run_process(["{0}vcf-concat".format(vcftools_path)] + current_sample_files + [previous_concat_file], file_for_output = current_concat_file)
    sort_file = "temp_data/temp_sort_file{0}.vcf".format(random.random())
    #once everything is concatenated, sort the SNPs, make a compressed version of the file and make an index for tabix
    gen.run_process(["vcf-sort".format(vcftools_path), current_concat_file], file_for_output = sort_file)
    gen.run_process(["bgzip", "-c", sort_file], file_for_output = output_file_name)
    gen.run_process(["tabix", "-f", "-p", "vcf", output_file_name])
    #clean up
    for sample_file in sample_files:
        gen.remove_file(sample_file)
    for concat_file in concat_files:
        gen.remove_file(concat_file)
    gen.remove_file(sort_file)
