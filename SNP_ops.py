'''
Authors: Rosina Savisaar and Liam Abrahams
Module containing functions for processing SNP data and SNP-related operations.
'''

import bam_ops as bmo
from cython_func import calc_density_for_concat_several_c
import generic as gen
import numpy as np
import os
import random
import string
import collections
import re
import itertools

def filter_by_snp_type(input_file, output_file, snp_type, set_seed=None):
    '''
    Filter a file of processed SNP reads by SNP type.
    snp_type: ptc, syn, non
    '''
    #get header and anything that contains snp type
    grep_args = "echrom\|{0}".format(snp_type)
    gen.run_process(["grep", grep_args, input_file], file_for_output = output_file)

def filter_motif_SNPs(fasta, SNPs, motifs, out_file, complement = False):
    '''
    Filter a SNPs file to only leave those SNPs that overlap a set of motifs.
    If complement, only leave SNPs that DON'T overlap the motifs.
    '''
    #read in data
    with open(motifs) as file:
        motifs = file.readlines()[1:]
        motifs = [i.rstrip("\n") for i in motifs]
    motif_lengths = [len(i) for i in motifs]
    motifs = gen.motif_to_regex(motifs)
    names, seqs = gen.read_fasta(fasta)

    with open(SNPs) as file, open(out_file, "w") as ofile:
        header = file.readline()
        ofile.write(header)
        for line_raw in file:
            line = line_raw.split("\t")
            ORF_pos = int(line[11])
            trans_name = line[3].split(".")[0]
            #check where the motif positions are
            motif_pos = calc_density_for_concat_several_c(motifs, seqs[names.index(trans_name)], motif_lengths)
            motif_pos = [i.flatten() for i in  motif_pos if len(i) > 0]
            if motif_pos:
                motif_pos = np.unique(np.concatenate(tuple(motif_pos)))
            else:
                motif_pos = []
            #write down those lines that overlap with motifs
            if complement:
                if not np.any(motif_pos == ORF_pos):
                    ofile.write(line_raw)
            else:
                if np.any(motif_pos == ORF_pos):
                    ofile.write(line_raw)

def get_allele_frequency(snp):
    '''
    Get the allele frequency of a SNP.
    '''
    alleles = []
    #I added in the extra condition just in case you have a SNP file that has gone through two overlaps
    #and thus has two overlap count columns at the end
    [alleles.extend(i.split('|')) for i in snp[15:] if i != "1"]
    alleles = [int(i) for i in alleles]

    # print(np.divide(sum(alleles), len(alleles)))
    return(np.divide(sum(alleles), len(alleles)))


def generate_pesudo_monomorphic_ptcs(ptc_file, index_fastas, output_file, seed=None, without_replacement=None, with_weighting=None):

    '''
    Generate a file of pseudo PTC mutations where infact the site is a monomorphic site.
    Give the 'new' PTC the same allele freqeuncies as the real PTC.
    with_weighting: give each chunk without a mutation a weighting dependingon how many of a particular nt there are in that chunk
    '''

    nts = ["A", "C", "G", "T"]
    index_files = {}
    exon_search = re.compile('^.+\.\d')
    if without_replacement:
        replace = False
    else:
        replace = True

    for nt in nts:
        names, indices = gen.read_fasta(index_fastas[nt])
        # clean up the exon name to remove (+) and (-)
        names = [re.search(exon_search, name).group(0) for name in names]
        indices = [x.split(',') for x in indices]
        if with_weighting:
            # get a ist of all indices
            concat_indices = list(itertools.chain.from_iterable(indices))
            # create the weightings
            weights = [len(x)/len(concat_indices) for x in indices]
        else:
            # give all chunks the same weighting
            weights = [1/len(indices) for x in indices]
        index_files[nt] = [names, indices, weights]

    # set the randomisation seed
    np.random.seed(seed)

    ptcs = gen.read_many_fields(ptc_file, "\t")
    with open(output_file, "w") as output:
        head = ptcs[0]
        head[7] = "sim_spos"
        head[11] = "sim_rel_chunk_pos"
        head[12] = "sim_status"
        output.write("{0}\n".format("\t".join(head)))

        # for each ptc
        for ptc in ptcs[1:]:
            aa = ptc[9]
            pseudo_ptc = ptc
            # choose a random exon chunk
            random_exon = np.random.choice(list(range(len(index_files[aa][0]))), 1, p=index_files[aa][2])[0]
            # choose a random position within that chunk
            random_pos = np.random.choice([p for p in index_files[aa][1][random_exon]], 1, replace=replace)[0]
            # outputs to file, keeping same allele frequencies
            pseudo_ptc[3] = index_files[aa][0][int(random_exon)]
            pseudo_ptc[11] = random_pos
            pseudo_ptc[12] = "pseudo_ptc_snp"
            output.write('{0}\n'.format("\t".join(pseudo_ptc)))


def generate_pseudo_ptc_snps(input_ptc_snps, input_other_snps, ptc_output_file, other_snps_file, without_replacement=None, match_allele_frequency=None, match_allele_frequency_window=None, group_by_gene=None, seed=None):
    '''
    Generate a new file of pseudo PTC snps that are instead snps of different type.
    For each PTC snp in input_ptc_snps, take a random snp from the alternative file
    ensuring the ancestral and derived alleles match.
    Need to also then remove that snp from the file
    replacement: random choice with/without replacement
    match_allele_frequency: match the allele frequencies (ptcs are likely rare whereas alternative snps may be more common)
    match_allele_frequency_window: the proporition size of the window around the allele frequency to match, e.g 0.05 for allele frequency of 0.2 is 0.15-0.25.
    seed: list of seeds (must be greater or equal to the number of simulations)
    '''

    if (match_allele_frequency and not match_allele_frequency_window) or (match_allele_frequency_window and not match_allele_frequency):
        print("_match_allele_frequency_ and _match_allele_frequency_window_ must both be set")
        raise Exception

    #set up a default dictionary to hold indices of positions in list of alternative snps,
    #grouped by gene, ancestral base, mutated base
    alternative_snp_indices = collections.defaultdict(lambda: (collections.defaultdict(lambda: collections.defaultdict(lambda: []))))
    alternative_snps = gen.read_many_fields(input_other_snps, "\t")
    ptc_snps = gen.read_many_fields(input_ptc_snps, "\t")

    #check for header of files
    alt_header = False
    ptc_header = False
    if alternative_snps[0][0] == "echrom":
        alt_header = "\t".join(alternative_snps[0])
    if ptc_snps[0][0] == "echrom":
        ptc_header = "\t".join(ptc_snps[0])

    #go through each of the alternative snps and add to dictionary
    for i,snp in enumerate(alternative_snps):

        #check if there is a header
        if i == 0 and alt_header:
            pass
        else:
            #if grouping by gene, set the gene id
            if group_by_gene:
                gene_id = snp[3].split('.')[0]
            else:
                gene_id = "all"

            #if match allele frequency, calculate the allele frequency for each ptc snp
            if match_allele_frequency:
                snp_index = [i, get_allele_frequency(snp)]
            else:
                snp_index = [i]

            #add to dictionary: alternative_snp_indices[gene_id][ancestral_base][derived_base] = [[snp_index, allele_freqency], [snp_index, allele_freqency],...]
            #index 9 = ancestral base, #index 10 = mutation base
            alternative_snp_indices[gene_id][alternative_snps[i][9]][alternative_snps[i][10]].append(snp_index)
            #building a second dictionary that doesn't distinguish by gene so that we could handle genes that don't have non-synonymous SNPs
            if group_by_gene:
                alternative_snp_indices["all"][alternative_snps[i][9]][alternative_snps[i][10]].append(snp_index)

    #set seed for randomisation
    np.random.seed(seed)

    #count how many genes have no non-synonymous SNPs
    empty_counter = 0
    total_counter = 0

    #create an empty list to hold the alternative snps chosen
    pseudo_ptc_indices = {}
    #get the real ptc snps
    ptc_snps = gen.read_many_fields(input_ptc_snps, "\t")
    #backwards logic but makes more sense in the flags
    replacement = not without_replacement
    for i, ptc in enumerate(ptc_snps):
        total_counter = total_counter + 1
        #check for header
        if i == 0 and ptc_header:
            pass
        else:

            ptc_id = ptc[14]

            #if group_by_gene, get the dict for the particular gene id
            if group_by_gene:
                gene_id = ptc[3].split('.')[0]
            else:
                gene_id = 'all'


            #if match_allele_frequency, only allow snps with similar allele frequency defined by window
            if match_allele_frequency:
                #get the allele frequency of the ptc
                ptc_allele_frequency = get_allele_frequency(ptc)
                #set the upper and lower bounds
                alt_snp_allele_frequency_lower_limit = ptc_allele_frequency - match_allele_frequency_window
                alt_snp_allele_frequency_upper_limit = ptc_allele_frequency + match_allele_frequency_window
                #get all alternative snps with allele frequenecies within those bounds
                alt_snp_choices = [i for i in alternative_snp_indices[gene_id][ptc[9]][ptc[10]] if alt_snp_allele_frequency_lower_limit <= i[1] and i[1] <= alt_snp_allele_frequency_upper_limit]
            else:
                alt_snp_choices = alternative_snp_indices[gene_id][ptc[9]][ptc[10]]

            #check if there are any alternative snps
            empty_gene = False
            if len(alt_snp_choices) == 0:
                empty_gene = True
                empty_counter = empty_counter + 1

                if match_allele_frequency:
                    alt_snp_choices = [i for i in alternative_snp_indices["all"][ptc[9]][ptc[10]] if alt_snp_allele_frequency_lower_limit <= i[1] and i[1] <= alt_snp_allele_frequency_upper_limit]
                else:
                    alt_snp_choices = alternative_snp_indices["all"][ptc[9]][ptc[10]]

            #have to do it this way round because we have a list of lists, not list of items
            #which numpy doesnt like
            #generate a list of indices for the snp choices
            alt_snp_choices_indices = range(len(alt_snp_choices))
            #pick one of those indices
            alt_snp_choice_index = np.random.choice(alt_snp_choices_indices, 1)[0]
            #add the index of the snp to the dictionary of chosen snps
            #it's a dictionary of lists because some SNPs might be sampled several times when sampling
            #with replacement
            chosen_abs_index = alt_snp_choices[alt_snp_choice_index][0]
            if chosen_abs_index not in pseudo_ptc_indices:
                pseudo_ptc_indices[chosen_abs_index] = []
            pseudo_ptc_indices[chosen_abs_index].append(ptc_id)
            #if no replacement, remove snp from the snp choices list
            if (not replacement) and (not empty_gene):
                del alternative_snp_indices[gene_id][ptc[9]][ptc[10]][alt_snp_choice_index]

##    print("Genes with no non-synonymous SNPs: {0}/{1} ({2}%).".format(empty_counter, total_counter, round(empty_counter / total_counter * 100, 3)))

    #open a pesudo ptc snps output file and other snps output file
    pseudo_ptc_output = open(ptc_output_file, "w")
    other_snps_output = open(other_snps_file, "w")

    #write header to both files
    pseudo_ptc_output.write("{0}\n".format(ptc_header))
    other_snps_output.write("{0}\n".format(alt_header))

    #for each alternative snp, write to ptc file if in the list or other file if not
    for i, alternative_snp in enumerate(alternative_snps):
        if i == 0 and alt_header:
            pass
        else:
            if i in pseudo_ptc_indices:
                #need this because if we are doing with replacement, the snp may appear more than once
                for ptc_id in pseudo_ptc_indices[i]:
                    snp_to_write = alternative_snps[i].copy()
                    snp_to_write[14] = ptc_id
                    pseudo_ptc_output.write("{0}\n".format("\t".join(snp_to_write)))
            else:
                other_snps_output.write("{0}\n".format("\t".join(alternative_snps[i])))

    pseudo_ptc_output.close()
    other_snps_output.close()

def get_codon_start(variant_pos, seq_len, shift = False):
    '''
    Given the (0-based) position of a SNP in a CDS, determine what (0-based) position is the 1st in the overlapping codon.
    seq_len: length of the CDS
    shift: whether or not the codon position should be shifted (see get_snp_type
    for further details)
    '''
    codon_start = (variant_pos//3) * 3
    if not shift:
        return(codon_start)
    else:
        #if the SNP overlaps the first position in the codon,
        #then you can't shift forward by a base or your new codon wouldn't
        #include the SNP
        if variant_pos == codon_start:
            codon_start = codon_start - 1
            if codon_start < 0:
                return("error")
            return(codon_start)
        codon_start = codon_start + 1
        if (codon_start + 3) >= seq_len:
            return("error")
        return(codon_start)

def get_snp_relative_cds_position(snp_exon_relative_positions, snp_cds_position_output, full_bed):
    '''
    Get the position of the snp within a CDS using the relative positions of snps in the features they are found
    '''
    #read in the coordinates from the full_bed
    bed_data = gen.read_many_fields(full_bed, "\t")
    cds_exons = {}
    #iterate through the bed entries, sorting by exon (needed in case bed is not sorted)
    #NB! check if two of the exons have the same end coordinates
    #if yes, delete the shorter one
    #if not, concatenate
    for i, bed_line in enumerate(bed_data):
        name_field = bed_line[3].split(".")
        trans_name = name_field[0]
        exon = int(name_field[1])
        #store the length of the exon
        length = int(bed_line[2]) - int(bed_line[1])
        if trans_name not in cds_exons:
            cds_exons[trans_name] = {}
        if exon in cds_exons[trans_name]:
            cds_exons[trans_name][exon].append([bed_line, length])
        else:
            cds_exons[trans_name][exon] = [[bed_line, length]]

    #when you have the same exon number for two features (CDS and stop),
    #sometimes the relative exon position will need to be increased to account for this
    #I will make a dicionary that contains such cases
    add_to_rel_coords = {}

    entry_regex = re.compile("(\w+)\.(\d+)(\..*)*")
    #set up dict to hold the feature positions relative to the cds
    cds_features_relative_positions = collections.defaultdict(lambda: collections.defaultdict())
    #get the relative start positions of each exon
    for cds in cds_exons:
        total = 0
        for exon in sorted(cds_exons[cds]):
            exon_info = cds_exons[cds][exon]
            if len(exon_info) == 1:
                to_add = exon_info[0][1]
            elif len(exon_info) == 2:
                strand = exon_info[0][0][5]
                if strand == "+":
                    ends = [int(i[0][2]) for i in exon_info]
                elif strand == "-":
                    ends = [int(i[0][1]) for i in exon_info]
                else:
                    print("Problem with strand information!")
                    print(exon_info)
                    raise Exception
                if ends[0] == ends[1]:
                    #this is the case where a stop has been split in two
                    #and so the first chunk overlaps with the last CDS feature (same region annotated twice)
                    #and the second chunk is annotated only as stop
                    #so you need to make sure you don't count the overlapping region twice
                    to_add = max([exon_info[0][1], exon_info[1][1]])
                    #this bit it because you don't want to count the SNPs in the overlapping bits twice
                    if strand == "+":
                        if int(exon_info[0][0][1]) > int(exon_info[1][0][1]):
                            #means the stop is the first element
                            #store both the coordinates of the stop bit
                            #and the length of the CDS bit
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[0][0], "remove"]
                        else:
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[1][0], "remove"]
                    elif strand == "-":
                        #means the stop is the first element
                        if int(exon_info[0][0][2]) < int(exon_info[1][0][2]):
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[0][0], "remove"]
                        else:
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[1][0], "remove"]
                else:
                    #this is the case of a normal stop that appears after
                    #a CDS so you add up the lengths for the total
                    #and store the length of the CDS chunk in the dictionary for coords
                    #where the relative exon position must be augmented
                    to_add = exon_info[0][1] + exon_info[1][1]
                    if strand == "+":
                        if ends[0] > ends[1]:
                            #means the stop is the first element
                            #store both the coordinates of the stop bit
                            #and the length of the CDS bit
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[0][0], exon_info[1][1]]
                        else:
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[1][0], exon_info[0][1]]
                    elif strand == "-":
                        if ends[0] < ends[1]:
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[0][0], exon_info[1][1]]
                        else:
                            add_to_rel_coords[".".join([cds, str(exon)])] = [exon_info[1][0], exon_info[0][1]]
            else:
                print("This makes no sense, there's more than two elements with the same exon number!")
                print(cds_exons[cds])
                raise Exception
            cds_features_relative_positions[cds][exon] = total

            total += to_add

    #now get the relative position of the snp within the cds
    with open(snp_cds_position_output, "w") as output:
        error_count = 0
        for snp in snp_exon_relative_positions:
            cds_id = snp[3]
            cds_id_meta = re.search(entry_regex, cds_id)
            snp_cds = cds_id_meta.group(1)
            snp_exon_id = cds_id_meta.group(2)
            snp_exon_position = int(snp[11])
            #snp position in cds is the position in the exon plus the exons position in the cds
            #plus handling stops
            full_ID = cds_id_meta.group(0)
            extra = 0
            if full_ID in add_to_rel_coords:
                #this check is to make sure that you're processing the stop and not the CDS bit from the same exon
                if (snp[1] == add_to_rel_coords[full_ID][0][1]) and (snp[2] == add_to_rel_coords[full_ID][0][2]):
                    extra = add_to_rel_coords[full_ID][1]
            if extra != "remove":
                try:
                    snp_cds_position = snp_exon_position + cds_features_relative_positions[snp_cds][int(snp_exon_id)] + extra
                    snp[11] = str(snp_cds_position)
                    output.write("{0}\n".format("\t".join(snp)))
                except KeyError:
                    error_count = error_count + 1
    print("Errors: {0}.".format(error_count))

def get_snp_relative_exon_position(intersect_file, snp_relative_exon_position_file):
    '''
    Get the relative position of a snp within the exon it is found. Used as an intermediate step
    before calculating the snp position in the cds using get_snp_cds_relative_position.
    '''
    #read the intersects file
    cds_strands = collections.defaultdict()
    relative_positions = []
    intersects = gen.read_many_fields(intersect_file, "\t")
    for intersect in intersects:
        #get strand of exon
        entry_regex = re.compile("(\w+)\.(\d+)(\..*)*")
        trans = re.search(entry_regex, intersect[3]).group(1)
        strand = intersect[5]
        cds_strands[trans] = strand
        #get the features
        feature_start = intersect[1]
        feature_end = intersect[2]
        feature = intersect[3]
        #-1 because VCF files are 1-based
        snp_start = int(intersect[7]) - 1
        #get the position of the snp compared with the feature
        if strand == "+":
            relative_position = int(snp_start) - int(feature_start)
        elif strand == "-":
            #-1 because we're in base 0
            relative_position = (int(feature_end) - 1) - snp_start
        else:
            print("Incorrect strand: {0}.".format(strand))
            print(intersect)
            raise Exception
        #replace . field with relative position
        intersect[11] = str(relative_position)
        relative_positions.append(intersect)

    # write relative snp positions to output file
    with open(snp_relative_exon_position_file, "w") as outfile:
        for snp in relative_positions:
            outfile.write("{0}\n".format("\t".join(snp)))

    return(relative_positions)

def get_snp_change_status(snp_cds_relative_positions, cds_fasta, ptcs_output_file, others_output_file, out_of_frame = False):
    '''
    For a set of SNPs, determine the effect of the PTC (nonsense, missense, synonymous).
    Store the nonsense ones in one file and all the rest in another.
    '''

    snps = gen.read_many_fields(snp_cds_relative_positions, "\t")
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    entry_regex = re.compile("(\w+)\.(\d+)(\..*)*")
    var_type_reg = re.compile("VT=([A-z]+)")
    ancestral_reg = re.compile("AA=([A-z]+)")

    ptc_id_counter = 1
    other_id_counter = 1
    with open(ptcs_output_file, "w") as ptc_outputs, open(others_output_file, "w") as other_outputs:
        refbase_error = 0
        snp_count = 0
        #the first line is the header
        header = "{0}\n".format("\t".join(snps[0]))
        ptc_outputs.write(header)
        other_outputs.write(header)

        print(snps[1:5])
        for snp in snps[1:]:
            cds_id = re.search(entry_regex, snp[3]).group(1)
            snp_index = int(snp[11])
            #get the strand
            strand = snp[5]
            #get ancestral base
            ref_base = snp[9]
            #get the information on the variant type
            var_base = snp[10].split(",")
            var_base_count = len(var_base)
            var_base = [i for i in var_base if i in ["A", "C", "G", "T"]]
            ancestral_allele = re.search(ancestral_reg, snp[13])

            snp_id = snp[8]

            #get the feature type
            var_type = re.search(var_type_reg, snp[13])
            if var_type:
                var_type = var_type.group(1)

            #check whether the cds is in the fasta (can be after the quality control filterings)
            if cds_id in cds_names:
                snp_count += 1
                #check that the snp is only one base
                if len(ref_base) == 1:
                    #from RS: filter out polymorphisms with more than 2 segregating alleles
                    #need to check the number both before and after filtering out non-canonical bases
                    #also check whether the variant type is annotated as a snp
                    if var_base_count == 1 and len(var_base) == 1 and var_type and var_type == "SNP":
                        var_base = var_base[0]

                        if strand == "-":
                            ref_base = gen.reverse_complement(ref_base)
                            var_base = gen.reverse_complement(var_base)
                        #get the base of reference cds where the snp occured
##                        print(snp)
##                        print(cds_seqs[cds_names.index(cds_id)])
##                        print(len(cds_seqs[cds_names.index(cds_id)]))
##                        print("\n")
                        cds_base = cds_seqs[cds_names.index(cds_id)][snp_index]

                        #check whether cds base and ref base are the same
                        if cds_base != ref_base:
                            refbase_error += 1
##                            print("Cds base and reference base not the same (id: {0})".format(snp[8]))
##                            print("Cds base: {0}".format(cds_base))
##                            print("Ref base: {0}".format(ref_base))
##                            print("Variant base: {0}".format(var_base))
##                            print("\n")
                            pass
                        else:
                            cds_codon, snp_codon, mutation_type = get_snp_type(cds_seqs[cds_names.index(cds_id)], [snp_index, var_base], snp_id, shift = out_of_frame)

                            if ancestral_allele:
                                aa = ancestral_allele.group(1)
                            else:
                                aa = "UNDEFINED"

                            snp[13] = "CDS_CODON={0}$SNP_CODON={1}$AA={2}".format(cds_codon, snp_codon, aa)
                            snp[12] = mutation_type
                            if(mutation_type == "ptc"):
                                snp[14] = str(ptc_id_counter)
                                ptc_id_counter = ptc_id_counter + 1
                                ptc_outputs.write("{0}\n".format("\t".join(snp)))
                            else:
                                snp[14] = str(other_id_counter)
                                other_id_counter = other_id_counter + 1
                                other_outputs.write("{0}\n".format("\t".join(snp)))

    if snp_count:
        print("Number of ref errors: {0}/{1} ({2}%)".format(refbase_error, snp_count, np.divide(refbase_error, snp_count)*100))
    else:
        print("No SNPs were extracted!")
        raise Exception

def get_snp_type(sequence, variant, snp_id, shift = False):
    '''
    Get the effect of a particular SNP (nonsense/missense/synonymous).
    If shift is True (when doing an out-of-frame simulation),
    the codon that each SNP is in will be shifted 3' by a single base
    (for the purposes of determining the type of SNP),
    except if this would move the SNP out of the codon or the codon further than the sequence end, in which case
    the codon is shifted one base 5' instead.
    '''

    codon_map = {
        "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
        "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
        "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
        "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGT":"s", "AGC":"s", "AGA":"r", "AGG":"r",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    }



    shift_amount = 0

    #get the sequence with the snp in position
    snp_sequence = sequence[:int(variant[0])] + variant[1] + sequence[int(variant[0])+1:]
    #extract the SNP codon both with the reference and the SNP allele
    codon_start = get_codon_start(int(variant[0]), len(sequence), shift = shift)
    if codon_start != "error":
        if snp_id == "rs74315366":
            print(sequence[:int(variant[0])], variant[1], sequence[int(variant[0])+1:])
        cds_codon = sequence[codon_start:codon_start + 3]
        snp_codon = snp_sequence[codon_start:codon_start + 3]

        if cds_codon != snp_codon:
            #determine the type of snp
            #if the snp generated an in frame stop
            if codon_map[snp_codon] == "*":
                mutation_type = "ptc"
            #if the snp generated a synonymous codon
            elif codon_map[cds_codon] == codon_map[snp_codon]:
                mutation_type = "syn"
            #if the snp generated a nonsynonymous codon
            elif codon_map[cds_codon] != codon_map[snp_codon] and codon_map[snp_codon] != "*":
                mutation_type = "non"
            #error calling the snp (shouldn't occur)
            else:
                mutation_type = "call_error"
        else:
            cds_codon, snp_codon, mutation_type = "error", "error", "error"
    else:
        cds_codon, snp_codon, mutation_type = "error", "error", "error"

    return cds_codon, snp_codon, mutation_type

def get_snps_in_cds(bed, full_bed, vcf_folder, panel_file, names, sample_file, output_file, out_prefix):
    '''
    Given a bed file of CDS regions (with the corresponding interval fasta), a .vcf file, a panel file and a set of sample identifiers from 1000Genomes,
    pick out SNPs that overlap with any of the bed intervals in any of the selected samples and calculate their relative
    position in the full ORF. full_bed is a bed_file that has all of the exons from the relevant transcripts,
    whereas the bed might only have some exons. Also write a filtered vcf to sample_file.
    '''
    #get the relevant SNPs
    tabix_samples(bed, sample_file, panel_file, vcf_folder, samples = names, chr_prefix = True, remove_empty = True, exclude_xy = True)
    #the tabix_samples and the intersec-bed are kind of redundant
    #however, this way we have a proper vcf as a result of tabix_samples that we can query using tabix
    #and we have intersect_bed, which has both the SNP and the exon information in a nice format
    intersect_file = "{0}_CDS_SNP_intersect.bed".format(out_prefix)
    bmo.intersect_bed(bed, sample_file, write_both = True, output_file = intersect_file, no_dups = False)


def get_snp_positions(sample_file, output_file, out_prefix):
    intersect_file = "{0}_CDS_SNP_intersect.bed".format(out_prefix)
    snp_relative_exon_position_file = "{0}_SNP_relative_exon_position.bed".format(out_prefix)
    exon_pos = get_snp_relative_exon_position(intersect_file, snp_relative_exon_position_file)
    #this last bit is just to add a header to the final output file
    #so you'd know which sample is which
    header_line = "echrom\testart\teend\teID\tfeature\tstrand\t#schr\tspos\tsID\taa\tma\trel_pos\tstatus\tinfo\tformat"
    samples_header = gen.run_process(["grep", "CHROM", "{0}_uncompressed.txt".format(sample_file)])
    samples_header = re.search("(?<=FORMAT)[\tA-Z0-9]*", samples_header).group(0)
    header_line = header_line + samples_header + "\toverlap_count\n"
    temp_file = "temp_data/temp{0}.txt".format(random.random())
    temp_file2 = "temp_data/temp{0}.txt".format(random.random())
    with open(temp_file, "w") as file:
        file.write(header_line)
    gen.run_process(["cat", temp_file, output_file], file_for_output = temp_file2)
    gen.run_process(["mv", temp_file2, output_file])
    gen.remove_file(temp_file)

def merge_and_header(file1, file2, out_file):
    '''
    Merge two SNP files and add a header so it'd be recognized as VCF.
    '''
    temp1 = "temp_data/temp{0}.txt".format(random.random())
    temp2 = "temp_data/temp{0}.txt".format(random.random())
    header = "##fileformat=VCF\n"
    with open(temp1, "w") as file:
        file.write(header)
    #remove header from second file
    gen.run_process(["tail", "-n", 2, file2], file_for_output = temp2)
    #put header, file1 and file2 together
    gen.run_process(["cat", temp1, file1, temp2], file_for_output = out_file)
    gen.remove_file(temp1)
    gen.remove_file(temp2)

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
    #this is how the output files have been named by the core tabix function
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

def tabix_samples(bed_file, output_file_name, panel_file, vcf_folder, superpop = None, subpop = None, samples = None, downsample_by = None, exclude_xy = False, chr_prefix = False, remove_empty = False):
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
    chr_prefix: if True, prefix "chr" to chromosome names in the final output file
    remove_empty: if True, don't report SNPs that don't appear in any of the selected samples or that appear in all of them (i.e. allele frequency 1)
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

    #get a list of all the files in the VCF folder
    vcf_files = os.listdir(vcf_folder)
    #just in case there is nonsense in the directory
    vcf_files = [i for i in vcf_files if "vcf" in i]

    # with open(bed_file) as file:
    #     sample_files = []
    #     counter = 0
    #     #loop over lines in bed file
    #     for line in file:
    #         #print out every 100th line number
    #         counter = gen.update_counter(counter, 500, "Bed lines processed: ")
    #         #parse line in bed file
    #         line = line.split("\t")
    #         chrom = line[0].lstrip("chr")
    #         if chrom in sex_chromosomes and exclude_xy:
    #             pass
    #         else:
    #             #add 1 to start coordinate because bed files are 0-based, whereas the vcf files are 1-based
    #             start = int(line[1]) + 1
    #             end = line[2]
    #             trans = line[3]
    #             #get the vcf file for the right chromosome, making sure not to get the index file instead
    #             current_vcf = ["{0}/{1}".format(vcf_folder, i) for i in vcf_files if "chr{0}.".format(chrom) in i and ".tbi" not in i]
    #             #check that you only got a single file
    #             if len(current_vcf) != 1:
    #                 print("Ambiguous or missing files in VCF folder!")
    #                 print(current_vcf)
    #                 print(chrom)
    #                 raise Exception
    #             else:
    #                 current_vcf = current_vcf[0]
    #             #generate temporary output file for all SNPs in interval
    #             temp_output_file = "temp_data/temp_vcf{0}.vcf".format(random.random())
    #             #get ALL SNPs (that is to say, for all samples) for current interval
    #             gen.run_process(["tabix", "-h", current_vcf, "{0}:{1}-{2}".format(chrom, start, end)], file_for_output = temp_output_file)
    #             #uncomment the following line for debug
    #             # gen.run_process(["cp", temp_output_file, "temp_data/{0}:{1}-{2}_tabix_slice.txt".format(chrom, start, end)])
    #             #generate temporary output file for SNPs from your seleceted samples
    #             sample_output_file = "temp_data/temp_sample_tabix{0}.txt".format(random.random())
    #             sample_files.append(sample_output_file)
    #             #filter the file you made with all the SNPs to only leave the SNPs that appear in your samples
    #             gen.run_process(["vcf-subset", "-c", samples, temp_output_file], file_for_output = sample_output_file)
    #             gen.remove_file(temp_output_file)

    sample_files = []
    for file in os.listdir('temp_data/'):
        if file.startswith('temp_sample_tabix'):
            path = 'temp_data/{0}'.format(file)
            sample_files.append(path)


    # you want to concatenate the sample files you made (one file per bed interval) but you can't in one go cause there's too many
    # therefore, you take the 10 last files, concatenate those
    # then concatenate the next 10 files (moving from the end of the list towards the beginning) to each-other and to the file you got in the previous step
    # etc.
    # you juggle the two temp concat file names just so you would be overwriting files rather than creating new ones
    print('Concatenating files...')
    concat_files = ["temp_data/temp_concat_file{0}.vcf".format(random.random()), "temp_data/temp_concat_file{0}.vcf".format(random.random())]
    current_sample_files = sample_files[-10:]
    del sample_files[-10:]
    gen.run_process(["vcf-concat"] + current_sample_files, file_for_output = concat_files[0])
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
        gen.run_process(["vcf-concat"] + current_sample_files + [previous_concat_file], file_for_output = current_concat_file)
    sort_file = "{0}_uncompressed.txt".format(output_file_name)
    #once everything is concatenated, sort the SNPs, prefix "chr" if needed, make a compressed version of the file and make an index for tabix
    print('Sort SNPs, prefix and compress for tabix...')
    gen.run_process(["vcf-sort", current_concat_file], file_for_output = sort_file)
    if chr_prefix or remove_empty:
        allele_regex = re.compile("[0-9]+\|[0-9]+")
        temp_file = "temp_data/temp{0}.txt".format(random.random())
        with open(sort_file) as infile, open(temp_file, "w") as outfile:
            for line in infile:
                dont_write = False
                if line[0] != "#":
                    if chr_prefix:
                        line = "chr" + line
                    if remove_empty:
                        alleles = "".join(re.findall(allele_regex, line))
                        alleles_ones = [i for i in alleles if i != "0" and i != "|"]
                        alleles_zeroes = [i for i in alleles if i != "1" and i != "|"]
                        if not alleles_ones or not alleles_zeroes:
                            dont_write = True
                if not dont_write:
                    outfile.write(line)
        gen.run_process(["mv", temp_file, sort_file])
    gen.run_process(["bgzip", "-c", sort_file], file_for_output = output_file_name)
    print('Run tabix...')
    gen.run_process(["tabix", "-f", "-p", "vcf", output_file_name])
    # #clean up
    # for sample_file in sample_files:
    #     gen.remove_file(sample_file)
    # for concat_file in concat_files:
    #     gen.remove_file(concat_file)
