'''
Auhor: Liam Abrahams.
Look at disease SNPs from ClinVar.
'''

import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os
import collections
import re
import copy
import numpy as np
from scipy.stats import chisquare

class Define_Exon(object):
    def __init__(self, exon):
        self.chr = exon[0]
        self.start = int(exon[1])
        self.stop = int(exon[2])
        self.transcript_id = exon[3].split('.')[0]
        self.exon_no = int(exon[3].split('.')[1])
        self.type = exon[4]
        self.strand = exon[5]

class EntryInfo(object):
    def __init__(self, entry):
        self.chr = entry[0]
        self.exon_start = int(entry[1])
        self.exon_stop = int(entry[2])
        self.transcript_id = entry[3].split('.')[0]
        self.exon_id = int(entry[3].split('.')[1])
        self.strand = entry[5]
        self.snp_pos = int(entry[7])
        self.snp_id = int(entry[8])
        self.aa = entry[9]
        self.ma = entry[10]
        self.rel_pos = int(entry[11])
        self.type = entry[12]
        self.given_id = entry[-1]

def disease_ptc_location_test(rel_pos_file, coding_exons_fasta, output_file):
    '''
    Chisquare test on ptc locations
    '''

    # get coding exons
    coding_exons = get_coding_exons(coding_exons_fasta)
    # get relative positions
    relative_positions = get_relative_position_list(rel_pos_file)

    ptc_locations = get_ptc_positions(relative_positions, coding_exons)
    exon_nts = get_exon_nts(relative_positions, coding_exons)

    expected = []
    for i in exon_nts:
        expected_ptcs = np.divide(i, sum(exon_nts))*sum(ptc_locations)
        expected.append(expected_ptcs)

    chisq = chisquare(ptc_locations, expected_ptcs)

    positions = ["0.2bp", "3.69bp", "70+bp"]
    with open(output_file, "w") as outfile:
        outfile.write('position,num_nts,num_ptcs\n')
        for i, pos in enumerate(positions):
            outfile.write("{0},{1},{2}\n".format(pos, exon_nts[i], ptc_locations[i]))
        outfile.write('sum,{0},{1}\n'.format(sum(exon_nts), sum(ptc_locations)))

        outfile.write('\n,nts%,ptcs%\n')
        for i, pos in enumerate(positions):
            outfile.write("{0},{1},{2}\n".format(pos, np.divide(exon_nts[i], sum(exon_nts))*100, np.divide(ptc_locations[i], sum(ptc_locations))*100))

        outfile.write('\n,expected,observed\n')
        for i, pos in enumerate(positions):
            outfile.write("{0},{1},{2}\n".format(pos, expected[i], ptc_locations[i]))

        outfile.write("df,chisq,p\n")
        outfile.write("{0},{1},{2}\n".format(len(ptc_locations)-1,chisq.statistic, chisq.pvalue))

def ese_hit_simulation(rel_pos_file, coding_exons_fasta, simulations, output_file, ese_file, window_start, window_end, exclude_cpg, clinvar=None):
    '''
    Simulate ese hits strictly within a region
    '''

    # get a list of the relative positions of the ptcs
    relative_positions_list = get_relative_position_list(rel_pos_file)
    # get a list of eses
    ese_list = get_eses_from_file(ese_file)
    # get the coding exons
    coding_exons = get_coding_exons(coding_exons_fasta)

    long_exons = get_long_exons(relative_positions_list, coding_exons, window_end*2)

    # get the ptcs that are in the 3-69 bp region for each exon of exon
    # this requires exons at least 128 bp in length for comparison
    window_ptcs = get_ptcs_in_window(long_exons, window_start, window_end, coding_exons)
    real_ese_hits = get_ese_hits(window_ptcs, coding_exons, ese_list)

    # simulate the hit counts for nt matched mutations
    simulation_list = list(range(simulations))
    # simulate_ese_hits(simulation_list, simulations, window_ptcs, coding_exons_fasta, ese_list, window_start, window_end)
    processes = gen.run_in_parallel(simulation_list, ["foo", simulations, window_ptcs, coding_exons_fasta, ese_list, window_start, window_end], simulate_ese_hits)
    #
    simulation_outputs = {}
    for process in processes:
        simulation_hits = process.get()
        simulation_outputs = {**simulation_outputs, **simulation_hits}


    with open(output_file, "w") as outfile:
        outfile.write("simulation,ese_hit_count,cant_count\n")
        outfile.write("real,{0},0\n".format(real_ese_hits))
        for simulation in sorted(simulation_outputs):
            outlist = [simulation+1, simulation_outputs[simulation][0], simulation_outputs[simulation][1]]
            outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))

def excess_test(unique_ptcs_rel_pos_file, coding_exons_fasta, output_file):
    '''
    Ask whether there is an excess of nonsense mutations in the exon ends
    compared with the exon cores
    '''

    # get coding exons
    coding_exons = get_coding_exons(coding_exons_fasta)
    # get relative positions
    relative_positions = get_relative_position_list(unique_ptcs_rel_pos_file)

    # keep only those with exons longer than 138
    long_exons = get_long_exons(relative_positions, coding_exons, 138)
    short_exons = {}
    for name in relative_positions:
        if name not in long_exons:
            short_exons[name] = relative_positions[name]

    # get the ptc locations and exon nts for both long and short exons
    long_ptc_locations = get_ptc_positions(long_exons, coding_exons)
    long_exon_nts = get_exon_nts(long_exons, coding_exons)
    short_ptc_locations = get_ptc_positions(short_exons, coding_exons)
    short_exon_nts = get_exon_nts(short_exons, coding_exons)

    exon_cores_ptc_per_nt = np.divide(long_ptc_locations[2], long_exon_nts[2])

    total_nts = [x + y for x, y in zip(long_exon_nts, short_exon_nts)]
    total_ptcs = [x + y for x, y in zip(long_ptc_locations, short_ptc_locations)]

    expected_0_2 = total_nts[0]*exon_cores_ptc_per_nt
    excess_0_2 = total_ptcs[0] - expected_0_2
    excess_percentage_0_2 = np.divide(excess_0_2, sum(total_ptcs))*100

    expected_3_69 = total_nts[1]*exon_cores_ptc_per_nt
    excess_3_69 = total_ptcs[1] - expected_3_69
    excess_percentage_3_69 = np.divide(excess_3_69, sum(total_ptcs))*100

    expected_large_0_2 = long_exon_nts[0]*exon_cores_ptc_per_nt
    excess_large_0_2 = long_ptc_locations[0] - expected_large_0_2
    excess_large_percentage_0_2 = np.divide(excess_large_0_2, sum(long_ptc_locations))*100

    expected_large_3_69 = long_exon_nts[1]*exon_cores_ptc_per_nt
    excess_large_3_69 = long_ptc_locations[1] - expected_large_3_69
    excess_large_percentage_3_69 = np.divide(excess_large_3_69, sum(long_ptc_locations))*100

    with open(output_file, "w") as outfile:
        outfile.write("nts in exons > 138 nt\n")
        outfile.write("0.2,{0}\n".format(long_exon_nts[0]))
        outfile.write("3.69,{0}\n".format(long_exon_nts[1]))
        outfile.write("70+,{0}\n".format(long_exon_nts[2]))
        outfile.write("\n")
        outfile.write("\nnts in exons <= 138\n")
        outfile.write("0.2,{0}\n".format(short_exon_nts[0]))
        outfile.write("3.69,{0}\n".format(short_exon_nts[1]))
        outfile.write("70+,{0}\n".format(short_exon_nts[2]))
        outfile.write("\n")
        outfile.write("\nTotal nts in exons\n")
        outfile.write("0.2,{0}\n".format(total_nts[0]))
        outfile.write("3.69,{0}\n".format(total_nts[1]))
        outfile.write("70+,{0}\n".format(total_nts[2]))

        outfile.write("\n")
        outfile.write("\n***\n")
        outfile.write("\n")

        outfile.write("\nnonsense mutations > 138 nt\n")
        outfile.write("0.2,{0}\n".format(long_ptc_locations[0]))
        outfile.write("3.69,{0}\n".format(long_ptc_locations[1]))
        outfile.write("70+,{0}\n".format(long_ptc_locations[2]))
        outfile.write("\n")
        outfile.write("\nnonsense mutations <= 138 nt\n")
        outfile.write("0.2,{0}\n".format(short_ptc_locations[0]))
        outfile.write("3.69,{0}\n".format(short_ptc_locations[1]))
        outfile.write("70+,{0}\n".format(short_ptc_locations[2]))
        outfile.write("\n")
        outfile.write("\nTotal nonsense mutations\n")
        outfile.write("0.2,{0}\n".format(total_ptcs[0]))
        outfile.write("3.69,{0}\n".format(total_ptcs[1]))
        outfile.write("70+,{0}\n".format(total_ptcs[2]))
        outfile.write("\n")
        outfile.write("\nTotal nonsense mutations,{0}\n\n".format(sum(total_ptcs)))

        outfile.write("\n")
        outfile.write("\n***\n")
        outfile.write("\n")

        outfile.write("nonsense per nt in cores (>138),{0}\n".format(exon_cores_ptc_per_nt))

        outfile.write("\n")
        outfile.write("\n***\n")
        outfile.write("\n")

        outfile.write("\nexcesses in exons > 138\n")
        outfile.write("\n0.2 expected,{0}\n".format(expected_large_0_2))
        outfile.write("0.2 difference,{0}\n".format(excess_large_0_2))
        outfile.write("0.2 excess,{0}%\n".format(excess_large_percentage_0_2))
        outfile.write("\n")
        outfile.write("\n3.69 expected,{0}\n".format(expected_large_3_69))
        outfile.write("3.69 difference,{0}\n".format(excess_large_3_69))
        outfile.write("3.69 excess,{0}%\n".format(excess_large_percentage_3_69))

        outfile.write("\n")
        outfile.write("\n***\n")
        outfile.write("\n")
        outfile.write("\nexcesses in all exons\n")
        outfile.write("\n0.2 expected,{0}\n".format(expected_0_2))
        outfile.write("0.2 difference,{0}\n".format(excess_0_2))
        outfile.write("0.2 excess,{0}%\n".format(excess_percentage_0_2))
        outfile.write("\n")
        outfile.write("\n3.69 expected,{0}\n".format(expected_3_69))
        outfile.write("3.69 difference,{0}\n".format(excess_3_69))
        outfile.write("3.69 excess,{0}%\n".format(excess_percentage_3_69))


# def (intersect_file, ptc_file, relative_exon_positions_file, exon_fasta, cds_fasta, output_file):
#     '''
#     Get the location of ptcs and whether or not they are disease causing
#     '''
#     exon_list = collections.defaultdict(lambda: [])
#     exon_names, exon_seqs = gen.read_fasta(exon_fasta)
#     for i, name in enumerate(exon_names):
#         transcript = name.split('.')[0]
#         exon = int(name.split('.')[1])
#         if len(exon_seqs[i]) > 3:
#             exon_list[transcript].append(exon)
#
#     cds_list = {}
#     cds_names, cds_seqs = gen.read_fasta(cds_fasta)
#     for i, name in enumerate(cds_names):
#         cds_list[name] = cds_seqs[i]
#
#     disease_ptcs = gen.read_many_fields(intersect_file, "\t")
#     disease_ptc_ids = []
#     for ptc in disease_ptcs:
#         disease_ptc_ids.append(int(ptc[11]))
#
#     relative_position_entries = gen.read_many_fields(relative_exon_positions_file, "\t")
#     relative_positions = {}
#     for entry in relative_position_entries:
#         entry = entry[:12]
#         snp_id = entry[8]
#         exon_length = int(entry[2]) - int(entry[1])
#         rel_pos = entry[-1]
#         relative_positions[snp_id] = [entry[3], int(entry[7]), int(rel_pos), exon_length]
#
#
#     ptcs = gen.read_many_fields(ptc_file, "\t")
#     with open(output_file, "w") as outfile:
#         outfile.write("snp_id\ttranscript\texon_no\ttotal_exons\texon_length\trel_exon_position\texon_dist\texon_end\trel_cds_position\tcds_length\tref_allele\tmut_allele\tdisease\tnmd\n")
#
#         for ptc in ptcs[1:]:
#             snp_pos = int(ptc[7])
#             snp_id = ptc[8]
#             transcript_exon = ptc[3]
#             transcript = ptc[3].split('.')[0]
#             exon = int(ptc[3].split('.')[1])
#             snp_ref = int(ptc[14])
#             cds_pos = ptc[11]
#             ref_allele = ptc[9]
#             mut_allele = ptc[10]
#
#             relative_position_entry = relative_positions[snp_id]
#             exon_rel_pos = relative_position_entry[2]
#             exon_length = relative_position_entry[3]
#             dist_to_end = exon_length - exon_rel_pos
#
#             distances = [exon_rel_pos, dist_to_end]
#             min_dist = min(distances)
#
#             if distances.index(min_dist) == 0:
#                 end = 5
#             else:
#                 end = 3
#
#             total_exons = len(exon_list[transcript])
#             cds_length = len(cds_list[transcript])
#
#             if exon == sorted(exon_list[transcript])[-1] or exon == sorted(exon_list[transcript])[-2] and dist_to_end < 50:
#                 nmd = 0
#             else:
#                 nmd = 1
#
#             if snp_ref in disease_ptc_ids:
#                 disease = 1
#             else:
#                 disease = 0
#
#             outlist = gen.stringify([snp_id, transcript, exon, total_exons, exon_length, exon_rel_pos, min_dist, end, cds_pos, cds_length, ref_allele, mut_allele, disease, nmd])
#             outfile.write("{0}\n".format("\t".join(outlist)))

def generate_possible_ptc_locations(full_bed, cds_fasta, output_directory):
    '''
    Get all the locations of positions that if mutated to a
    particular nt, could generate a ptc
    '''

    stop_bases = ["A", "G", "T"]    # a mutation to c cant generate a stop so its not included
    stop_codons = ["TAA", "TAG", "TGA"]

    # get the list of coding sequnces
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {}
    for i, name in enumerate(cds_names):
        cds_list[name] = cds_seqs[i]

    # get all the indicies of mutable positions in the cds that could generate an in frame stop
    mutation_index_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for x, name in enumerate(cds_list):
        seq = cds_list[name]
        inframe_codons = re.findall('.{3}', seq[:-3])
        for i, codon in enumerate(inframe_codons):
            nts = list(codon)
            for j, nt in enumerate(nts):
                for mut in stop_bases:
                    if mut != nt:                       # we dont need 'mutations' that are the same as the ref base
                        mut_nts = copy.deepcopy(nts)
                        mut_nts[j] = mut
                        mut_codon = "".join(mut_nts)
                        if mut_codon in stop_codons:
                            mut_index = (i*3) + j
                            mutation_index_list[name][nt][mut].append(mut_index)

    # now also get the exon lengths
    exons = gen.read_many_fields(full_bed, "\t")
    exon_list_lengths = collections.defaultdict(lambda: collections.defaultdict())
    for exon in exons:
        start = int(exon[1])
        stop = int(exon[2])
        transcript_id = exon[3].split('.')[0]
        exon_id = int(exon[3].split('.')[1])
        type = exon[4]
        if type == "CDS":
            exon_list_lengths[transcript_id][exon_id] = stop-start

    # get all the possible cds indicies but grouped into exons
    exon_list_indices = collections.defaultdict(lambda: collections.defaultdict())
    for transcript in exon_list_lengths:
        current_end = 0
        current_start = 0
        for i, exon in enumerate(sorted(exon_list_lengths[transcript])):
            # get the length of the exon, add to the current length
            exon_length = exon_list_lengths[transcript][exon]
            current_end += exon_length
            # generate the indicies for that exon
            indices = list(range(current_start, current_end))
            exon_list_indices[transcript][exon] = [indices, current_start]
            # now move the start
            current_start += exon_length

    mutation_matched_indices = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    # for each of the possible ptc generating mutations
    for transcript in mutation_index_list:
        for ra in mutation_index_list[transcript]:
            for ma in mutation_index_list[transcript][ra]:
                # get the list of cds indices
                cds_index_list = mutation_index_list[transcript][ra][ma]
                for cds_index in cds_index_list:
                    for exon in exon_list_indices[transcript]:
                        if cds_index in exon_list_indices[transcript][exon][0]:
                            start_index = exon_list_indices[transcript][exon][1]
                            mutation_matched_indices[transcript][ra][ma].append("eidx{0}:cidx{1}:eid{2}".format(cds_index-start_index, cds_index, exon))

    nts = ["A", "C", "G", "T"]
    possible_mutation_files = {}
    for nt in nts:
        output_file = "{0}/possible_ptc_locations_{1}.fasta".format(output_directory, nt)
        with open(output_file, "w") as outfile:
            for transcript in mutation_matched_indices:
                for ma in sorted(mutation_matched_indices[transcript][nt]):
                    outfile.write(">{0}:{1}-{2}\n".format(transcript, nt, ma))
                    outfile.write("{0}\n".format(",".join(mutation_matched_indices[transcript][nt][ma])))


def generate_pseudo_snps(snp_file, possible_locations, exon_list, output_file):
    snps = gen.read_many_fields(snp_file, "\t")

    with open(output_file, "w") as outfile:
        np.random.seed()
        for snp in snps[1:]:
            snp_info = EntryInfo(snp)
            if len(possible_locations[snp_info.transcript_id][snp_info.aa][snp_info.ma]):
                possible_ptcs = possible_locations[snp_info.transcript_id][snp_info.aa][snp_info.ma][0].split(',')
                choices = list(range(len(possible_ptcs)))
                choice = np.random.choice(choices)
                pseudo_ptc_choice = possible_ptcs[choice]
                exon_index = int(re.findall('\d+', pseudo_ptc_choice.split(':')[0])[0])
                cds_index = re.findall('\d+', pseudo_ptc_choice.split(':')[1])[0]
                exon_number = re.findall('\d+', pseudo_ptc_choice.split(':')[2])[0]

                id = "{0}.{1}".format(snp_info.transcript_id, exon_number)
                exon_length = exon_list[id]
                min_dist = min(exon_index, int(exon_length) - int(exon_index))

                snp.extend([id, "{0}".format(exon_index), "{0}".format(cds_index), str(min_dist)])
                outfile.write("{0}\n".format("\t".join(snp)))

def get_coding_exons(coding_exons_fasta):

    coding_exons = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    exon_names, exon_seqs = gen.read_fasta(coding_exons_fasta)

    for i, name in enumerate(exon_names):
        transcript = name.split('(')[0].split('.')[0]
        exon = int(name.split('(')[0].split('.')[1])
        seq = exon_seqs[i]
        coding_exons[transcript][exon] = seq

    return coding_exons


def get_coding_exons_indices_no_ptcs(coding_exons, relative_positions_list, exclude_cpg=None):
    '''
    Get the index of each nt in each coding exon that isnt in the ptc list
    '''

    # get the locations of all ptcs grouped by transcript, exon, nt
    relative_ptc_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for ptc in relative_positions_list:
        transcript = relative_positions_list[ptc][0]
        exon_id = relative_positions_list[ptc][1]
        rel_pos = relative_positions_list[ptc][2]
        strand = relative_positions_list[ptc][3]
        ref_allele = relative_positions_list[ptc][4]
        if strand == "-":
            ref_allele = gen.reverse_complement(ref_allele)
        relative_ptc_locations[transcript][exon_id][ref_allele].append(rel_pos)

    # now get a list of all positions in an exon without a ptc
    coding_exon_nt_positions_no_ptc = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for transcript in coding_exons:
        for exon in coding_exons[transcript]:
            exon_seq = list(coding_exons[transcript][exon])
            for i, nt in enumerate(exon_seq):
                if i not in relative_ptc_locations[transcript][exon][nt]:
                    # if we want to exclude cpg regions
                    if exclude_cpg and i < len(exon_seq)-1:
                        if nt == "G" and exon_seq[i+1] != "C" or nt == "C" and exon_seq[i+1] != "G":
                            coding_exon_nt_positions_no_ptc[transcript][exon][nt].append(i)
                    else:
                        coding_exon_nt_positions_no_ptc[transcript][exon][nt].append(i)

    return coding_exon_nt_positions_no_ptc


def get_eses_from_file(ese_file):
    '''
    Get a list of eses from a file
    '''
    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t") if ese[0][0] != "#"]
    return ese_list

def get_ese_hits(ptc_list, coding_exons, ese_list):
    '''
    Get the number of ese hits in a given set
    '''
    hit_count = 0

    for ptc_id in ptc_list:
        transcript = ptc_list[ptc_id][0]
        exon_id = ptc_list[ptc_id][1]
        rel_pos = ptc_list[ptc_id][2]
        exon_seq = coding_exons[transcript][exon_id]
        possible_eses = get_possible_eses(exon_seq, rel_pos)
        overlaps = list(set(ese_list) & set(possible_eses))
        if len(overlaps):
            hit_count += 1
    return hit_count

def get_exon_list(file):
    '''
    Get a list of exons
    '''
    exons = gen.read_many_fields(file, "\t")
    exon_list = {}
    for exon in exons:
        exon_info = EntryInfo(exon)
        exon_list[exon_info.snp_id] = exon
    return(exon_list)

def get_exon_nts(exon_list, coding_exons):
    '''
    Get the count of nts in each window
    '''
    nts = [0,0,0]
    for id in exon_list:
        transcript = exon_list[id][0]
        exon_id = exon_list[id][1]
        exon_seq = coding_exons[transcript][exon_id]
        exon_length = len(exon_seq)

        if exon_length <= 4:
            nts[0] += 4
        elif exon_length > 4 and exon_length <= 138:
            nts[0] += 4
            nts[1] += 134
        else:
            nts[0] += 4
            nts[1] += 134
            nts[2] += (exon_length - 138)

    return nts

def get_introns(exon_bed, output_file):
    '''
    Get the introns between exons in a file
    '''

    exons = gen.read_many_fields(exon_bed, "\t")

    # get a dictionary of exons split by transcript and number
    exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for item in exons:
        exon = Define_Exon(item)
        if exon.type == "stop_codon":
            exon.exon_no = 999999
        exon_list[exon.transcript_id][exon.exon_no] = item

    # now get the introns and write to file
    with open(output_file, "w") as outfile:
        for transcript in exon_list:
            for exon_no in sorted(exon_list[transcript]):
                exon = Define_Exon(exon_list[transcript][exon_no])
                # check that the next exon exists, assuming its id will not be higher than 999999
                if exon.exon_no + 1 in exon_list[transcript]:
                    next_exon = Define_Exon(exon_list[exon.transcript_id][exon.exon_no + 1])

                    if exon.strand == "-":
                        intron_start = next_exon.stop
                        intron_stop = exon.start
                    else:
                        intron_start = exon.stop
                        intron_stop = next_exon.start

                    outlist = gen.stringify([exon.chr, intron_start, intron_stop, "{0}.{1}-{2}".format(exon.transcript_id, exon.exon_no, next_exon.exon_no), ".", exon.strand])
                    outfile.write("{0}\n".format("\t".join(outlist)))

def get_long_exons(relative_positions, coding_exons, length):
    '''
    Get only exons with length greater than defined length
    '''
    long_exons = {}
    for ptc_id in relative_positions:
        transcript = relative_positions[ptc_id][0]
        exon_id = relative_positions[ptc_id][1]
        exon_seq = coding_exons[transcript][exon_id]
        exon_length = len(exon_seq)

        # only consider exons with length > 138
        if exon_length > length:
            long_exons[ptc_id] = relative_positions[ptc_id]
    return long_exons

def get_min_dist(rel_pos, exon_length):
    '''
    Get the distance to the end of the exon
    rel_pos is in base 1
    '''

    # five prime dist is the relative position minus 1,
    # because a snp in the first nucleotide (rel_pos = 1)
    # has a distance of 0 to the end
    five_prime_dist = rel_pos - 1
    three_prime_dist = exon_length - rel_pos
    min_dist = min(five_prime_dist, three_prime_dist)

    return min_dist

def get_non_mut_window_nts(ptc_list, coding_exons, window_start, window_end):
    '''
    Get the incdices of non nonsense mutants in the given window
    '''

    # get the locations of all ptcs grouped by transcript, exon, nt
    relative_ptc_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for ptc in ptc_list:
        transcript = ptc_list[ptc][0]
        exon_id = ptc_list[ptc][1]
        rel_pos = ptc_list[ptc][2]
        strand = ptc_list[ptc][3]
        ref_allele = ptc_list[ptc][4]
        if strand == "-":
            ref_allele = gen.reverse_complement(ref_allele)
        relative_ptc_locations[transcript][exon_id][ref_allele].append(rel_pos)

    # now get a list of all positions in an exon without a ptc
    coding_exon_nt_positions_no_ptc = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for ptc in ptc_list:
        transcript = ptc_list[ptc][0]
        exon_id = ptc_list[ptc][1]
        strand = ptc_list[ptc][3]
        ref_allele = ptc_list[ptc][4]
        if strand == "-":
            ref_allele = gen.reverse_complement(ref_allele)

        exon_seq = list(coding_exons[transcript][exon_id])
        interested_indices = []
        five_prime_start = window_start - 1  # have to account for it being base 1
        five_prime_end = window_end
        interested_indices.extend(list(range(five_prime_start, five_prime_end)))
        three_prime_start = len(exon_seq) - window_end - 1
        three_prime_end = len(exon_seq) - window_start
        interested_indices.extend(list(range(three_prime_start, three_prime_end)))

        for id in interested_indices:
            if id not in relative_ptc_locations[transcript][exon_id][ref_allele]:
                coding_exon_nt_positions_no_ptc[transcript][exon_id][ref_allele].append(id)

    return coding_exon_nt_positions_no_ptc


def get_possible_eses(exon_seq, rel_pos):
    possible_eses = []
    if rel_pos - 5 >= 0:
        possible_eses.append(exon_seq[rel_pos-5:rel_pos+1])
    if rel_pos - 4 >= 0 and rel_pos + 1 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-4:rel_pos+2])
    if rel_pos - 3 >= 0 and rel_pos + 2 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-3:rel_pos+3])
    if rel_pos - 2 >= 0 and rel_pos + 3 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-2:rel_pos+4])
    if rel_pos - 1 >= 0 and rel_pos + 4 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-1:rel_pos+5])
    if rel_pos >= 0 and rel_pos + 5 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos:rel_pos+6])

    return possible_eses

def get_ptcs_in_window(relative_positions_list, window_start, window_end, coding_exons):
    '''
    Filter PTCs to only those in the given window (inclusive).
    window_start: start of the window (in base 1)
    window_end: end of the window (in base 1)
    '''

    ptcs_in_window = {}

    for ptc_id in relative_positions_list:
        transcript = relative_positions_list[ptc_id][0]
        exon_id = relative_positions_list[ptc_id][1]
        rel_pos = relative_positions_list[ptc_id][2]
        rel_pos_base_1 = rel_pos + 1

        exon_seq = coding_exons[transcript][exon_id]
        exon_length = len(exon_seq)

        min_dist = get_min_dist(rel_pos_base_1, exon_length)

        # need to account for the fact that min_dist will be 0 for ptc with relative
        # position 1 (in base 1) or 0
        # i.e. position 3, will have min_dist of 2
        if min_dist >= window_start - 1 and min_dist <= window_end -1:
            ptcs_in_window[ptc_id] = relative_positions_list[ptc_id]

    return ptcs_in_window


def get_ptc_ese_overlap(relative_position_list, coding_exons, ese_list):
    '''
    Get the number of snps in each region that overlap with
    something that looks like an ese
    '''

    ese_overlaps = [0,0,0]

    for ptc_id in relative_position_list:
        transcript = relative_position_list[ptc_id][0]
        exon_id = relative_position_list[ptc_id][1]
        rel_pos = relative_position_list[ptc_id][2]
        exon_seq = coding_exons[transcript][exon_id]
        exon_length = len(exon_seq)

        # need to get the location in base 1
        rel_pos_base_1 = rel_pos + 1
        # get the distance to the exon end
        min_dist = get_min_dist(rel_pos_base_1, exon_length)

        # get the possible ese sequences, use the relative position in base
        # 0 here becuase we want to extract sequences
        possible_eses = get_possible_eses(exon_seq, rel_pos)
        ese_overlap = list(set(ese_list) & set(possible_eses))

        # relative positions are in base 0
        # 0-2 base pairs correspond to min_dist of 0 or 1
        if min_dist < 2 and len(ese_overlap):
            ese_overlaps[0] += 1
        # 3-69 base pairs correspond to min_dist 2-68
        elif min_dist >= 2 and min_dist <= 68 and len(ese_overlap):
            ese_overlaps[1] += 1
        elif len(ese_overlap):
            ese_overlaps[2] += 1

    return ese_overlaps


def get_ptc_overlaps(disease_ptc_file, ptc_file, output_file):
    '''
    Get the overlap between disease ptc file and 1000 genomes ptcs
    '''

    disease_ptcs = gen.read_many_fields(disease_ptc_file, "\t")
    ptcs = gen.read_many_fields(ptc_file, "\t")

    disease_ptc_list = {}
    for ptc in disease_ptcs[1:]:
        ptc_info = EntryInfo(ptc)
        disease_ptc_list[ptc_info.snp_pos] = ptc

    ptc_list = {}
    for ptc in ptcs[1:]:
        ptc_list[int(ptc[7])] = ptc

    with open(output_file, "w") as outfile:
        for ptc_pos in ptc_list:
            if ptc_pos in disease_ptc_list:
                outfile.write('{0}\t**\t{1}\n'.format("\t".join(ptc_list[ptc_pos]), "\t".join(disease_ptc_list[ptc_pos])))




# def compare_distances(clinvar_ptc_file, clinvar_relative_exon_positions_file, ptc_file, ptc_relative_exon_positions_file, output_file):
#
#
#     clinvar_relative_position_entries = gen.read_many_fields(clinvar_relative_exon_positions_file, "\t")
#     clinvar_relative_positions = {}
#     for entry in clinvar_relative_position_entries:
#         # print(entry)
#         entry = entry[:12]
#         snp_id = entry[8]
#         exon_length = int(entry[2]) - int(entry[1])
#         rel_pos = entry[-1]
#         pos = int(entry[7])
#         clinvar_relative_positions[snp_id] = [snp_id, entry[3], int(entry[7]), int(rel_pos), exon_length]
#
#     # print(sorted([pos for pos in clinvar_relative_positions]))
#
#     ptc_relative_position_entries = gen.read_many_fields(ptc_relative_exon_positions_file, "\t")
#     ptc_relative_positions = {}
#     for entry in ptc_relative_position_entries:
#         entry = entry[:12]
#         snp_id = entry[8]
#         rel_pos = entry[-1]
#         pos = int(entry[7])
#         exon_length = int(entry[2]) - int(entry[1])
#         ptc_relative_positions[pos] = [snp_id, entry[3], int(entry[7]), int(rel_pos), exon_length]
#
#     clinvar_ptcs = gen.read_many_fields(clinvar_ptc_file, "\t")
#     clinvar_ptc_locations = []
#     clinvar_ptc_list = {}
#     for ptc in clinvar_ptcs:
#         pos = int(ptc[7])
#         clinvar_ptc_locations.append(pos)
#         clinvar_ptc_list[pos] = ptc
#     ptcs = gen.read_many_fields(ptc_file, "\t")
#     ptc_locations = []
#     ptc_list = {}
#     for ptc in ptcs[1:]:
#         pos = int(ptc[7])
#         ptc_locations.append(pos)
#         ptc_list[pos] = ptc
#
#     overlap = [loc for loc in ptc_locations if loc in clinvar_ptc_locations]
#
#     count = 0
#     with open(output_file, "w") as outfile:
#         outfile.write("transcript,exon,exon_length,rel_exon_pos,min_dist,exon_end,1000_genomes,clinvar\n")
#         for pos in ptc_list:
#             if pos not in overlap:
#                 ptc = ptc_list[pos]
#                 ptc_id = ptc[8]
#                 exon = ptc[3]
#                 t = exon.split('.')[0]
#                 e = exon.split('.')[1]
#                 e_start = int(ptc[1])
#                 e_end = int(ptc[2])
#                 gpos = int(ptc[7])
#                 rel_pos_info = ptc_relative_positions[gpos]
#                 rel_pos = rel_pos_info[3]
#                 exon_length = rel_pos_info[4]
#                 clinvar = 0
#                 kgenomes = 1
#
#                 dist_to_end = exon_length - rel_pos
#
#                 distances = [rel_pos, dist_to_end]
#                 min_dist = min(distances)
#
#                 if distances.index(min_dist) == 0:
#                     exon_end = 5
#                 else:
#                     exon_end = 3
#
#                 outlist = gen.stringify([t, e, exon_length, rel_pos, min_dist, exon_end, kgenomes, clinvar])
#                 outfile.write("{0}\n".format(",".join(outlist)))
#         for pos in clinvar_ptc_list:
#             if pos not in overlap:
#                 ptc = clinvar_ptc_list[pos]
#                 snp_id = ptc[8]
#                 if snp_id in clinvar_relative_positions:
#
#                     rel_pos_info = clinvar_relative_positions[snp_id]
#
#                     ptc_id = ptc[8]
#                     exon = ptc[3]
#                     t = exon.split('.')[0]
#                     e = exon.split('.')[1]
#                     e_start = int(ptc[1])
#                     e_end = int(ptc[2])
#                     gpos = int(ptc[7])
#
#                     rel_pos = rel_pos_info[3]
#                     exon_length = rel_pos_info[4]
#                     clinvar = 1
#                     kgenomes = 0
#
#
#
#                     dist_to_end = exon_length - rel_pos
#
#                     distances = [rel_pos, dist_to_end]
#                     min_dist = min(distances)
#
#                     if distances.index(min_dist) == 0:
#                         exon_end = 5
#                     else:
#                         exon_end = 3
#
#                     outlist = gen.stringify([t, e, exon_length, rel_pos, min_dist, exon_end, kgenomes, clinvar])
#                     outfile.write("{0}\n".format(",".join(outlist)))





def get_ptc_positions(relative_positions_list, coding_exons):

    '''
    Get the positions of nonsense mutations in the real dataset
    '''

    location_counts = [0,0,0]

    for ptc_id in relative_positions_list:
        transcript = relative_positions_list[ptc_id][0]
        exon_id = relative_positions_list[ptc_id][1]
        rel_pos = relative_positions_list[ptc_id][2]
        exon_seq = coding_exons[transcript][exon_id]
        exon_length = len(exon_seq)
        rel_pos = rel_pos + 1   # because the relative position is in base 0

        # have to minus 1 from rel pos because the distance from
        # nt 1 and the start is 0
        min_dist = min(rel_pos-1, exon_length - rel_pos)

        # 3-69bp region is between min_dist = 2 and min_dist = 68
        # becuase a mutation in the first position has a min_dist = 0,
        # mutation in second position has min_dist = 1 etc
        if min_dist < 2:
            location_counts[0] += 1
        elif min_dist >= 2 and min_dist <= 68:
            location_counts[1] += 1
        else:
            location_counts[2] += 1

    return location_counts

def get_real_ese_overlap(clean_ptc_list, relative_positions_list, regions, coding_exons, ese_list):
    '''
    Get the number of real eses overlaps
    '''
    ese_overlaps = {}
    for region in regions:
        ese_overlaps[region] = 0

    for ptc_id in clean_ptc_list:
        ptc = clean_ptc_list[ptc_id]
        rel_pos_info = relative_positions_list[ptc_id]
        transcript = ptc[4]
        exon = ptc[5]
        exon_seq = coding_exons[transcript][exon]
        rel_pos = rel_pos_info[1]
        exon_length = len(exon_seq)

        rel_pos_base_1 = rel_pos + 1

        distances = [rel_pos_base_1, exon_length - rel_pos_base_1]
        min_dist = min(distances)

        possible_eses = get_possible_eses(exon_seq, rel_pos)
        ese_overlap = list(set(ese_list) & set(possible_eses))

        # relative positions are in base 0
        if min_dist < 2 and len(ese_overlap):
            ese_overlaps[regions[0]] += 1
        elif min_dist >= 2 and min_dist <= 68 and len(ese_overlap):
            ese_overlaps[regions[1]] += 1
        elif len(ese_overlap):
            ese_overlaps[regions[2]] += 1

    return ese_overlaps

def get_ref_allele_matched_simulants(relative_positions_list, coding_exons_indices_no_ptcs):
    '''
    Get a set of random positions in exon with same matched ref allele
    '''
    simulant_list = {}
    cant_count = 0

    for ptc in relative_positions_list:
        transcript = relative_positions_list[ptc][0]
        exon_id = relative_positions_list[ptc][1]
        rel_pos = relative_positions_list[ptc][2]
        strand = relative_positions_list[ptc][3]
        ref_allele = relative_positions_list[ptc][4]
        mut_allele = relative_positions_list[ptc][5]
        if strand == "-":
            ref_allele = gen.reverse_complement(ref_allele)

        # get a list of positions to choose from
        possible_positions = coding_exons_indices_no_ptcs[transcript][exon_id][ref_allele]

        # if a position exists, pick a random one
        if len(possible_positions):
            simulant_position = np.random.choice(possible_positions, 1)[0]
        # otherwise, use relative position
        else:
            cant_count += 1
            simulant_position = rel_pos

        simulant_list[ptc] = [transcript, exon_id, simulant_position, strand, ref_allele, mut_allele]

    return simulant_list, cant_count

def get_relative_position_list(relative_position_file):
    '''
    Return relative positions of file, with ptc positions as key
    '''
    relative_positions = gen.read_many_fields(relative_position_file, "\t")
    relative_positions_list = {}
    for ptc in relative_positions:
        ptc_pos = ptc[8]
        exon = ptc[3]
        transcript = exon.split('.')[0]
        exon_id = int(exon.split('.')[1])
        rel_pos = int(ptc[11])
        strand = ptc[5]
        ref_allele = ptc[9]
        mut_allele = ptc[10]
        relative_positions_list[ptc_pos] = [transcript, exon_id, rel_pos, strand, ref_allele, mut_allele]

    return relative_positions_list


def get_unique_ptcs(disease_ptc_file, ptc_file, ptc_intersect_file, disease_output_file, kgenomes_output_file):
    '''
    Get unique disease ptcs that arent in the 1000 genomes set
    '''
    disease_ptcs = gen.read_many_fields(disease_ptc_file, "\t")
    ptcs = gen.read_many_fields(ptc_file, "\t")
    intersect_ptcs = gen.read_many_fields(ptc_intersect_file, "\t")

    intersect_list = []
    for ptc in intersect_ptcs:
        intersect_list.append(int(ptc[1]))

    with open(disease_output_file, "w") as outfile:
        for ptc in disease_ptcs:
            ptc_pos = int(ptc[7])
            if ptc_pos not in intersect_list:
                outfile.write('{0}\n'.format("\t".join(ptc)))

    with open(kgenomes_output_file, "w") as outfile:
        for ptc in ptcs[1:]:
            ptc_pos = int(ptc[7])
            if ptc_pos not in intersect_list:
                outfile.write('{0}\n'.format("\t".join(ptc)))

# def write_to_file(real_positions, process_list, regions, output_file):
#
#     # write all to file
#     header_list = ["sim_id"]
#     for region in regions:
#         header_list.append("{0}".format(region))
#
#     with open(output_file, "w") as outfile:
#         outfile.write("{0}\n".format(",".join(header_list)))
#         outlist = ["real"]
#         # for region in regions:
#         #     outlist.append(sum(real_positions[region].values()))
#         # for end in ends:
#         for region in regions:
#             outlist.append(real_positions[region])
#         outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))
#
#         for simulation in sorted(process_list):
#             outlist = [simulation]
#             # for region in regions:
#             #     outlist.append(sum(process_list[simulation][region].values()))
#             # for end in ends:
#             for region in regions:
#                 outlist.append(process_list[simulation][region])
#             outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))
# def get_clean_ptc_list(ptc_file, disease_ptcs_file):
#     '''
#     Get a list of PTCs from the dataset
#     '''
#     ptcs = gen.read_many_fields(ptc_file, "\t")
#     ptc_list = {}
#     ptc_positions = []
#     for ptc in ptcs[1:]:
#         ptc_id = ptc[8]
#         ptc_pos = int(ptc[7])
#         ptc_positions.append(ptc_pos)
#         start = int(ptc[1])
#         end = int(ptc[2])
#         strand = ptc[5]
#         ref_allele = ptc[9]
#         name = ptc[3]
#         t = ptc[3].split('.')[0]
#         e = int(ptc[3].split('.')[1])
#         ptc_list[ptc_id] = [name, start, end, ptc_pos, t, e, strand, ref_allele]
#
#     disease_ptcs = gen.read_many_fields(disease_ptcs_file, "\t")
#     disease_ptc_list = {}
#     disease_ptc_positions = []
#     for ptc in disease_ptcs:
#         ptc_id = ptc[8]
#         ptc_pos = int(ptc[7])
#         disease_ptc_positions.append(ptc_pos)
#         start = int(ptc[1])
#         end = int(ptc[2])
#         strand = ptc[5]
#         ref_allele = ptc[9]
#         name = ptc[3]
#         t = ptc[3].split('.')[0]
#         e = int(ptc[3].split('.')[1])
#         disease_ptc_list[ptc_id] = [name, start, end, ptc_pos, t, e, strand, ref_allele]
#
#     overlaps = [i for i in ptc_positions if i in disease_ptc_positions]
#
#     clean_genomes_ptc_list = {}
#     clean_clinvar_ptc_list = {}
#
#
#     for ptc_id in disease_ptc_list:
#         if disease_ptc_list[ptc_id][3] not in overlaps:
#             clean_clinvar_ptc_list[ptc_id] = disease_ptc_list[ptc_id]
#
#     for ptc_id in ptc_list:
#         if ptc_list[ptc_id][3] not in overlaps:
#             clean_genomes_ptc_list[ptc_id] = ptc_list[ptc_id]
#
#     return clean_genomes_ptc_list, clean_clinvar_ptc_list



def get_unique_rel_pos(unique_ptcs, disease_snps_relative_exon_positions, kgenomes_ptcs_file, kgenomes_ptcs_exon_positions, unique_ptcs_rel_pos_file, kgenomes_ptcs_rel_pos_file):
    '''
    Get the relative positions of the unique ptcs
    '''
    snps = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
    snp_list = collections.defaultdict(lambda: collections.defaultdict())
    for snp in snps:
        snp_pos = int(snp[7])
        rel_pos = int(snp[11])
        snp_list[snp_pos] = rel_pos

    ptcs = gen.read_many_fields(unique_ptcs, "\t")
    with open(unique_ptcs_rel_pos_file, "w") as outfile:
        for ptc in ptcs:
            ptc_pos = int(ptc[7])
            ptc[11] = snp_list[ptc_pos]
            outfile.write("{0}\n".format("\t".join(gen.stringify(ptc))))

    kgenomes_ptc_positions = gen.read_many_fields(kgenomes_ptcs_exon_positions, "\t")
    kgenomes_ptc_list = collections.defaultdict(lambda: collections.defaultdict())
    for ptc in kgenomes_ptc_positions[1:]:
        snp_pos = int(ptc[7])
        rel_pos = int(ptc[11])
        kgenomes_ptc_list[snp_pos] = rel_pos

    kgenomes_ptcs = gen.read_many_fields(kgenomes_ptcs_file, "\t")
    with open(kgenomes_ptcs_rel_pos_file, "w") as outfile:
        for ptc in kgenomes_ptcs:
            ptc_pos = int(ptc[7])
            ptc[11] = kgenomes_ptc_list[ptc_pos]
            outfile.write("{0}\n".format("\t".join(gen.stringify(ptc))))



def ptc_location_simulation(rel_pos_file, coding_exons_fasta, simulations, output_file, ese_overlap_output_file, ese_file=None, only_ese=None, exclude_cpg=None):
    '''
    Simulation mutation locations of PTCs.
    Take the exon in which each PTC is location and randomly pick a site with the
    same nt composition.
    Locations of these matched mutations are used for null.
    '''

    # get a list of the relative positions of the ptcs
    relative_positions_list = get_relative_position_list(rel_pos_file)
    # get a list of eses
    ese_list = get_eses_from_file(ese_file)
    # get the coding exons
    coding_exons = get_coding_exons(coding_exons_fasta)

    # get the positions of the ptcs
    real_positions = get_ptc_positions(relative_positions_list, coding_exons)
    # get the number of ptcs with ese overlaps
    real_ese_overlap = get_ptc_ese_overlap(relative_positions_list, coding_exons, ese_list)

    # now do the simulations
    simulant_list = list(range(1, simulations+1))
    processes = gen.run_in_parallel(simulant_list, ["foo", simulations, relative_positions_list, coding_exons_fasta, ese_list, exclude_cpg], simulate_mutation_locations)

    position_list = {}
    ese_overlap_list = {}
    for process in processes:
        result = process.get()
        position_list = {**position_list, **result[0]}
        ese_overlap_list = {**ese_overlap_list, **result[1]}

    # ignore writing this to file if we just want the ese overlap
    if not only_ese:
        with open(output_file, "w") as outfile:
            outfile.write('simulation,0.2,3.69,70+\n')
            outfile.write('real,{0}\n'.format(",".join(gen.stringify(real_positions))))
            for simulant in position_list:
                outfile.write("{0},{1}\n".format(simulant, ",".join(gen.stringify(position_list[simulant]))))

    with open(ese_overlap_output_file, "w") as outfile:
        outfile.write('simulation,0.2,3.69,70+\n')
        outfile.write('real,{0}\n'.format(",".join(gen.stringify(real_ese_overlap))))
        for simulant in ese_overlap_list:
            outfile.write("{0},{1}\n".format(simulant, ",".join(gen.stringify(ese_overlap_list[simulant]))))


# def ptc_location_simulation(snp_file, full_bed, cds_fasta, possible_positions_dir, output_directory, required_simulations, coding_exons_file):
    '''
    Simulate the snp location.
    For each snp, pick another site that has the same reference allele and that would generate a ptc with the mutated allele.
    Repeat n times.
    '''

    # return all the possible_locations
    possible_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    nts = ["A", "C", "G", "T"]
    for nt in nts:
        location_file = "{0}/possible_ptc_locations_{1}.fasta".format(possible_positions_dir, nt)
        entry_names, entry_locations = gen.read_fasta(location_file)
        for i, name in enumerate(entry_names):
            exon = name.split(':')[0]
            aa = name.split(':')[1][0]
            ma = name.split(':')[1][-1]
            possible_locations[exon][aa][ma].append(entry_locations[i])

    # get a list of exons and their lengths
    exons = gen.read_many_fields(coding_exons_file, "\t")
    exon_list = {}
    for exon in exons:
        exon_list[exon[3]] = int(exon[2]) - int(exon[1])

    # create a list of required simulations
    simulations = list(range(1, int(required_simulations) + 1))
    run_location_simulations(simulations, snp_file, possible_locations, exon_list, output_directory)

def refactor_ptc_file(input_file, output_file, header=None):
    '''
    Refactor PTC file ready for intersect
    '''
    entries = gen.read_many_fields(input_file, "\t")
    if header:
        start = 1
    else:
        start = 0

    with open(output_file, "w") as outfile:
        for entry in entries[start:]:
            exon_info = entry[:6]
            ptc_info = entry[6:]
            ptc_end = str(int(ptc_info[1]) + 1)
            strand = exon_info[5]
            ptc_info.insert(2, ptc_end)
            ptc_info.insert(3, ".")
            ptc_info.insert(4, strand)
            ptc_info.extend(exon_info)
            ptc_info[0] = ptc_info[0].strip('chr')
            outfile.write("{0}\n".format("\t".join(ptc_info)))

def run_location_simulations(simulations, snp_file, possible_locations, exon_list, output_directory):
    for simulation in simulations:
        print('{0}/{1}...'.format(simulation, len(simulations)))
        output_file = "{0}/disease_location_simulation_{1}.bed".format(output_directory, simulation)
        generate_pseudo_snps(snp_file, possible_locations, exon_list, output_file)


def simulate_ese_hits(simulation_list, simulations, ptc_list, coding_exons_fasta, ese_list, window_start, window_end):
    '''
    Simulations for ese hits in window
    '''

    hit_counts = {}

    # have to get the coding exons here during parallelisation
    coding_exons = get_coding_exons(coding_exons_fasta)

    # # get indices of each nt in that region that dont contain a ptc mutation
    window_ptcs_no_mut_indices = get_non_mut_window_nts(ptc_list, coding_exons, window_start, window_end)

    for simulation in simulation_list:
        np.random.seed()
        print("Simulation {0}/{1}".format(simulation+1, simulations))

        hit_count = 0
        cant_count = 0

        ref_allele_matched_simulants, cant_count = get_ref_allele_matched_simulants(ptc_list, window_ptcs_no_mut_indices)
        simulant_ese_hits = get_ese_hits(ref_allele_matched_simulants, coding_exons, ese_list)
        hit_counts[simulation] = [simulant_ese_hits, cant_count]

    return hit_counts


def simulate_mutation_locations(simulant_list, simulations, relative_positions_list, coding_exons_fasta, ese_list, exclude_cpg):
    '''
    Run the simulations for nonsense mutation locations
    '''

    # This needs to be done in here for parallelisation for some reason
    coding_exons = get_coding_exons(coding_exons_fasta)

    # get a list of all exon indices without a ptc
    coding_exons_indices_no_ptcs = get_coding_exons_indices_no_ptcs(coding_exons, relative_positions_list, exclude_cpg)

    location_outputs = {}
    ese_overlap_outputs = {}

    for simulation in simulant_list:
        np.random.seed()
        print("Running simulation {0}/{1}".format(simulation, simulations))

        ref_allele_matched_simulants, cant_counts = get_ref_allele_matched_simulants(relative_positions_list, coding_exons_indices_no_ptcs)
        simulant_positions = get_ptc_positions(ref_allele_matched_simulants, coding_exons)
        simulant_ese_overlaps = get_ptc_ese_overlap(ref_allele_matched_simulants, coding_exons, ese_list)
        location_outputs[simulation] = simulant_positions
        ese_overlap_outputs[simulation] = simulant_ese_overlaps

    return location_outputs, ese_overlap_outputs
