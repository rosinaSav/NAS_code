'''
Author: Liam Abrahams
Miscellaneous tests
'''


import disease_snps_ops as dso
import generic as gen
import os
import numpy as np
from scipy.stats import chisquare, ranksums
import collections
import itertools

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


def get_filtered_skipped_exons(final_output_file):
    lines = gen.read_many_fields(final_output_file, "\t")

    max_sample_count = max([float(line[2]) for line in lines[1:]])
    filtered_list = []

    # Exclude exons in which none of the individuals that we could quantify splicing contained a PTC,
    # and exons for which we could quantify splicing in less than half of samples
    for line in lines[1:]:
        sample_count = float(line[2])
        if float(line[1]) > 0 and sample_count > np.divide(max_sample_count, 2):
            line[3], line[4], line[5], line[6], line[7], line[8] = float(line[3])*100, float(line[4])*100, float(line[5])*100, float(line[6])*1000000, float(line[7])*1000000, float(line[8])*1000000
            filtered_list.append(line)

    return filtered_list

def get_large_effect_overlaps(filtered_list, psi_diff, rpmskip_diff):
    psi_ptc_large_skip = []
    for exon in filtered_list[:]:
        psi_het_ptc = exon[4]
        psi_no_ptc = exon[5]
        diff = psi_no_ptc - psi_het_ptc

        # we want those where psi is lower with ptc
        if diff > psi_diff:
            psi_ptc_large_skip.append(exon[0])

    rpmskim_ptc_large_skip = []
    for exon in filtered_list:
        rpmskip_het_ptc = exon[7]
        rpmskip_no_ptc = exon[8]
        diff = rpmskip_no_ptc - rpmskip_het_ptc

        # here we want those where rpmskip is higher with ptc
        if diff < -rpmskip_diff:
            rpmskim_ptc_large_skip.append(exon[0])


    overlap = [exon for exon in psi_ptc_large_skip if exon in rpmskim_ptc_large_skip]
    non_large = [exon[0] for exon in filtered_list if exon[0] not in overlap]

    return overlap, non_large

def get_relative_exon_positions(relative_positions_file, large_effect_exons, ptc_file):

    snp_pos = {}
    relative_positions = gen.read_many_fields(relative_positions_file, "\t")
    for snp in relative_positions:
        id = snp[8]
        pos = int(snp[11])
        snp_pos[id] = pos

    with open("results/clean_run_2/clean_run_PTC_relative_exon_positions.bed", "w") as outfile:
        ptcs = gen.read_many_fields(ptc_file, "\t")
        outfile.write('{0}'.format("\t".join(ptcs[0])))

        for ptc in ptcs[1:]:
            id = ptc[8]
            ptc[11] = snp_pos[id]
            outfile.write('{0}\n'.format("\t".join(gen.stringify(ptc))))

def list_from_fasta(file):

    fasta_list = {}
    names, seq = gen.read_fasta(file)
    for i, name in enumerate(names):
        name = name.split('(')[0]
        fasta_list[name] = seq[i]
    return fasta_list

def get_ese_overlap_count(ptcs, exon_seqs, ese_list, real=None):
    overlap_count = 0
    for ptc in ptcs:
        exon = ptc[3]
        rel_pos = int(ptc[11])
        possible_eses = dso.get_possible_eses(exon_seqs[exon], rel_pos)
        ese_overlap = list(set(ese_list) & set(possible_eses))
        if len(ese_overlap):
            overlap_count += 1

    return overlap_count

def simulate_exon_ese_overlaps(sims, large_effect_info, exon_seqs, ese_list, non_large_effect_info):
    '''
    Pick a random set of exons from the non-large effect group, then count how many
    hit an ese from the ese list
    '''
    choices = list(range(len(non_large_effect_info)))

    outputs = []
    for simulation in sims:
        np.random.seed()
        # choose the same number of non large effect cases as there are large effect cases
        random_choices = np.random.choice(choices, len(large_effect_info))
        random_exons = [non_large_effect_info[i] for i in random_choices]
        ese_overlap_count = get_ese_overlap_count(random_exons, exon_seqs, ese_list)
        outputs.append(ese_overlap_count)
    return outputs

def simulate_exon_ese_hits(simulations, relative_positions, large_effects, non_large_effects, exon_seqs, ese_list):

    # get the information on the ptc and add to the two lists
    large_effect_info = []
    non_large_effect_info = []
    for ptc in relative_positions:
        exon = ptc[3]
        if exon in large_effects:
            large_effect_info.append(ptc)
        elif exon in non_large_effects:
            non_large_effect_info.append(ptc)

    real_ese_overlap = get_ese_overlap_count(large_effect_info, exon_seqs, ese_list, real=True)

    sims = list(range(simulations))
    ese_overlaps = simulate_exon_ese_overlaps(sims, large_effect_info, exon_seqs, ese_list, non_large_effect_info)
    # print(ese_overlaps)
    processes = gen.run_in_parallel(sims, ["foo", large_effect_info, exon_seqs, ese_list, non_large_effect_info], simulate_exon_ese_overlaps)

    outputs = []
    for process in processes:
        output = process.get()
        outputs.extend(output)

    pval = np.divide(len([i for i in outputs if i >= real_ese_overlap]) + 1, simulations + 1)
    print("Number of PTCs that hit an ESE in the real cases: {0}/{1}".format(real_ese_overlap, len(large_effects)))
    print("Is this a significant number when picking non-large effect cases: {0}".format(pval))

def ese_hits_simulation(ese_file):
    output_prefix = "results/clean_run_2/clean_run"
    ptc_file = "{0}_ptc_file.txt".format(output_prefix)
    relative_positions_file = "{0}_PTC_relative_exon_positions.bed".format(output_prefix)
    coding_exons_file = "{0}_coding_exons.fasta".format(output_prefix)
    final_output_file = "{0}__analysis_final_output.txt".format(output_prefix)

    filtered_list = get_filtered_skipped_exons(final_output_file)
    large_effects, non_large_effects = get_large_effect_overlaps(filtered_list, 5, 0.025)

    relative_positions = gen.read_many_fields(relative_positions_file, "\t")
    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t")[1:]]
    exon_seqs = list_from_fasta(coding_exons_file)

    simulations = 10000
    simulate_exon_ese_hits(simulations, relative_positions, large_effects, non_large_effects, exon_seqs, ese_list)


def get_ptc_regions(ptc_list, rel_pos_list):

    region_counts = [0,0,0]
    ends = [0,0]

    for ptc in ptc_list:
        rel_pos = rel_pos_list[ptc][0]
        start = rel_pos_list[ptc][1]
        stop = rel_pos_list[ptc][2]
        exon_length = stop-start
        five_prime_dist = rel_pos
        three_prime_dist = exon_length - rel_pos
        min_dist = min(five_prime_dist, three_prime_dist)
        if min_dist == five_prime_dist:
            ends[0] += 1
        else:
            ends[1] += 1

        if min_dist <= 2:
            region_counts[0] += 1
        elif min_dist > 2 and min_dist <= 69:
            region_counts[1] += 1
        elif min_dist > 69:
            region_counts[2] += 1

    return region_counts, ends

def large_effects_locations_sim(sims, large_effects, non_large_effects, rel_pos_list):

    choices = list(range(len(non_large_effects)))
    outputs = []
    ends = []
    for simulation in sims:
        np.random.seed()
        # choose the same number of non large effect cases as there are large effect cases
        random_choices = np.random.choice(choices, len(large_effects))
        random_exons = [non_large_effects[i] for i in random_choices]
        sim_region_count, sim_ends = get_ptc_regions(random_exons, rel_pos_list)
        outputs.append(sim_region_count)
        ends.append(sim_ends)
    return outputs, ends

def large_effect_locations_sim():
    '''
    Test where the large effect cases are
    '''
    output_prefix = "results/clean_run_2/clean_run"
    ptc_file = "{0}_ptc_file.txt".format(output_prefix)
    relative_positions_file = "{0}_PTC_relative_exon_positions.bed".format(output_prefix)
    final_output_file = "{0}__analysis_final_output.txt".format(output_prefix)

    filtered_list = get_filtered_skipped_exons(final_output_file)
    large_effects, non_large_effects = get_large_effect_overlaps(filtered_list, 5, 0.025)

    relative_positions = gen.read_many_fields(relative_positions_file, "\t")

    rel_pos_list = {}
    for ptc in relative_positions[1:]:
        start = int(ptc[1])
        stop = int(ptc[2])
        exon = ptc[3]
        rel_pos = int(ptc[11])
        rel_pos_list[exon] = [rel_pos, start, stop]

    real_regions, real_ends = get_ptc_regions(large_effects, rel_pos_list)

    simulations = 10000
    sims = list(range(simulations))
    # outputs = large_effects_locations_sim(sims, large_effects, non_large_effects, rel_pos_list)
    processes = gen.run_in_parallel(sims, ["foo", large_effects, non_large_effects, rel_pos_list], large_effects_locations_sim)

    regions = []
    ends = []
    for process in processes:
        output = process.get()
        region = output[0]
        end = output[1]
        for i in region:
            regions.append(i)
        for i in end:
            ends.append(i)


    ese_region_pval = np.divide(len([1 for region in regions if region[1] >= real_regions[1]]) + 1, len(regions) + 1)
    ends_pval = np.divide(len([1 for end in ends if end[0] >= real_ends[0]]) + 1, len(ends) + 1)
    # ese_region_pval = 1

    print("PTCs in regions (0-2,3-69,70+): {0}".format(real_regions))
    print("Ends (5', 3'): {0}".format(real_ends))
    print("ESE region compared with sims: {0}".format(ese_region_pval))
    print("ESE exon ends with sims: {0}".format(ends_pval))

def get_exon_length_info(ptc_list, rel_pos_list):

    periods = [0,0,0]
    lengths = []

    for ptc in ptc_list:
        start = rel_pos_list[ptc][1]
        stop = rel_pos_list[ptc][2]
        exon_length = stop-start
        lengths.append(exon_length)
        periods[exon_length % 3] += 1
    return lengths, periods

def sim_exon_length_info(sims, large_effects, non_large_effects, rel_pos_list):

    lengths = []
    periods = []
    choices = list(range(len(non_large_effects)))
    for simulation in sims:
        np.random.seed()
        # choose the same number of non large effect cases as there are large effect cases
        random_choices = np.random.choice(choices, len(large_effects))
        random_exons = [non_large_effects[i] for i in random_choices]
        sim_lengths, sim_periods = get_exon_length_info(random_exons, rel_pos_list)
        lengths.append(sim_lengths)
        periods.append(sim_periods)
    return lengths, periods




def large_effects_lengths_sim():
    '''
    Test the lengths of the large effects exons
    and whether they are biased of length 3
    '''

    output_prefix = "results/clean_run_2/clean_run"
    ptc_file = "{0}_ptc_file.txt".format(output_prefix)
    relative_positions_file = "{0}_PTC_relative_exon_positions.bed".format(output_prefix)
    final_output_file = "{0}__analysis_final_output.txt".format(output_prefix)

    filtered_list = get_filtered_skipped_exons(final_output_file)
    large_effects, non_large_effects = get_large_effect_overlaps(filtered_list, 5, 0.025)

    relative_positions = gen.read_many_fields(relative_positions_file, "\t")

    rel_pos_list = {}
    for ptc in relative_positions[1:]:
        start = int(ptc[1])
        stop = int(ptc[2])
        exon = ptc[3]
        rel_pos = int(ptc[11])
        rel_pos_list[exon] = [rel_pos, start, stop]

    real_lengths, real_periods = get_exon_length_info(large_effects, rel_pos_list)

    simulations = 10000
    sims = list(range(simulations))
    # sim_exon_length_info(sims, large_effects, non_large_effects, rel_pos_list)
    processes = gen.run_in_parallel(sims, ["foo", large_effects, non_large_effects, rel_pos_list], sim_exon_length_info)

    lengths = []
    periods = []
    for process in processes:
        output = process.get()
        length = output[0]
        period = output[1]
        for i in length:
            lengths.append(i)
        for i in period:
            periods.append(i)

    length_pval = np.divide(len([1 for length in lengths if np.mean(length) <= np.mean(real_lengths)]) + 1, len(lengths) + 1)
    # for those exons of length 3
    period_pval = np.divide(len([1 for period in periods if period[0] >= real_periods[0]]) + 1, len(periods) + 1)
    # ese_region_pval = 1

    print("Mean large effect exon length: {0}".format(np.mean(real_lengths)))
    print("Real periodicity (0,1,2): {0}".format(real_periods))
    print("Lengths compared with sims: {0}".format(length_pval))
    print("Periodicity with sims: {0}".format(period_pval))

def get_filtered_exons():
    '''
    Create a ptc file containing only the filtered final ptcs
    and large effect cases that overlap disease
    '''
    output_prefix = "results/clean_run_2/clean_run"
    disease_prefix = "results/clinvar"
    ptc_file = "{0}_ptc_file.txt".format(output_prefix)
    relative_positions_file = "{0}_PTC_relative_exon_positions.bed".format(output_prefix)
    final_output_file = "{0}__analysis_final_output.txt".format(output_prefix)

    intersect_file = "{0}/ptc_intersect.bed".format(disease_prefix)
    intersect_ptcs = gen.read_many_fields(intersect_file, "\t")
    intersect_list = []
    for ptc in intersect_ptcs:
        exon = ptc[-4]
        intersect_list.append(exon)

    filtered_list = get_filtered_skipped_exons(final_output_file)
    large_effects, non_large_effects = get_large_effect_overlaps(filtered_list, 5, 0.025)
    filtered_list = [exon[0] for exon in get_filtered_skipped_exons(final_output_file)]

    ptcs = gen.read_many_fields(ptc_file, "\t")
    ptc_list = []
    le_overlap = []
    for ptc in ptcs:
        exon = ptc[3]
        if exon in filtered_list and exon not in intersect_list:
            ptc_list.append(ptc)
        if exon in large_effects and exon in intersect_list:
            le_overlap.append(ptc)

    large_effect_disease_overlap = "{0}_large_effect_disease_overlap.txt".format(output_prefix)
    with open(large_effect_disease_overlap, "w") as outfile:
        for ptc in le_overlap:
            exon = ptc[3]
            outfile.write('{0}\n'.format(exon))

    filtered_no_overlaps = "{0}_ptc_file_filtered_no_disease.txt".format(output_prefix)
    with open(filtered_no_overlaps, "w") as outfile:
        for ptc in ptc_list:
            outfile.write("{0}\n".format("\t".join(gen.stringify(ptc))))



def main():

    description = "Miscellaneous tests."
    args = gen.parse_arguments(description,  ["results_prefix", "disease_output_dir", "ese_file", "get_filtered", "get_info", "disease_locations_chisquare", "large_effect_ese_hits_simulation", "large_effect_locations", "large_effect_lengths"], flags = [3,4,5,6,7,8])
    results_prefix, disease_output_dir, ese_file, get_filtered, get_info, disease_locations_chisquare, large_effect_ese_hits_simulation, large_effect_locations, large_effect_lengths = args.results_prefix, args.disease_output_dir, args.ese_file, args.get_filtered, args.get_info, args.disease_locations_chisquare, args.large_effect_ese_hits_simulation, args.large_effect_locations, args.large_effect_lengths

    if get_filtered:
        get_filtered_exons()

    # tests on the large effect cases
    if large_effect_ese_hits_simulation:
        ese_hits_simulation(ese_file)
    if large_effect_locations:
        large_effect_locations_sim()
    if large_effect_lengths:
        large_effects_lengths_sim()


if __name__ == "__main__":
    main()
