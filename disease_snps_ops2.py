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


def ptc_locations_test(coding_exons_fa, pathogenic_ptc_file, likely_pathogenic_ptc_file, output_file):

    coding_exon_names, coding_exons_seqs = gen.read_fasta(coding_exons_fa)
    coding_exons = {}
    for i, name in enumerate(coding_exon_names):
        name = name.split("(")[0]
        coding_exons[name] = coding_exons_seqs[i]

    pathogenic = gen.read_many_fields(pathogenic_ptc_file, "\t")
    likely_pathogenic = gen.read_many_fields(likely_pathogenic_ptc_file, "\t")

    # pathogenic = []
    # likely_pathogenic = []
    #
    # for ptc in ptcs:
    #     if ptc[-2] == "pathogenic":
    #         pathogenic.append(ptc)
    #     if ptc[-2] == "likely_pathogenic":
    #         likely_pathogenic.append(ptc)
    #
    # with open(output_file1, "w") as outfile:
    #     [outfile.write("{0}\n".format("\t".join(i))) for i in pathogenic]
    # with open(output_file2, "w") as outfile:
    #     [outfile.write("{0}\n".format("\t".join(i))) for i in likely_pathogenic]

    threshold = 138

    pathogenic_chi = region_chisquare(pathogenic, coding_exons)
    pathogenic_chi_restrict = region_chisquare(pathogenic, coding_exons, restrict_length = threshold)
    likely_pathogenic_chi = region_chisquare(likely_pathogenic, coding_exons)
    likely_pathogenic_chi_restrict = region_chisquare(likely_pathogenic, coding_exons, restrict_length = threshold)

    with open(output_file, "w") as outfile:
        write_chi_results(pathogenic_chi, outfile, title = "Tag: Pathogenic")
        write_chi_results(likely_pathogenic_chi, outfile, title = "Tag: Likely Pathogenic")

        write_chi_results(pathogenic_chi_restrict, outfile, title = "Tag: Pathogenic (exon length > {0})".format(threshold))
        write_chi_results(likely_pathogenic_chi_restrict, outfile, title = "Tag: Likely Pathogenic (exon length > {0})".format(threshold))



def write_chi_results(result, outfile, title = None):

    if title:
        outfile.write("{0}\n".format(title))

    nts, nt_props, observed, expected, chi = result
    observed_props = {i: np.divide(observed[i], sum(observed.values())) for i in observed}

    outfile.write("exon position,nucleotides,nucleotide proportion,,ptcs,ptc proportion,,expected ptcs,o/e\n")
    for region in ["splice", "flank", "core"]:
        outfile.write("{0},{1},{2},,{3},{4},,{5},{6}\n".format(region, nts[region], nt_props[region], observed[region], observed_props[region], expected[region], np.divide(observed[region], expected[region])))
    outfile.write(",{0},{1},,{2},{3},,{4}\n".format(sum(nts.values()), sum(nt_props.values()), sum(observed.values()), sum(observed_props.values()), sum(expected.values())))
    outfile.write("\n")
    outfile.write("chisquare\n")
    outfile.write("statistic:,{0}\n".format(chi.statistic))
    outfile.write("p value:,{0}\n".format(chi.pvalue))
    outfile.write("\n")
    outfile.write("\n")


def region_chisquare(ptc_list, exon_list, restrict_length = None):

    seen = []
    nts = collections.Counter()

    kept_count = 0

    groups = collections.defaultdict(lambda: [])

    for ptc in ptc_list:
        # need to add 1 here for 0 indexing
        rel_pos = int(ptc[11]) + 1
        exon_length = len(exon_list[ptc[3]])
        exon_id = ptc[3]

        if not restrict_length or restrict_length and exon_length > restrict_length:
            kept_count += 1
            if exon_id not in seen:
                seen.append(exon_id)
                if 4 <= exon_length:
                    nts["splice"] += 4
                if 4 < exon_length <= 138:
                    nts["flank"] += (exon_length - 4)
                if 138 < exon_length:
                    nts["flank"] += 134
                    nts["core"] += (exon_length - 138)

            if rel_pos <= 2 or exon_length - rel_pos <= 2:
                groups["splice"].append(ptc)
            elif 2 < rel_pos <= 69 or 2 < exon_length - rel_pos <= 69:
                groups["flank"].append(ptc)
            else:
                groups["core"].append(ptc)

    group_counts = {i: len(groups[i]) for i in groups}

    nt_props = {i: np.divide(nts[i], sum(nts.values())) for i in nts}

    expected = {i: nt_props[i]*kept_count for i in nt_props}

    observed_counts = [group_counts[i] for i in sorted(group_counts)]
    expected_counts = [expected[i] for i in sorted(expected)]

    chi = chisquare(observed_counts, expected_counts)
    return [nts, nt_props, group_counts, expected, chi]


def ptc_location_sim(coding_exons_fa, ptc_file, output_file):

    coding_exon_names, coding_exons_seqs = gen.read_fasta(coding_exons_fa)
    coding_exons = {}
    for i, name in enumerate(coding_exon_names):
        name = name.split("(")[0]
        coding_exons[name] = coding_exons_seqs[i]

    ptcs = gen.read_many_fields(ptc_file, "\t")

    pathogenic = []
    likely_pathogenic = []

    for ptc in ptcs:
        if ptc[-2] == "pathogenic":
            pathogenic.append(ptc)
        if ptc[-2] == "likely_pathogenic":
            likely_pathogenic.append(ptc)

    ptc_positions = collections.defaultdict(lambda: [])

    for ptc in pathogenic:
        exon_id = ptc[3]
        rel_pos = int(ptc[11])
        ptc_positions[exon_id].append(rel_pos)

    for ptc in likely_pathogenic:
        exon_id = ptc[3]
        rel_pos = int(ptc[11])
        ptc_positions[exon_id].append(rel_pos)

    choices = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    for id in coding_exons:
        exon_seq = coding_exons[id]
        exon_length = len(exon_seq)
        for i in range(exon_length):
            if i not in ptc_positions[id]:
                choices[id][exon_seq[i]].append(i)

    test_choices = {}
    for id in choices:
        test_choices[id] = {}
        for nt in choices[id]:
            test_choices[id][nt] = choices[id][nt]
    choices = test_choices


    real_groups = collections.defaultdict(lambda: [])
    for i, ptc in enumerate(pathogenic):
        exon_id = ptc[3]
        ref_allele = ptc[9]
        rel_pos = int(ptc[11])
        exon_seq = coding_exons[exon_id]
        exon_length = len(exon_seq)

        rel_pos = rel_pos + 1

        if rel_pos <= 2 or exon_length - rel_pos <= 2:
            real_groups["splice"].append(ptc)
        elif 2 < rel_pos <= 69 or 2 < exon_length - rel_pos <= 69:
            real_groups["flank"].append(ptc)
        else:
            real_groups["core"].append(ptc)

    for region in ["splice", "flank", "core"]:
        if region not in real_groups:
            real_groups[region] = 0



    real_groups = {i: len(real_groups[i]) for i in real_groups}
    # print(real_groups)

    simulant_list = list(range(1, 20+1))
    outputs = gen.run_simulation_function(simulant_list, [pathogenic, choices, coding_exons], sim_ptc_position, sim_run = False)



    with open(output_file, "w") as outfile:
        outfile.write("id,pathogenic_splice,pathogenic_flank,pathogenic_core\n")
        outfile.write("real,{0},{1},{2}\n".format(real_groups["splice"], real_groups["flank"], real_groups["core"]))
        for i in sorted(outputs):
            outfile.write("{0},{1},{2},{3}\n".format(i, outputs[i]["splice"], outputs[i]["flank"], outputs[i]["core"]))

def sim_ptc_position(iterations, ptc_list, choices, coding_exons):

    outputs = {}

    if len(iterations) > 0:
        np.random.seed()
        for iteration in iterations:
            print(iteration)
            groups = collections.defaultdict(lambda: [])
            seen = collections.defaultdict(lambda: [])
            #
            for i, ptc in enumerate(ptc_list):
                exon_id = ptc[3]
                ref_allele = ptc[9]
                rel_pos = int(ptc[11])

                exon_seq = coding_exons[exon_id]
                exon_length = len(exon_seq)

                sim_choices = [i for i in choices[exon_id][ref_allele] if i not in seen[exon_id]]
                choice = np.random.choice(sim_choices)
                seen[exon_id].append(choice)

                choice = choice + 1

                if choice <= 2 or exon_length - choice <= 2:
                    groups["splice"].append(ptc)
                elif 2 < choice <= 69 or 2 < exon_length - choice <= 69:
                    groups["flank"].append(ptc)
                else:
                    groups["core"].append(ptc)

            for region in ["splice", "flank", "core"]:
                if region not in groups:
                    groups[region] = 0

            groups = {i: len(groups[i]) for i in groups}
            outputs[iteration] = groups

    return outputs
