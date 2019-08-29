import generic as gen
import collections
import re
import requests
import json
import numpy as np
import scipy.stats

ptc_file = "results/clean_run_2/clean_run_ptc_file.txt"
others_file = "results/clean_run_2/clean_run_syn_nonsyn_file.txt"
final_outputs = "results/clean_run_2/clean_run__analysis_final_output_filtered.txt"

final_ptc_file = "results/clean_run_2/clean_run__final_ptc_list.txt"
final_other_mutations = "results/clean_run_2/clean_run__final_others_same_exons_list.txt"

ptcs = gen.read_many_fields(ptc_file, "\t")[1:]
print(len(ptcs))
exons = [i[0] for i in gen.read_many_fields(final_outputs, "\t")[1:]]
exons = gen.read_many_fields(final_outputs, "\t")
others = gen.read_many_fields(others_file, "\t")[1:]
max_sample_count = max([float(i[2]) for i in exons])
exons = [i for i in exons if float(i[1]) != 0 and float(i[2]) > 0.5*max_sample_count]

exon_list = [i[0] for i in exons]

kept = []
kept_ptcs = []
for i in ptcs:
    if i[3] not in kept and i[3] in exon_list:
        kept.append(i[3])
        kept_ptcs.append(i)

with open(final_ptc_file, "w") as outfile:
    [outfile.write("{0}\n".format("\t".join(i))) for i in kept_ptcs]

exon_ids = [i[3] for i in ptc_list]

other_mutations = gen.read_many_fields(others_file, "\t")[1:]
other_mutations_list = [i for i in other_mutations if i[3] in exon_ids]

with open(final_other_mutations, "w") as outfile:
    [outfile.write("{0}\n".format("\t".join(i))) for i in other_mutations_list]


ptc_list = gen.read_many_fields(final_ptc_file, "\t")
other_mutations = gen.read_many_fields(final_other_mutations, "\t")

def calc_frequencies(set):
    # mutations = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    mutations = {"A":0,"C":0,"G":0,"T":0}
    for ptc in set:
        if ptc[5] == "-":
            ptc[10] = gen.reverse_complement(ptc[10])
            ptc[9] = gen.reverse_complement(ptc[9])
        mutations[ptc[10]] += 1
        # mutations[ptc[9]][ptc[10]] += 1
    return mutations

ptc_frequencies = calc_frequencies(ptc_list)
other_frequencies = calc_frequencies(other_mutations)




ptc_proportions = {i: np.divide(ptc_frequencies[i], sum(ptc_frequencies.values())) for i in ptc_frequencies}
other_proportions = {i: np.divide(other_frequencies[i], sum(other_frequencies.values())) for i in other_frequencies}

expected_frequencies = [other_proportions[i]*len(ptc_list) for i in sorted(other_proportions)]
observed_frequencies = [ptc_frequencies[i] for i in sorted(ptc_frequencies)]
chi = scipy.stats.chisquare(observed_frequencies, f_exp = expected_frequencies)


ptc_count_no_c = sum([ptc_frequencies[i] for i in ptc_frequencies if i != "C"])
other_count_no_c = sum([other_frequencies[i] for i in other_frequencies if i != "C"])


ptc_proportions_no_c = {i: np.divide(ptc_frequencies[i], ptc_count_no_c) for i in ptc_frequencies if i != "C"}
other_proportions_no_c = {i: np.divide(other_frequencies[i], other_count_no_c) for i in other_frequencies if i != "C"}


expected_frequencies_no_c = [other_proportions_no_c[i] * len(ptc_list) for i in sorted(other_proportions_no_c) if i != "C"]
observed_frequencies_no_c = [ptc_frequencies[i] for i in sorted(ptc_frequencies) if i != "C"]
chi_no_c = scipy.stats.chisquare(observed_frequencies_no_c, f_exp = expected_frequencies_no_c)




output_file = "results/clean_run_2/mutation_types_chisq_all.csv"
with open(output_file, "w") as outfile:
    outfile.write("nucleotide,non_ptc,prop_non_ptc,ptc,prop_ptc,,expected_ptc,o/e\n")
    for i, nt in enumerate(sorted(ptc_frequencies)):
        outfile.write("{0},{1},{2},{3},{4},,{5},{6}\n".format(nt, other_frequencies[nt], other_proportions[nt], ptc_frequencies[nt], ptc_proportions[nt], expected_frequencies[i], np.divide(ptc_frequencies[nt], expected_frequencies[i])))

    outfile.write(",{0},{1},{2},{3},,{4}\n".format(sum(other_frequencies.values()), sum(other_proportions.values()), sum(ptc_frequencies.values()), sum(ptc_proportions.values()), sum(expected_frequencies) ))

    outfile.write("\n")
    outfile.write("Chisq\n")
    outfile.write("Statistic:,{0}\n".format(chi.statistic))
    outfile.write("P value:,{0}\n".format(chi.pvalue))

output_file = "results/clean_run_2/mutation_types_chisq_no_c.csv"
with open(output_file, "w") as outfile:
    outfile.write("nucleotide,non_ptc,prop_non_ptc,ptc,prop_ptc,,expected_ptc,o/e\n")
    for i, nt in enumerate(sorted(ptc_proportions_no_c)):
        outfile.write("{0},{1},{2},{3},{4},,{5},{6}\n".format(nt, other_frequencies[nt], other_proportions_no_c[nt], ptc_frequencies[nt], ptc_proportions_no_c[nt], expected_frequencies_no_c[i], np.divide(observed_frequencies_no_c[i], expected_frequencies_no_c[i])))

    outfile.write(",{0},{1},{2},{3},,{4}\n".format(other_count_no_c, sum(other_proportions_no_c.values()), ptc_count_no_c, sum(ptc_proportions_no_c.values()), sum(expected_frequencies_no_c) ))

    outfile.write("\n")
    outfile.write("Chisq\n")
    outfile.write("Statistic:,{0}\n".format(chi_no_c.statistic))
    outfile.write("P value:,{0}\n".format(chi_no_c.pvalue))
