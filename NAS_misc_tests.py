import generic as gen
import os
import numpy as np
from scipy.stats import chisquare, ranksums
import collections


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


def simulate_ese_hits(simulations, exon_list, ese_hits):

    hit_counts = []

    for i, simulation in enumerate(simulations):
        print('Simulation {0}/{1}'.format(i+1, len(simulations)))

        np.random.seed()
        random_exons = np.random.choice(exon_list, 31)
        ese_hit_count = sum([ese_hits[exon] for exon in random_exons])
        hit_counts.append(ese_hit_count)

    return(hit_counts)




def check_ese_hits(results_prefix, ese_file):

    final_output_file = "{0}__analysis_final_output.txt".format(results_prefix)
    ptc_file = "{0}_ptc_file.txt".format(results_prefix)
    relative_positions_file = "{0}_SNP_relative_exon_position.bed".format(results_prefix)
    coding_exons_fasta = "{0}_coding_exons.fasta".format(results_prefix)
    # coding_exons_bed = "{0}_coding_exons.bed".format(results_prefix)

    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t")[1:]]

    filtered_list = get_filtered_skipped_exons(final_output_file)
    filtered_exon_list = [exon[0] for exon in filtered_list]

    overlap = get_large_effect_overlaps(filtered_list, 5, 0.025)

    ptcs = gen.read_many_fields(ptc_file, "\t")

    relative_positions = gen.read_many_fields(relative_positions_file, "\t")
    relative_positions_list = {}
    for rel in relative_positions:
        relative_positions_list[int(rel[7])] = int(rel[11])

    ptc_positions = {}
    for ptc in ptcs[1:]:
        exon = ptc[3]
        rel_pos = int(ptc[11])
        position = int(ptc[7])
        ptc_positions[exon] = relative_positions_list[position]

    coding_exon_names, coding_exon_seqs = gen.read_fasta(coding_exons_fasta)
    coding_exon_list = {}
    for i, name in enumerate(coding_exon_names):
        coding_exon_list[name.split('(')[0]] = coding_exon_seqs[i]

    ese_hits = {}
    exon_list = []
    for exon in ptc_positions:
        seq = coding_exon_list[exon]
        rel_pos = ptc_positions[exon]
        possible_eses = get_possible_eses(seq, rel_pos)
        ese_overlap = list(set(ese_list) & set(possible_eses))
        exon_list.append(exon)
        if len(ese_overlap):
            ese_hits[exon] = 1
        else:
            ese_hits[exon] = 0

    real_hit_count = 0
    for exon in overlap:
        real_hit_count += ese_hits[exon]


    simulation_list = list(range(1000))

    sim_exon_list = [exon for exon in exon_list if exon not in filtered_exon_list]

    hit_counts = []
    processes = gen.run_in_parallel(simulation_list, ["foo", sim_exon_list, ese_hits], simulate_ese_hits)
    for process in processes:
        hit_counts.extend(process.get())


    comparison = [x for x in hit_counts if x >= real_hit_count]
    print(np.divide(len(comparison), len(hit_counts)))



def get_large_effect_overlaps(filtered_list, psi_diff, rpmskip_diff):
    psi_ptc_large_skip = []
    for exon in filtered_list[:]:
        psi_het_ptc = float(exon[4])*100
        psi_no_ptc = float(exon[5])*100
        diff = psi_no_ptc - psi_het_ptc

        # we want those where psi is lower with ptc
        if diff > psi_diff:
            psi_ptc_large_skip.append(exon[0])

    rpmskim_ptc_large_skip = []
    for exon in filtered_list:
        rpmskip_het_ptc = float(exon[7])*1000000
        rpmskip_no_ptc = float(exon[8])*1000000
        diff = rpmskip_no_ptc - rpmskip_het_ptc

        # here we want those where rpmskip is higher with ptc
        if diff < -rpmskip_diff:
            rpmskim_ptc_large_skip.append(exon[0])

    overlap = [exon for exon in psi_ptc_large_skip if exon in rpmskim_ptc_large_skip]
    return overlap
#
# def check_distances(results_prefix):
#
#     final_output_file = "{0}__analysis_final_output.txt".format(results_prefix)
#     ptc_file = "{0}_ptc_file.txt".format(results_prefix)
#     relative_positions_file = "{0}_SNP_relative_exon_position.bed".format(results_prefix)
#     coding_exons_fasta = "{0}_coding_exons.fasta".format(results_prefix)
#
#     filtered_list = get_filtered_skipped_exons(final_output_file)
#     overlap = get_large_effect_overlaps(filtered_list, 5, 0.025)
#     ptcs = gen.read_many_fields(ptc_file, "\t")
#
#     relative_positions = gen.read_many_fields(relative_positions_file, "\t")
#     relative_positions_list = {}
#     for rel in relative_positions:
#         relative_positions_list[int(rel[7])] = int(rel[11])
#
#     ptc_positions = {}
#     for ptc in ptcs[1:]:
#         exon = ptc[3]
#         rel_pos = int(ptc[11])
#         position = int(ptc[7])
#         ptc_positions[exon] = relative_positions_list[position]
#
#     coding_exon_names, coding_exon_seqs = gen.read_fasta(coding_exons_fasta)
#     coding_exon_list = {}
#     for i, name in enumerate(coding_exon_names):
#         coding_exon_list[name.split('(')[0]] = coding_exon_seqs[i]
#
#     min_distances = {}
#     exon_list = []
#     for exon in ptc_positions:
#         seq = coding_exon_list[exon]
#         rel_pos = ptc_positions[exon]
#         min_dist = min(rel_pos, len(seq) - rel_pos)
#         min_distances[exon] = min_dist
#
#         if exon not in overlap:
#             exon_list.append(exon)
#
#     big_effects = [min_distances[exon] for exon in overlap]
#     other_effects = [min_distances[exon] for exon in exon_list]
#     print(ranksums(big_effects, other_effects))
#



def get_disease_coding_exons_makeup(results_prefix, disease_output_dir):
    '''
    Get the makeup of coding exons
    '''
    disease_ptc_file = "{0}_ptcs.txt".format(disease_output_dir)
    overlap_file = "{0}__analysis_overlaps.txt".format(disease_output_dir)
    coding_exons_bed = "{0}_coding_exons.bed".format(results_prefix)
    coding_exons = gen.read_many_fields(coding_exons_bed, "\t")

    coding_exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for exon in coding_exons:
        t = exon[3].split('.')[0]
        e = int(exon[3].split('.')[1])
        length = int(exon[2]) - int(exon[1])
        coding_exon_list[t][e] = length

    overlap_locations = []
    overlaps = gen.read_many_fields(overlap_file, "\t")
    for overlap in overlaps:
        overlap_locations.append(int(overlap[7]))

    disease_ptcs = gen.read_many_fields(disease_ptc_file, "\t")

    nts = ["A", "C", "G", "T"]
    nt_counts = [0,0,0]
    nt_identities = {}
    for nt in nts:
        nt_identities[nt] = 0

    for ptc in disease_ptcs:
        location = int(ptc[7])
        length = int(ptc[2]) - int(ptc[1])

        if location not in overlap_locations:
            strand = ptc[5]
            ref = ptc[9]
            if strand == "-":
                ref = gen.reverse_complement(ref)

            nt_identities[ref] += 1


            if length >= 6:
                nt_counts[0] += 6
                if length >= 138:
                    nt_counts[1] += 132
                    nt_counts[2] += (length-138)

                else:
                    nt_counts[1] += (length-6)
            else:
                nt_counts[0] += length

    regions = ["0-2","3-69","70+"]
    output_file = "{0}_ptc_info.txt".format(disease_output_dir)
    with open(output_file, "w") as outfile:
        outfile.write("PTC containing exon info\n\n")
        outfile.write("Region,Nucleotides,Proportion\n")
        for i, count in enumerate(nt_counts):
            outfile.write("{0},{1},{2}\n".format(regions[i], count, np.divide(count, sum(nt_counts))))
        outfile.write("Total nts,{0}\n".format(sum(nt_counts)))

        outfile.write("\Reference nucleotide counts\n")
        outfile.write("Nucleotide,Count,Proportion\n")
        for nt in sorted(nt_identities):
            outfile.write("{0},{1},{2}\n".format(nt, nt_identities[nt], np.divide(nt_identities[nt], sum(nt_identities.values()))))

def get_filtered_skipped_exons(final_output_file):

    lines = gen.read_many_fields(final_output_file, "\t")

    max_sample_count = max([float(line[2]) for line in lines[1:]])
    filtered_list = []

    # Exclude exons in which none of the individuals that we could quantify splicing contained a PTC,
    # and exons for which we could quantify splicing in less than half of samples
    for line in lines[1:]:
        sample_count = float(line[2])
        if float(line[1]) > 0 and sample_count > np.divide(max_sample_count, 2):
            filtered_list.append(line)

    return filtered_list

def chisquare_skipped_non_PTC_lengths(results_prefix, disease_output_dir):
    '''
    Chisqaure test of exon lengths for exons that are skipped despite no PTC.
    '''

    print("\nChisqaure test of exon lengths for exons that are skipped despite no PTC.\n")

    final_output_file = "{0}__analysis_final_output.txt".format(results_prefix)
    coding_exons_bed = "{0}_coding_exons.bed".format(results_prefix)

    filtered_list = get_filtered_skipped_exons(file_output_file)

    skipped_non_ptc_exons = []
    all_skipped_exons = []
    non_skipped_non_ptc_exons = []

    psi = []
    psi_het = []

    for line in filtered_list:
        exon = line[0]
        psi_non_ptc = float(line[5])
        if psi_non_ptc < 1:
            skipped_non_ptc_exons.append(exon)
            psi.append(float(line[5]))
            psi_het.append(float(line[4]))
        else:
            non_skipped_non_ptc_exons.append(exon)
        all_skipped_exons.append(exon)

    coding_exons = gen.read_many_fields(coding_exons_bed, "\t")
    coding_exons_list = {}
    for exon in coding_exons:
        length = int(exon[2]) - int(exon[1])
        coding_exons_list[exon[3]] = length

    exon_frames = []
    skipped_exon_frames = []
    not_skipped_exon_frames = []

    for exon in non_skipped_non_ptc_exons:
        length = coding_exons_list[exon]
        not_skipped_exon_frames.append(length%3)
    for exon in skipped_non_ptc_exons:
        length = coding_exons_list[exon]
        skipped_exon_frames.append(length%3)

    obs = [skipped_exon_frames.count(0), skipped_exon_frames.count(1), skipped_exon_frames.count(2)]
    exp = [np.divide(not_skipped_exon_frames.count(0), len(not_skipped_exon_frames)) * len(skipped_exon_frames), np.divide(not_skipped_exon_frames.count(1), len(not_skipped_exon_frames)) * len(skipped_exon_frames), np.divide(not_skipped_exon_frames.count(2), len(not_skipped_exon_frames)) * len(skipped_exon_frames)]

    print("counts, frame 0, frame 1, frame 2")
    print("obs, {0}".format(", ".join(gen.stringify(obs))))
    print("exp, {0}".format(", ".join(gen.stringify(exp))))

    print("Chisquare test: ", chisquare(obs, exp))



def main():

    description = "Miscellaneous tests."
    args = gen.parse_arguments(description,  ["results_prefix", "disease_output_dir", "ese_file", "skipped_lengths_chisquare", "disease_coding_exons_makeup", "ese_hits", "distance"], flags = [3,4,5,6])
    results_prefix, disease_output_dir, ese_file, skipped_lengths_chisquare, disease_coding_exons_makeup, ese_hits, distance = args.results_prefix, args.disease_output_dir, args.ese_file, args.skipped_lengths_chisquare, args.disease_coding_exons_makeup, args.ese_hits, args.distance

    if skipped_lengths_chisquare:
        chisquare_skipped_non_PTC_lengths(results_prefix, disease_output_dir)
    if disease_coding_exons_makeup:
        get_disease_coding_exons_makeup(results_prefix, disease_output_dir)
    if ese_hits:
        check_ese_hits(results_prefix, ese_file)
    if distance:
        check_distances(results_prefix)

if __name__ == "__main__":
    main()
