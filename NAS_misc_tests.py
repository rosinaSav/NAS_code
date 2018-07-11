import generic as gen
import os
import numpy as np
from scipy.stats import chisquare
import collections



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

def chisquare_skipped_non_PTC_lengths(results_prefix, disease_output_dir):
    '''
    Chisqaure test of exon lengths for exons that are skipped despite no PTC.
    '''

    print("\nChisqaure test of exon lengths for exons that are skipped despite no PTC.\n")

    final_output_file = "{0}__analysis_final_output.txt".format(results_prefix)
    coding_exons_bed = "{0}_coding_exons.bed".format(results_prefix)

    lines = gen.read_many_fields(final_output_file, "\t")

    skipped_non_ptc_exons = []
    all_skipped_exons = []

    for line in lines[1:]:
        exon = line[0]
        psi_non_ptc = float(line[5])
        if psi_non_ptc < 1:
            skipped_non_ptc_exons.append(exon)
        all_skipped_exons.append(exon)

    coding_exons = gen.read_many_fields(coding_exons_bed, "\t")

    exon_frames = []
    skipped_exon_frames = []

    for exon in coding_exons:
        start = int(exon[1])
        end = int(exon[2])
        length = end - start

        if exon[3] in skipped_non_ptc_exons:
            skipped_exon_frames.append(length % 3)
        elif exon[3] not in all_skipped_exons:
            exon_frames.append(length % 3)


    obs = [skipped_exon_frames.count(0), skipped_exon_frames.count(1), skipped_exon_frames.count(2)]
    exp = [np.divide(exon_frames.count(0), len(exon_frames)) * len(skipped_exon_frames), np.divide(exon_frames.count(1), len(exon_frames)) * len(skipped_exon_frames), np.divide(exon_frames.count(2), len(exon_frames)) * len(skipped_exon_frames)]

    print("counts, frame 0, frame 1, frame 2")
    print("obs, {0}".format(", ".join(gen.stringify(obs))))
    print("exp, {0}".format(", ".join(gen.stringify(exp))))

    print("Chisquare test: ", chisquare(obs, exp))



def main():

    description = "Miscellaneous tests."
    args = gen.parse_arguments(description,  ["results_prefix", "disease_output_dir", "skipped_lengths_chisquare", "disease_coding_exons_makeup"], flags = [2,3])
    results_prefix, disease_output_dir, skipped_lengths_chisquare, disease_coding_exons_makeup = args.results_prefix, args.disease_output_dir, args.skipped_lengths_chisquare, args.disease_coding_exons_makeup

    if skipped_lengths_chisquare:
        chisquare_skipped_non_PTC_lengths(results_prefix, disease_output_dir)
    if disease_coding_exons_makeup:
        get_disease_coding_exons_makeup(results_prefix, disease_output_dir)

if __name__ == "__main__":
    main()
