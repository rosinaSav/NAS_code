'''
Author: Rosina Savisaar.
Take an output file from prepare_FANTOM.py and make a file with the expression data for each gene.
'''

import bam_ops as bmo
import bed_ops as bo
import generic as gen
import numpy as np
import re
import scipy.stats

def main():

    description = "Take an output file from prepare_FANTOM.py and make a file with the expression data for each gene."
    args = gen.parse_arguments(description, ["clean_fasta", "promoters_file_name", "cage_file_name", "out_prefix", "TPM_threshold"], ints = [4])
    [clean_fasta, promoters_file_name, cage_file_name, out_prefix, TPM_threshold] = [args.clean_fasta, args.promoters_file_name, args.cage_file_name, args.out_prefix, args.TPM_threshold]

    #extract transcript coordinates
    transcripts_file = "{0}_transcripts_clean.bed".format(out_prefix)
    bo.extract_features("../source_data/Homo_sapiens.GRCh37.87.gtf", transcripts_file, ["transcript"])

    #get the names of the transcripts you're interested
    names = gen.read_fasta(clean_fasta)[0]

    #write the coordinates of the promoter regions of those transcripts to file
    with open(promoters_file_name, "w") as out_file, open(transcripts_file, "r") as in_file:
        for line in in_file:
            parsed = (line.rstrip("\n")).split("\t")
            #parse out the transcript name
            name = parsed[3].split(".")[0]
            #skip transcripts that aren't among your transcripts of interest
            if name in names:
                #determine the coordinates of a 1001 bp region centered on the TSS (the supposed promoter region)
                if parsed[5] == "+":
                    current_line = ["chr" + parsed[0], int(parsed[1]) - 500, int(parsed[1]) + 500 + 1, name, ".", parsed[5]]
                elif parsed[5] == "-":
                    current_line = ["chr" + parsed[0], int(parsed[2]) - 500 - 1, int(parsed[2]) + 500, name, ".", parsed[5]]
                else:
                    RuntimeError("Invalid strand information!")
                out_file.write("\t".join([str(i) for i in current_line]))
                out_file.write("\n")

    #check which CAGE peaks overlap which promoters
    overlapping_peaks_file = "{0}_FANTOM_overlap_peaks.bed".format(out_prefix)
    bmo.intersect_bed(cage_file_name, promoters_file_name, output_file = overlapping_peaks_file, force_strand = True, write_both = True, no_dups = False)

    #for each transcript, get all overlapping peaks
    #(store only the expression information)
    peaks_dict = {name: [] for name in names}
    with open(overlapping_peaks_file, "r") as peaks:
        for peak in peaks:
            peak = peak.split("\t")
            name = peak[9]
            peaks_dict[name].append(peak[3])
            
    #for each transcript,
    #store the mean TPM within each tissue (averaged over the different peaks
    #associated to that transcript)
    mean_dict = {}
    np.set_printoptions(suppress = True)
    for name in peaks_dict:
        if len(peaks_dict[name]) > 0:
            current_mat = np.array([[float(j) for j in i.split("|")] for i in peaks_dict[name]])
            means = np.mean(current_mat, axis = 0)
            mean_dict[name] = means

    #calculate expression parameters
    final_dict = {}
    for gene in mean_dict:
        expressed = len([i for i in mean_dict[gene] if i > TPM_threshold])
        fraction = expressed/len(mean_dict[gene])
        maximum = np.max(mean_dict[gene])
        median_expr = np.median(mean_dict[gene])
        median_if_expressed = np.median([i for i in mean_dict[gene] if i > TPM_threshold])
        final_dict[gene] = [fraction, maximum, median_expr, median_if_expressed]


    output_file_name = "{0}_FANTOM_expression_per_transcript.txt".format(out_prefix)      
    with open(output_file_name, "w") as file:
        file.write("gene\tbreadth\tmax\tmedian\tmedian_expr\n")
        for i in sorted(list(final_dict.keys())):
            if final_dict[i] != None:
                file.write("\t".join([i] + [str(j) for j in final_dict[i]]))
                file.write("\n")               
     

if __name__ == "__main__":
    main()
