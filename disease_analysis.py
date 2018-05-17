'''
Author: Rosina Savisaar and Liam Abrahams
Check whether PTCs are associated with greater rates of exon skipping.
'''

import bed_ops as bo
import bam_ops as bmo
import disease_ops as do
import generic as gen
import os
import random
import SNP_ops as so
import time
import numpy as np

# To do:
# intersect the snp files with the coding exons
# check the mutation status of the snps
# then check the overlaps

def main():

    description = "Analysis of disease mutations from firehose broad."
    args = gen.parse_arguments(description,  ["output_dir", "mutations_dir", "results_prefix", "process_mutations", "intersect_mutations"], flags = [3,4], ints = [])
    output_dir, mutations_dir, results_prefix, process_mutations, intersect_mutations = args.output_dir, args.mutations_dir, args.results_prefix, args.process_mutations, args.intersect_mutations

    start = time.time()
    gen.create_directory('temp_data/')

    # create folder containing processed mutation files
    processed_dir = "{0}/processed_mutation_files".format(output_dir)
    filename_prefix = "processed_mutations"
    full_mutation_file = "{0}/mutations.txt".format(output_dir)
    gen.create_output_directories(processed_dir)
    if process_mutations:
        print("Processing mutation files...")
        entrylimit = 10000
        do.refactor_files(mutations_dir, processed_dir, filename_prefix, full_mutation_file, entrylimit, clean_directory = True)

    coding_exons = "{0}_coding_exons.bed".format(results_prefix)
    # if intersect_mutations:
    #     print("Intersecting mutations with coding exons...")


if __name__ == "__main__":
    main()
