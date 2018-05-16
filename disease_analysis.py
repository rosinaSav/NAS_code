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
    args = gen.parse_arguments(description,  ["output_dir", "mutations_dir", "process_mutations"], flags = [2], ints = [])
    output_dir, mutations_dir, process_mutations = args.output_dir, args.mutations_dir, args.process_mutations

    start = time.time()
    gen.create_directory('temp_data/')

    # create folder containing processed mutation files
    processed_dir = "{0}/processed_mutation_files".format(output_dir)
    filename_prefix = "processed_mutations"
    gen.create_output_directories(processed_dir)
    if process_mutations:
        print("Processing mutation files...")
        entrylimit = 10000
        do.refactor_files(mutations_dir, processed_dir, filename_prefix, entrylimit, clean_directory = True)


if __name__ == "__main__":
    main()
