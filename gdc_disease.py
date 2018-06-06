'''
Author: Rosina Savisaar and Liam Abrahams
'''

import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import generic as gen
import gdc_ops as gdco
import os
import random
import SNP_ops as so
import time
import numpy as np
import shutil

def main():

    description = "Analysis of disease mutations from GDC."
    args = gen.parse_arguments(description,  ["output_dir", "mutations_dir"], flags = [])
    output_dir, mutations_dir = args.output_dir, args.mutations_dir

    start = time.time()
    # create the output directories
    gen.create_output_directories('temp_data/')
    gen.create_output_directories(output_dir)

    # prcoess the mutation files
    processed_mutations = "{0}/processed_mutations.vcf".format(output_dir)
    gdco.process_mutation_files(mutations_dir, output_dir, processed_mutations)


if __name__ == "__main__":
    main()
