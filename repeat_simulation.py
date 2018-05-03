import generic as gen
import os
import collections
import numpy as np

def prepare_results():

    psi_values = collections.defaultdict(lambda: collections.defaultdict())

    results_dir = "./results/neg_control_outputs/"
    for file in os.listdir(results_dir):
        if "simulation" in file:
            simulation_number = int(file.split('.')[0].split('_')[-1])
            exons = gen.read_many_fields(results_dir + file, "\t")[1:]
            for exon in exons:
                exon_id = exon[0]
                psi_het_ptc = float(exon[4])
                psi_values[exon_id][simulation_number] = psi_het_ptc

    return psi_values

def main():
    psi_values = prepare_results()


if __name__ == "__main__":
    main()
