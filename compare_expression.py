'''
Author: Rosina Savisaar.
Compare expression parameters for transcripts that contain true PTCs vs transcripts that contain pseudo-PTCs.
'''

import generic as gen
import numpy as np
import os

def display_comparison(parameter, true_data, sim_data, column):
    print("True {0}: {1}".format(parameter, true_data[column]))
    print("Median simulated {0}: {1}".format(parameter, np.median(sim_data[:, column])))
    print("One-tailed empirical p: {0}".format(gen.calc_eff_p(true_data[column], sim_data[:, column], greater = False)))
    print("\n")

def PTCs_to_expression(file, expression):
    '''
    Given a list of lists of expression data and a PTCs file, calculate median value
    for each of the expression parameters across the relevant transcripts.
    '''
    #make a list of all the transcripts that appear in the PTCs file
    #if there are several exons from the same transcript, then we will count that transcript once
    names = []
    with open(file, "r") as ptcs:
        for ptc in ptcs:
            ptc = ptc.split("\t")
            names.append(ptc[3].split(".")[0])

    #only retain relevant rows from expression data (transcripts that appear in the PTCs file)
    expression = np.array([i[1:] for i in expression if i[0] in names], dtype = "float")
    col_means = np.median(expression, axis = 0)

    return(col_means)

def main():
    description = "Compare expression parameters for transcripts that contain true PTCs vs transcripts that contain pseudo-PTCs."
    args = gen.parse_arguments(description, ["PTCs_file", "pseudo_PTCs_directory", "expression_file"])
    PTCs_file, pseudo_PTCs_directory, expression_file = args.PTCs_file, args.pseudo_PTCs_directory, args.expression_file

    expression = gen.read_many_fields(expression_file, "\t")

    #get median expression parameters for true PTCs
    true_values = PTCs_to_expression(PTCs_file, expression)

    #do the same for each of the simulant PTC files
    sim_files = os.listdir(pseudo_PTCs_directory)
    #I'm doing the first one separately so I could easily stack the outputs
    sim_values = PTCs_to_expression("{0}/{1}".format(pseudo_PTCs_directory, sim_files[0]), expression)
    for sim_file in sim_files[1:]:
        curr_sim_values = PTCs_to_expression("{0}/{1}".format(pseudo_PTCs_directory, sim_file), expression)
        sim_values = np.vstack((sim_values, curr_sim_values))

    display_comparison("breadth", true_values, sim_values, 0)
    display_comparison("maximum TPM", true_values, sim_values, 1)
    display_comparison("median TPM", true_values, sim_values, 2)
    display_comparison("median TPM (if expressed)", true_values, sim_values, 3)
 

if __name__ == "__main__":
    main()
