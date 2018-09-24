'''
Author: Liam Abrahams
Operations used to simulate eses.
'''

import generic as gen
import itertools as it
import re
import numpy as np
import os

def make_simulants(motifs, n_sim, output_file_name = None, retro = False, mono = False, no_duplicates = False, remove_stops = False, remove_existing = False, cap_runs = False, exclude = None, seed = None, concat = True):
    '''
    Given a set of motifs, generate n_sim sets of simulants with the same over-all dinucleotide frequencies.
    motifs: a list of strings (the motifs for which you want to generate simulants)
    n_sim: number of simulant sets required
    mono: use mono- rather than dinucleotide frequencies
    no_duplicates: don't allow duplicates within simulant sets
    remove_stops: don't allow simulants that contain stop codons
    remove_existing: don't allow simulants to coincide with motifs in the true set
    cap_runs: don't allow mononucleotide runs whose length exceeds that of the longst run of that base in the true motifs
    exclude: don't allow any of the motifs within this list
    seed: supply a seed for the PRNG
    concat: concatenate the true motifs before extracting dinucleotides
    '''

    if cap_runs:
        longest_runs_numbers = get_longest_run(motifs)
        longest_runs_strings = ["".join([base for i in range(longest_runs_numbers[base] + 1)]) for base in longest_runs_numbers]
    motifs = [list(i) for i in motifs]
    nts = flatten(motifs)
    if mono:
        dints = flatten(motifs)
    else:
        if concat:
            dints = []
            #grab all the dinucleotides in the two reading frames
            for i in range(0, len(nts) - 1, 2):
                dints.append([nts[i], nts[i+1]])
            for i in range(1, len(nts) - 1, 2):
                dints.append([nts[i], nts[i+1]])
        else:
            dints = []
            for motif in motifs:
                for i in range(0, len(motif) - 1, 2):
                    dints.append([motif[i], motif[i+1]])
                for i in range(1, len(motif) - 1, 2):
                    dints.append([motif[i], motif[i+1]])
    # print(dints)
    # #right now you have a list of lists. Turn that into a list of strings where each string is one dinucleotide.
    dints = ["".join(i) for i in dints]
    print(dints)




def get_dinucleotides(motifs, concat_motifs=None):

	'''
	Generate a list of dinucleotides for a given set of motifs.
	concat_motifs: join all motifs in string and then get dinucelotides.
	'''

	dinucleotides = []

	if concat_motifs:
		#concat the motifs to a string
		motif_string = "".join([i for i in motifs])
		#compile the regex to find all instances of two characters
		dint_regex = re.compile(".{2}")
		#search for all dnts in the two reading frames
		dinucleotides.extend(re.findall(dint_regex, motif_string))
		dinucleotides.extend(re.findall(dint_regex, motif_string[1:]))
	else:
		motifs = [list(i) for i in motifs]
		for motif in motifs:
			#search for all dnts in the two reading frames
			for i in range(0, len(motif)-1, 2):
				dinucleotides.append([motif[i], motif[i+1]])
			for i in range(1, len(motif)-1, 2):
				dinucleotides.append([motif[i], motif[i+1]])
		#join the nts comprising a dint back together
		dinucleotides = ["".join(i) for i in dinucleotides]

	return(dinucleotides)

def generate_motifs(simulations, motifs, dinucleotides, seed=None):

	'''
	Generate a random set of motifs given a set of motifs.
	Generate using a set of dinucleotides from the set of motifs.
	'''

	#create any empty set of motifs
	simulant_set = []
	#get a list of all the nucleotides in the set
	motif_nts = list(it.chain(*motifs))

	for i, simulation in enumerate(simulations):
		#set the randomisation seed
		if seed:
			#chunk seeds based on processes
			seed_chunks = [seed[i] for i in simulations]
			np.random.seed(seed_chunks[i])
		else:
			np.random.seed()

		created_simulants = []
		for i, motif in enumerate(motifs):
			#detemine whether the motif is of odd length
			if len(motif) % 2 == 0:
				odd = False
			else:
				odd = True

			#get the number of required dinucleotides to create simulant
			required_dinucleotides = len(motif) // 2
			created = False
			new_simulant = []
			while not created:
				#pick the number of required dinucleotides
				new_simulant.extend(np.random.choice(dinucleotides, required_dinucleotides))
				#if odd length, add one random nt from the nt list
				if odd:
					new_simulant.extend(np.random.choice(motif_nts, 1))
				#if the newly created simulant isnt in the set of simulated motifs, keep
				if new_simulant not in created_simulants:
					created = True
					created_simulants.append(new_simulant)
				else:
					#otherwise reset new simulant
					new_simulant = []
		simulants = ["".join(i) for i in created_simulants]
		simulant_set.append(simulants)
	return(simulant_set)

def generate_motifs_sets(motifs, simulations_to_run, output_file, seed_list=None, onebyone=None):

	'''
	Generate n sets of motifs based on the set of motifs provided.
	seed_list: a list of seeds to use (must be of length greater or equal to the number of simulations)
	'''

	#check that there are enough seeds if the seed is set
	if seed_list and simulations_to_run > len(seed_list):
		print('The number of seeds must be at least equal to the number of simulations!')
		raise Exception

	#get dinucleotides
	dinucleotides = get_dinucleotides(motifs)
	#create a list of processes
	input_list = [i for i in range(simulations_to_run)]
	#build processes
	processes = gen.run_in_parallel(input_list, ["foo", motifs, dinucleotides, seed_list], generate_motifs, onebyone)
	#run processes and output to output_file
	output = open(output_file, "w")
	for process in processes:
		simulants = process.get()
		if simulants:
			for simulant in simulants:
				output.write('{0}\n'.format("|".join(simulant)))
	output.close()

def get_stop_codon_count(motifs):
	'''
	Get the number of stop codons in any frame in a list of motifs.
	'''
	stop_regex = re.compile('TAA|TAG|TGA')
	stops = [len(re.findall(stop_regex, i)) for i in motifs]
	return(sum(stops))
