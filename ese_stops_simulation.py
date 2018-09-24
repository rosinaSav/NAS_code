'''
Author: Liam Abrahams.
Check whether stop codons are depleted in motif sets by simulating the motif set.
'''

import simulate_eses_ops as se
import generic as gen
import numpy as np

def run_simulations(simulation_sets, required_simulations):

    for motif_set in simulation_sets:

        motif_file = motif_set[0]
        simulation_output_file = motif_set[1]
        stops_count_output_file = motif_set[2]

        # clean up and previous simulations
        gen.remove_file(simulation_output_file)
        gen.remove_file(stops_count_output_file)

        motif_list = gen.read_many_fields(motif_file, ",")
        # get motifs, avoid header if there is one
        motifs = [i[0] for i in motif_list if i[0][0] != "#"]

        # get the number of stop codons found in the real set
        real_count = se.get_stop_codon_count(motifs)

        # generate simulated motifs using motif set
        print('Simulating {0}...'.format(motif_file))
        se.generate_motifs_sets(motifs, required_simulations, output_file = simulation_output_file)
        simulated_motif_sets = gen.read_many_fields(simulation_output_file, "|")
        with open(stops_count_output_file, "w") as output:
        	output.write('id,stop_count\n')
        	output.write('real,{0}\n'.format(real_count))
        	for i, simulated_set in enumerate(simulated_motif_sets):
        		stop_count = se.get_stop_codon_count(simulated_set)
        		output.write('{0},{1}\n'.format(i+1, stop_count))

def get_possible_hexamers(exon):
    hexamers = []
    for i in range(len(exon)):
        if i+6 <= len(exon):
            hexamers.append(exon[i:i+6])
    return hexamers


def generate_motifs(exon_seqs, required_motifs):
    '''
    Get random hexamers from the exon seqs
    '''
    simulant_choices = []
    exon_choices = list(range(len(exon_seqs)))
    while len(simulant_choices) < required_motifs:
        exon_choice = exon_seqs[np.random.choice(exon_choices)]
        hexamers = get_possible_hexamers(exon_choice)
        random_hexamer = np.random.choice(hexamers)
        if random_hexamer not in simulant_choices:
            simulant_choices.append(random_hexamer)

    return simulant_choices


def simulate_motifs(simulations, exon_seqs, motifs):

    outputs = []

    for simulation in simulations:
        print("Simulation {0}".format(simulation+1))
        seed = np.random.randint(0, 10000)
        np.random.seed(seed)
        motifs = generate_motifs(exon_seqs, len(motifs))
        stop_count = se.get_stop_codon_count(motifs)
        outputs.append(stop_count)

    return outputs

def run_exon_simulation(motif_file, exon_fasta, output_dir, required_simulations, output_file):
    '''
    Run simulation that picks hexamers from the exon sequences
    '''

    exon_names, exon_seqs = gen.read_fasta(exon_fasta)
    #exons needs to be >= 16 to get the two exon ends
    exon_seqs = [exon for exon in exon_seqs if len(exon) >= 16]

    # get motifs, avoid header if there is one
    motif_list = gen.read_many_fields(motif_file, ",")
    motifs = [i[0] for i in motif_list if i[0][0] != "#"]

    real_count = se.get_stop_codon_count(motifs)

    simulations = list(range(required_simulations))
    # simulated_counts = simulate_motifs(simulations, exon_seqs, motifs)
    processes = gen.run_in_parallel(simulations, ["foo", exon_seqs, motifs], simulate_motifs)

    outputs = []
    for process in processes:
        outputs.extend(process.get())

    with open(output_file, "w") as outfile:
        outfile.write('sim,count\n')
        outfile.write('real,{0}\n'.format(real_count))
        for i, count in enumerate(outputs):
            outfile.write('{0},{1}\n'.format(i+1,count))




def main():

    description = "Check whether stop codons are depleted in motif sets by simulating the motif set."
    args = gen.parse_arguments(description, ["motif_file", "output_dir", "results_dir", "required_simulations", "motif_simulation", "exon_simulation"], flags = [4,5], ints = [3])
    motif_file, output_dir, results_dir, required_simulations, motif_simulation, exon_simulation = args.motif_file,  args.output_dir, args.results_dir, args.required_simulations, args.motif_simulation, args.exon_simulation

    if not required_simulations:
        print('You must specify the number of simulations you require.')
        raise Exception

    gen.create_output_directories(output_dir)

    if motif_simulation:
        simulation_sets = []

        #create the output directory for the particular motif set
        motif_output_dir = "{0}/{1}".format(output_dir, ".".join(motif_file.split('.')[:-1]).split('/')[-1])
        gen.create_output_directories(motif_output_dir)


        simulated_motifs_output = "{0}/simulations_{1}.txt".format(motif_output_dir, required_simulations)
        output_file = "{0}/stop_counts_{1}.txt".format(motif_output_dir, required_simulations)

        # add the files to the required list
        simulation_sets.append([motif_file, simulated_motifs_output, output_file])

        # run the simulations
        run_simulations(simulation_sets, required_simulations)

    exon_hexamer_simulation = "{0}/region_hexamer_sim.csv".format(output_dir)
    if exon_simulation:
        exon_fasta = "{0}_CDS_intervals.fasta".format(results_dir)
        run_exon_simulation(motif_file, exon_fasta, output_dir, required_simulations, exon_hexamer_simulation)


if __name__ == "__main__":
    main()
