import simulate_eses_ops as se
import generic as gen

ese_sets = {
    "ESR": "./ese_data/ESR/ESR.txt",
    "Ke400_ESEs": "./ese_data/Ke400_ESEs/Ke400_ESEs.txt",
    "PESE": "./ese_data/PESE/PESE.txt",
    "RESCUE": "./ese_data/RESCUE/RESCUE.txt",
    "INT3": "./ese_data/CaceresHurstESEs_INT3/CaceresHurstESEs_INT3.txt",
    "RBP_motifs": "./ese_data/RBP_motifs/RBP_motifs.txt",
    "RBP_motifs_filtered": "./ese_data/RBP_motifs/filtered_RBP_motifs_nd.txt",
    "RBP_motifs_filtered_nonCDS": "./ese_data/RBP_motifs/filtered_RBP_motifs_nonCDS.txt",
}

def run_simulations(simulation_sets, required_simulations):

    for ese_set in simulation_sets:
        motifs = []
        motif_list = []

        gen.remove_file(ese_set[1])
        gen.remove_file(ese_set[2])

        #get motifs
        if ese_set[0] == "RBP_motifs_filtered":
            motifs = []
            motif_list = []
            #read the motif ids and seqs from fasta
            raw_motif_ids, raw_motifs = gen.read_fasta(ese_sets[ese_set[0]])
            index_list = []
            if len(ese_set) > 3 and ese_set[3] > 0:
                # restrict to rbp motifs with nd > 0
                index_list = [i for i, motif_id in enumerate(raw_motif_ids) if float(motif_id.split('|')[1]) > 0]
            elif len(ese_set) > 3 and ese_set[3] < 0:
                # restrict to rbp motifs with nd > 0
                index_list = [i for i, motif_id in enumerate(raw_motif_ids) if float(motif_id.split('|')[1]) < 0]
            else:
                index_list = [i for i, motif_id in enumerate(raw_motif_ids)]
            #extract the motifs
            [motif_list.extend(raw_motifs[i].split('|') for i in index_list)]
            [motifs.extend(i) for i in motif_list]

        elif ese_set[0] == "RBP_motifs" or ese_set[0] == "RBP_motifs_filtered_nonCDS":
            raw_motif_ids, raw_motifs = gen.read_fasta(ese_sets[ese_set[0]])
            motif_list.extend(i.split('|') for i in raw_motifs)
            [motifs.extend(i) for i in motif_list]

        else:

            motif_list = gen.read_many_fields(ese_sets[ese_set[0]], ",")
            # get motifs, avoid header
            motifs = [i[0] for i in motif_list if i[0][0] != "#"]


        #if the rbps are split
        if len(ese_set) > 3:

            print('Starting {0} {1} set...'.format(ese_set[0], ese_set[4]))
        else:
            print('Starting {0} set...'.format(ese_set[0]))

        # get the number of stop codons found in the real set
        real_count = se.get_stop_codon_count(motifs)
        # generate simulated motifs using motif set
        print('Simulating...'.format(ese_set[0]))
        se.generate_motifs_sets(motifs, required_simulations, output_file=ese_set[1])
        simulated_motif_sets = gen.read_many_fields(ese_set[1], "|")
        with open(ese_set[2], "w") as output:
        	output.write('id,stop_count\n')
        	output.write('real,{0}\n'.format(real_count))
        	for i, simulated_set in enumerate(simulated_motif_sets):
        		stop_count = se.get_stop_codon_count(simulated_set)
        		output.write('{0},{1}\n'.format(i+1, stop_count))

def main():

    description = "Check whether stop codons are depleted in motif sets by simulating the motif set."
    args = gen.parse_arguments(description, ["required_simulations", "all_sets", "ESR", "Ke", "PESE", "RESCUE", "INT3", "RBP_motifs", "filter_RBPs", "split_RBPs", "RBP_motifs_filtered_nonCDS"], flags = [1,2,3,4,5,6,7,8,9,10])
    required_simulations, all_sets, ESR, Ke, PESE, RESCUE, INT3, RBP_motifs, filter_RBPs, split_rbps, RBP_motifs_filtered_nonCDS = args.required_simulations, args.all_sets, args.ESR, args.Ke, args.PESE, args.RESCUE, args.INT3, args.RBP_motifs, args.filter_RBPs, args.split_RBPs, args.RBP_motifs_filtered_nonCDS


    if split_rbps and not filter_RBPs:
        print('You must specify the filtered RBPs if you want to split by ND.')
        raise Exception

    if not required_simulations:
        print('You must specify the number of simulations you require.')
        raise Exception

    #create the output_directory
    output_directory = "ese_stops_simulations_output"
    gen.create_directory(output_directory)

    #set up the simulations we want
    required_sets = []
    if all_sets:
        required_sets.extend([i for i in ese_sets])
    else:
        if ESR:
            required_sets.append("ESR")
        if Ke:
            required_sets.append("Ke400_ESEs")
        if PESE:
            required_sets.append("PESE")
        if RESCUE:
            required_sets.append("RESCUE")
        if INT3:
            required_sets.append("INT3")
        if RBP_motifs and not filter_RBPs:
            required_sets.append("RBP_motifs")
        if RBP_motifs and filter_RBPs:
            required_sets.append("RBP_motifs_filtered")
        if RBP_motifs_filtered_nonCDS:
            required_sets.append("RBP_motifs_filtered_nonCDS")

    #check whether any sets have been chosen
    if len(required_sets) == 0:
        print("\nPlease choose a motif set to analyse:\n")
        [print("--{0}".format(i)) for i in sorted(ese_sets)]
        print("\n")
        raise Exception

    #create the necessary files
    simulation_sets = []
    for ese_set in required_sets:
        if ese_set == "RBP_motifs_filtered":
            dir_name = "RBP_motifs"
        else:
            dir_name = ese_set
        #create the output directory for the particular motif set
        motif_output_directory = "{0}/{1}".format(output_directory, dir_name)
        gen.create_directory(motif_output_directory)
        if split_rbps:
            #if we want to split the rbp motifs based on nd, need to create 2 lots of outputs
            simulated_set_output_pos_nd = "{0}/{1}_simulants_pos_nd_{2}.txt".format(motif_output_directory, dir_name, required_simulations)
            output_file_pos_nd = "{0}/{1}_stop_counts_pos_nd_{2}.csv".format(motif_output_directory, dir_name, required_simulations)
            simulation_sets.append([ese_set, simulated_set_output_pos_nd, output_file_pos_nd, 1, "Positive ND"])
            simulated_set_output_neg_nd = "{0}/{1}_simulants_neg_nd_{2}.txt".format(motif_output_directory, dir_name, required_simulations)
            output_file_neg_nd = "{0}/{1}_stop_counts_neg_nd_{2}.csv".format(motif_output_directory, dir_name, required_simulations)
            simulation_sets.append([ese_set, simulated_set_output_neg_nd, output_file_neg_nd, -1, "Negative ND"])
        else:
            #create simulated set output, analysis output file
            simulated_set_output = "{0}/{1}_simulants_{2}.txt".format(motif_output_directory, dir_name, required_simulations)
            output_file = "{0}/{1}_stop_counts_{2}.csv".format(motif_output_directory, dir_name, required_simulations)
            simulation_sets.append([ese_set, simulated_set_output, output_file])

    run_simulations(simulation_sets, int(required_simulations))

if __name__ == "__main__":
    main()
