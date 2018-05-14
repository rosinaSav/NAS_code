'''
Auhor: Liam Abrahams.
Look at disease SNPs from ClinVar.
'''

import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os
import collections
import re
import copy
import numpy as np

def ptc_location_simulation(snp_file, full_bed, cds_fasta, output_directory, required_simulations, coding_exons_file):
    '''
    Simulate the snp location.
    For each snp, pick another site that has the same reference allele and that would generate a ptc with the mutated allele.
    Repeat n times.
    '''

    # return all the possible_locations
    possible_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    nts = ["A", "C", "G", "T"]
    for nt in nts:
        location_file = "{0}/possible_ptc_locations_{1}.fasta".format(output_directory, nt)
        entry_names, entry_locations = gen.read_fasta(location_file)
        for i, name in enumerate(entry_names):
            exon = name.split(':')[0]
            aa = name.split(':')[1][0]
            ma = name.split(':')[1][-1]
            possible_locations[exon][aa][ma].append(entry_locations[i])

    exons = gen.read_many_fields(coding_exons_file, "\t")
    exon_list = {}
    for exon in exons:
        exon_list[exon[3]] = int(exon[2]) - int(exon[1])

    # create a list of required simulations
    simulations = list(range(1, int(required_simulations) + 1))
    run_location_simulations(simulations, snp_file, possible_locations, exon_list, output_directory)

def run_location_simulations(simulations, snp_file, possible_locations, exon_list, output_directory):
    for simulation in simulations:
        print('{0}/{1}...'.format(simulation, len(simulations)))
        output_file = "{0}/disease_location_simulation_{1}.bed".format(output_directory, simulation)
        generate_pseudo_snps(snp_file, possible_locations, exon_list, output_file)


def generate_pseudo_snps(snp_file, possible_locations, exon_list, output_file):
    snps = gen.read_many_fields(snp_file, "\t")

    with open(output_file, "w") as outfile:
        np.random.seed()
        for snp in snps[1:]:
            snp_info = EntryInfo(snp)
            if len(possible_locations[snp_info.transcript_id][snp_info.aa][snp_info.ma]):
                possible_ptcs = possible_locations[snp_info.transcript_id][snp_info.aa][snp_info.ma][0].split(',')
                choices = list(range(len(possible_ptcs)))
                choice = np.random.choice(choices)
                pseudo_ptc_choice = possible_ptcs[choice]
                exon_index = int(re.findall('\d+', pseudo_ptc_choice.split(':')[0])[0])
                cds_index = re.findall('\d+', pseudo_ptc_choice.split(':')[1])[0]
                exon_number = re.findall('\d+', pseudo_ptc_choice.split(':')[2])[0]

                id = "{0}.{1}".format(snp_info.transcript_id, exon_number)
                exon_length = exon_list[id]
                min_dist = min(exon_index, int(exon_length) - int(exon_index))

                snp.extend([id, "{0}".format(exon_index), "{0}".format(cds_index), str(min_dist)])
                outfile.write("{0}\n".format("\t".join(snp)))


def generate_possible_ptc_locations(full_bed, cds_fasta, output_directory):

    stop_bases = ["A", "G", "T"]    # a mutation to c cant generate a stop so its not included
    stop_codons = ["TAA", "TAG", "TGA"]

    # get the list of coding sequnces
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    cds_list = {}
    for i, name in enumerate(cds_names):
        cds_list[name] = cds_seqs[i]

    # get all the indicies of mutable positions in the cds that could generate an in frame stop
    mutation_index_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for x, name in enumerate(cds_list):
        seq = cds_list[name]
        inframe_codons = re.findall('.{3}', seq[:-3])
        for i, codon in enumerate(inframe_codons):
            nts = list(codon)
            for j, nt in enumerate(nts):
                for mut in stop_bases:
                    if mut != nt:                       # we dont need 'mutations' that are the same as the ref base
                        mut_nts = copy.deepcopy(nts)
                        mut_nts[j] = mut
                        mut_codon = "".join(mut_nts)
                        if mut_codon in stop_codons:
                            mut_index = (i*3) + j
                            mutation_index_list[name][nt][mut].append(mut_index)

    # now also get the exon lengths
    exons = gen.read_many_fields(full_bed, "\t")
    exon_list_lengths = collections.defaultdict(lambda: collections.defaultdict())
    for exon in exons:
        start = int(exon[1])
        stop = int(exon[2])
        transcript_id = exon[3].split('.')[0]
        exon_id = int(exon[3].split('.')[1])
        type = exon[4]
        if type == "CDS":
            exon_list_lengths[transcript_id][exon_id] = stop-start

    # get all the possible cds indicies but grouped into exons
    exon_list_indices = collections.defaultdict(lambda: collections.defaultdict())
    for transcript in exon_list_lengths:
        current_end = 0
        current_start = 0
        for i, exon in enumerate(sorted(exon_list_lengths[transcript])):
            # get the length of the exon, add to the current length
            exon_length = exon_list_lengths[transcript][exon]
            current_end += exon_length
            # generate the indicies for that exon
            indices = list(range(current_start, current_end))
            exon_list_indices[transcript][exon] = [indices, current_start]
            # now move the start
            current_start += exon_length

    mutation_matched_indices = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    # for each of the possible ptc generating mutations
    for transcript in mutation_index_list:
        for ra in mutation_index_list[transcript]:
            for ma in mutation_index_list[transcript][ra]:
                # get the list of cds indices
                cds_index_list = mutation_index_list[transcript][ra][ma]
                for cds_index in cds_index_list:
                    for exon in exon_list_indices[transcript]:
                        if cds_index in exon_list_indices[transcript][exon][0]:
                            start_index = exon_list_indices[transcript][exon][1]
                            mutation_matched_indices[transcript][ra][ma].append("eidx{0}:cidx{1}:eid{2}".format(cds_index-start_index, cds_index, exon))

    nts = ["A", "C", "G", "T"]
    possible_mutation_files = {}
    for nt in nts:
        output_file = "{0}/possible_ptc_locations_{1}.fasta".format(output_directory, nt)
        with open(output_file, "w") as outfile:
            for transcript in mutation_matched_indices:
                for ma in sorted(mutation_matched_indices[transcript][nt]):
                    outfile.write(">{0}:{1}-{2}\n".format(transcript, nt, ma))
                    outfile.write("{0}\n".format(",".join(mutation_matched_indices[transcript][nt][ma])))


class EntryInfo(object):
    def __init__(self, entry):
        self.chr = entry[0]
        self.exon_start = int(entry[1])
        self.exon_stop = int(entry[2])
        self.transcript_id = entry[3].split('.')[0]
        self.exon_id = int(entry[3].split('.')[1])
        self.strand = entry[5]
        self.snp_pos = int(entry[7])
        self.snp_id = int(entry[8])
        self.aa = entry[9]
        self.ma = entry[10]
        self.rel_pos = int(entry[11])
        self.type = entry[12]
        self.given_id = entry[-1]

def get_exon_list(file):
    exons = gen.read_many_fields(file, "\t")
    exon_list = {}
    for exon in exons:
        exon_info = EntryInfo(exon)
        exon_list[exon_info.snp_id] = exon
    return(exon_list)

def get_ptc_info(ptc_file, relative_exon_positions_file, output_file):
    '''
    Get basic info about the location of snp.
    '''
    exon_list = get_exon_list(relative_exon_positions_file)

    ptcs = gen.read_many_fields(ptc_file, "\t")
    with open(output_file, "w") as outfile:
        header_list = ["ptc_id", "transcript_id", "exon_id", "aa", "ma", "exon_length", "5_dist", "3_dist", "min_dist"]
        outfile.write("{0}\n".format(",".join(header_list)))
        for ptc in ptcs:
            ptc_info = EntryInfo(ptc)
            ptc = EntryInfo(exon_list[ptc_info.snp_id])
            # check the type
            if ptc_info.type != ".":
                exon_length = ptc.exon_stop - ptc.exon_start
                three_prime_dist = exon_length - ptc.rel_pos
                output_list = gen.stringify([ptc_info.given_id, ptc.transcript_id, ptc.exon_id, ptc.aa, ptc.ma, exon_length, ptc.rel_pos, three_prime_dist, min(ptc.rel_pos, three_prime_dist)])
                outfile.write("{0}\n".format(",".join(output_list)))

def main():

    description = "Look at disease snps."
    arguments = ["intersect_snps", "get_relative_positions", "get_snp_status", "get_info", "simulate_ptc_location", "get_possible_ptc_locations", "required_simulations"]
    args = gen.parse_arguments(description, arguments, flags = [0,1,2,3,4,5])
    intersect_snps, get_relative_positions, get_snp_status, get_info, simulate_ptc_location, get_possible_ptc_locations, required_simulations = args.intersect_snps, args.get_relative_positions, args.get_snp_status, args.get_info, args.simulate_ptc_location, args.get_possible_ptc_locations, args.required_simulations

    results_prefix = "./results/clean_run_2/clean_run"

    output_directory = "results/disease_snps"
    gen.create_output_directories(output_directory)

    disease_snps_file = "./source_data/clinvar_20180429.vcf.gz"
    disease_snps_index_file = "{0}.tbi".format(disease_snps_file)

    if not os.path.isfile(disease_snps_file) or not os.path.isfile(disease_snps_index_file):
        print("\nERROR: Please provide the required disease SNPs files.\n")
        raise Exception

    # intersect the coding exons with the disease snps
    exon_bed = "{0}_coding_exons.bed".format(results_prefix)
    disease_snp_intersect_file_vcf = "{0}/disease_snp_intersect.vcf".format(output_directory)
    disease_snp_intersect_file_bed = "{0}/disease_snp_intersect.bed".format(output_directory)
    if intersect_snps or not os.path.isfile(disease_snp_intersect_file_vcf) or not os.path.isfile(disease_snp_intersect_file_bed):
        print("Intersecting snps with exons")
        # so.intersect_snps_parallel(exon_bed, disease_snps_file, disease_snp_intersect_file_vcf)
        so.intersect_vcf_to_bed(exon_bed, disease_snp_intersect_file_vcf, disease_snp_intersect_file_bed, change_names = True)

    # get relative positions of the snps in cds and exons
    full_bed = "{0}_CDS.bed".format(results_prefix)
    disease_snps_relative_exon_positions = "{0}/disease_snp_relative_exon_positions.bed".format(output_directory)
    disease_snps_relative_cds_positions = "{0}/disease_snp_relative_cds_positions.bed".format(output_directory)
    if get_relative_positions or not os.path.isfile(disease_snps_relative_exon_positions) or not os.path.isfile(disease_snps_relative_cds_positions):
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(disease_snp_intersect_file_bed, disease_snps_relative_exon_positions)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, disease_snps_relative_cds_positions, full_bed)

    # get the change status of the snps to check them
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    disease_ptcs_file = "{0}/disease_ptcs.txt".format(output_directory)
    disease_other_file = "{0}/disease_other_snps.txt".format(output_directory)
    if get_snp_status or not os.path.isfile(disease_ptcs_file) or not os.path.isfile(disease_other_file) or not os.path.isfile(cds_fasta):
        print("Getting snp status...")
        so.get_snp_change_status(disease_snps_relative_cds_positions, cds_fasta, disease_ptcs_file, disease_other_file)

    # get the info about the ptcs
    output_file_ptc_info = "{0}/disease__analysis_ptc_info.txt".format(output_directory)
    output_file_other_info = "{0}/disease__analysis_other_info.txt".format(output_directory)
    if get_info:
        print("Getting PTC information...")
        get_ptc_info(disease_ptcs_file, disease_snps_relative_exon_positions, output_file_ptc_info)
        get_ptc_info(disease_other_file, disease_snps_relative_exon_positions, output_file_other_info)

    # simulation to see if disease ptcs occur at exon ends more conmonly than by chance
    location_simulation_output_directory = "{0}/ptc_location_simulation".format(output_directory)
    coding_exons_file = "{0}_exons.bed".format(results_prefix)
    if simulate_ptc_location:
        if not required_simulations:
            print("\nERROR: please specify the number of simulations required.\n")
            # create the simulation_folder
        gen.create_output_directories(output_directory)

        if get_possible_ptc_locations:
            print("Getting possible PTC mutation locations...")
            generate_possible_ptc_locations(full_bed, cds_fasta, location_simulation_output_directory)
        ptc_location_simulation(disease_ptcs_file, full_bed, cds_fasta, location_simulation_output_directory, required_simulations, coding_exons_file)


if __name__ == "__main__":
    main()
