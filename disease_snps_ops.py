import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os
import collections
import re
import copy
import numpy as np

def get_ptc_overlaps(disease_ptc_file, ptc_file, output_file):
    '''
    Get the overlap between disease ptc file and 1000 genomes ptcs
    '''

    disease_ptcs = gen.read_many_fields(disease_ptc_file, "\t")
    ptcs = gen.read_many_fields(ptc_file, "\t")

    disease_ptc_list = {}
    for ptc in disease_ptcs[1:]:
        ptc_info = EntryInfo(ptc)
        disease_ptc_list[ptc_info.snp_pos] = ptc

    ptc_list = {}
    for ptc in ptcs[1:]:
        ptc_list[int(ptc[7])] = ptc

    with open(output_file, "w") as outfile:
        for ptc_pos in ptc_list:
            if ptc_pos in disease_ptc_list:
                outfile.write('{0}\t**\t{1}\n'.format("\t".join(ptc_list[ptc_pos]), "\t".join(disease_ptc_list[ptc_pos])))


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

def refactor_ptc_file(input_file, output_file, header=None):
    '''
    Refactor PTC file ready for intersect
    '''
    entries = gen.read_many_fields(input_file, "\t")
    if header:
        start = 1
    else:
        start = 0

    with open(output_file, "w") as outfile:
        for entry in entries[start:]:
            exon_info = entry[:6]
            ptc_info = entry[6:]
            ptc_end = str(int(ptc_info[1]) + 1)
            strand = exon_info[5]
            ptc_info.insert(2, ptc_end)
            ptc_info.insert(3, ".")
            ptc_info.insert(4, strand)
            ptc_info.extend(exon_info)
            ptc_info[0] = ptc_info[0].strip('chr')
            outfile.write("{0}\n".format("\t".join(ptc_info)))

def compare_ptcs(intersect_file, ptc_file, relative_exon_positions_file, exon_fasta, cds_fasta, output_file):
    '''
    Get the location of ptcs and whether or not they are disease causing
    '''
    exon_list = collections.defaultdict(lambda: [])
    exon_names, exon_seqs = gen.read_fasta(exon_fasta)
    for i, name in enumerate(exon_names):
        transcript = name.split('.')[0]
        exon = int(name.split('.')[1])
        if len(exon_seqs[i]) > 3:
            exon_list[transcript].append(exon)

    cds_list = {}
    cds_names, cds_seqs = gen.read_fasta(cds_fasta)
    for i, name in enumerate(cds_names):
        cds_list[name] = cds_seqs[i]

    disease_ptcs = gen.read_many_fields(intersect_file, "\t")
    disease_ptc_ids = []
    for ptc in disease_ptcs:
        disease_ptc_ids.append(int(ptc[11]))

    relative_position_entries = gen.read_many_fields(relative_exon_positions_file, "\t")
    relative_positions = {}
    for entry in relative_position_entries:
        entry = entry[:12]
        snp_id = entry[8]
        exon_length = int(entry[2]) - int(entry[1])
        rel_pos = entry[-1]
        relative_positions[snp_id] = [entry[3], int(entry[7]), int(rel_pos), exon_length]


    ptcs = gen.read_many_fields(ptc_file, "\t")
    with open(output_file, "w") as outfile:
        outfile.write("snp_id\ttranscript\texon_no\ttotal_exons\texon_length\trel_exon_position\texon_dist\texon_end\trel_cds_position\tcds_length\tref_allele\tmut_allele\tdisease\tnmd\n")

        for ptc in ptcs[1:]:
            snp_pos = int(ptc[7])
            snp_id = ptc[8]
            transcript_exon = ptc[3]
            transcript = ptc[3].split('.')[0]
            exon = int(ptc[3].split('.')[1])
            snp_ref = int(ptc[14])
            cds_pos = ptc[11]
            ref_allele = ptc[9]
            mut_allele = ptc[10]

            relative_position_entry = relative_positions[snp_id]
            exon_rel_pos = relative_position_entry[2]
            exon_length = relative_position_entry[3]
            dist_to_end = exon_length - exon_rel_pos

            distances = [exon_rel_pos, dist_to_end]
            min_dist = min(distances)

            if distances.index(min_dist) == 0:
                end = 5
            else:
                end = 3

            total_exons = len(exon_list[transcript])
            cds_length = len(cds_list[transcript])

            if exon == sorted(exon_list[transcript])[-1] or exon == sorted(exon_list[transcript])[-2] and dist_to_end < 50:
                nmd = 0
            else:
                nmd = 1

            if snp_ref in disease_ptc_ids:
                disease = 1
            else:
                disease = 0

            outlist = gen.stringify([snp_id, transcript, exon, total_exons, exon_length, exon_rel_pos, min_dist, end, cds_pos, cds_length, ref_allele, mut_allele, disease, nmd])
            outfile.write("{0}\n".format("\t".join(outlist)))


def compare_distances(clinvar_ptc_file, clinvar_relative_exon_positions_file, ptc_file, ptc_relative_exon_positions_file, output_file):

    clinvar_relative_position_entries = gen.read_many_fields(clinvar_relative_exon_positions_file, "\t")
    clinvar_relative_positions = {}
    for entry in clinvar_relative_position_entries:
        # print(entry)
        entry = entry[:12]
        snp_id = entry[8]
        exon_length = int(entry[2]) - int(entry[1])
        rel_pos = entry[-1]
        pos = int(entry[7])
        clinvar_relative_positions[snp_id] = [snp_id, entry[3], int(entry[7]), int(rel_pos), exon_length]

    # print(sorted([pos for pos in clinvar_relative_positions]))

    ptc_relative_position_entries = gen.read_many_fields(ptc_relative_exon_positions_file, "\t")
    ptc_relative_positions = {}
    for entry in ptc_relative_position_entries:
        entry = entry[:12]
        snp_id = entry[8]
        rel_pos = entry[-1]
        pos = int(entry[7])
        exon_length = int(entry[2]) - int(entry[1])
        ptc_relative_positions[pos] = [snp_id, entry[3], int(entry[7]), int(rel_pos), exon_length]

    clinvar_ptcs = gen.read_many_fields(clinvar_ptc_file, "\t")
    clinvar_ptc_locations = []
    clinvar_ptc_list = {}
    for ptc in clinvar_ptcs:
        pos = int(ptc[7])
        clinvar_ptc_locations.append(pos)
        clinvar_ptc_list[pos] = ptc
    ptcs = gen.read_many_fields(ptc_file, "\t")
    ptc_locations = []
    ptc_list = {}
    for ptc in ptcs[1:]:
        pos = int(ptc[7])
        ptc_locations.append(pos)
        ptc_list[pos] = ptc

    overlap = [loc for loc in ptc_locations if loc in clinvar_ptc_locations]

    count = 0
    with open(output_file, "w") as outfile:
        outfile.write("transcript,exon,exon_length,rel_exon_pos,min_dist,exon_end,1000_genomes,clinvar\n")
        for pos in ptc_list:
            if pos not in overlap:
                ptc = ptc_list[pos]
                ptc_id = ptc[8]
                exon = ptc[3]
                t = exon.split('.')[0]
                e = exon.split('.')[1]
                e_start = int(ptc[1])
                e_end = int(ptc[2])
                gpos = int(ptc[7])
                rel_pos_info = ptc_relative_positions[gpos]
                rel_pos = rel_pos_info[3]
                exon_length = rel_pos_info[4]
                clinvar = 0
                kgenomes = 1

                dist_to_end = exon_length - rel_pos

                distances = [rel_pos, dist_to_end]
                min_dist = min(distances)

                if distances.index(min_dist) == 0:
                    exon_end = 5
                else:
                    exon_end = 3

                outlist = gen.stringify([t, e, exon_length, rel_pos, min_dist, exon_end, kgenomes, clinvar])
                outfile.write("{0}\n".format(",".join(outlist)))
        for pos in clinvar_ptc_list:
            if pos not in overlap:
                ptc = clinvar_ptc_list[pos]
                snp_id = ptc[8]
                if snp_id in clinvar_relative_positions:

                    rel_pos_info = clinvar_relative_positions[snp_id]

                    ptc_id = ptc[8]
                    exon = ptc[3]
                    t = exon.split('.')[0]
                    e = exon.split('.')[1]
                    e_start = int(ptc[1])
                    e_end = int(ptc[2])
                    gpos = int(ptc[7])

                    rel_pos = rel_pos_info[3]
                    exon_length = rel_pos_info[4]
                    clinvar = 1
                    kgenomes = 0



                    dist_to_end = exon_length - rel_pos

                    distances = [rel_pos, dist_to_end]
                    min_dist = min(distances)

                    if distances.index(min_dist) == 0:
                        exon_end = 5
                    else:
                        exon_end = 3

                    outlist = gen.stringify([t, e, exon_length, rel_pos, min_dist, exon_end, kgenomes, clinvar])
                    outfile.write("{0}\n".format(",".join(outlist)))

def get_introns(exon_bed, output_file):
    '''
    Get the introns between exons in a file
    '''

    class Define_Exon(object):
        def __init__(self, exon):
            self.chr = exon[0]
            self.start = int(exon[1])
            self.stop = int(exon[2])
            self.transcript_id = exon[3].split('.')[0]
            self.exon_no = int(exon[3].split('.')[1])
            self.type = exon[4]
            self.strand = exon[5]


    exons = gen.read_many_fields(exon_bed, "\t")

    # get a dictionary of exons split by transcript and number
    exon_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for item in exons:
        exon = Define_Exon(item)
        if exon.type == "stop_codon":
            exon.exon_no = 999999
        exon_list[exon.transcript_id][exon.exon_no] = item

    # now get the introns and write to file
    with open(output_file, "w") as outfile:
        for transcript in exon_list:
            for exon_no in sorted(exon_list[transcript]):
                exon = Define_Exon(exon_list[transcript][exon_no])
                # check that the next exon exists, assuming its id will not be higher than 999999
                if exon.exon_no + 1 in exon_list[transcript]:
                    next_exon = Define_Exon(exon_list[exon.transcript_id][exon.exon_no + 1])

                    if exon.strand == "-":
                        intron_start = next_exon.stop
                        intron_stop = exon.start
                    else:
                        intron_start = exon.stop
                        intron_stop = next_exon.start

                    outlist = gen.stringify([exon.chr, intron_start, intron_stop, "{0}.{1}-{2}".format(exon.transcript_id, exon.exon_no, next_exon.exon_no), ".", exon.strand])
                    outfile.write("{0}\n".format("\t".join(outlist)))

def distance_intervals(coding_exons_file, disease_ptcs_file, disease_snps_relative_exon_positions, ptc_file, relative_exon_positions_file, output_file):

    # coding_exons = gen.read_many_fields(coding_exons_file, "\t")
    # coding_exon_list = collections.defaultdict(lambda: collections.defaultdict())
    # for i in coding_exons:
    #     t = i[3].split('.')[0]
    #     e = int(i[3].split('.')[1])
    #     coding_exon_list[t][e] = [i[-1], int(i[1]), int(i[2])]

    ptcs = gen.read_many_fields(ptc_file, "\t")
    ptc_list = {}
    ptc_positions = []
    for ptc in ptcs[1:]:
        ptc_id = ptc[8]
        ptc_positions.append(int(ptc[7]))
        ptc_list[ptc_id] = [ptc[3], int(ptc[1]), int(ptc[2]), int(ptc[7])]

    disease_ptcs = gen.read_many_fields(disease_ptcs_file, "\t")
    disease_ptc_list = {}
    disease_ptc_positions = []
    for ptc in disease_ptcs:
        ptc_id = ptc[8]
        disease_ptc_positions.append(int(ptc[7]))
        disease_ptc_list[ptc_id] = [ptc[3], int(ptc[1]), int(ptc[2]), int(ptc[7])]

    overlaps = [i for i in ptc_positions if i in disease_ptc_positions]

    # relative_positions = gen.read_many_fields(relative_exon_positions_file, "\t")
    # relative_positions_list = {}
    # for snp in relative_positions:
    #     snp_id = snp[8]
    #     relative_positions_list[snp_id] = [snp[3], int(snp[11])]
    #
    # non_disease_nt_counts = [0,0,0]
    # non_disease_ptc_counts = [0,0,0]
    # non_disease_ptc_counts_end = {}
    # for end in [5, 3]:
    #     non_disease_ptc_counts_end = [0,0,0]
    # non_disease_count = 0
    # non_disease_sample_list = []
    #
    # for ptc_id in ptc_list:
    #     ptc = ptc_list[ptc_id]
    #     exon_length = ptc[2] - ptc[1]
    #     transcript_exon = ptc[0]
    #     exon_start = ptc[1]
    #     exon_end = ptc[2]
    #     snp_pos = ptc[3]
    #     if snp_pos not in overlaps:
    #         if ptc_id in relative_positions_list and relative_positions_list[ptc_id][0] == ptc[0]:
    #             non_disease_count += 1
    #             rel_pos = relative_positions_list[ptc_id][1]
    #
    #             distances = [rel_pos, exon_length - rel_pos]
    #             min_dist = min(distances)
    #             if distances.index(min_dist) == 0:
    #                 end = 5
    #             else:
    #                 end = 3
    #
    #             if rel_pos <= 3:
    #                 non_disease_ptc_counts[0] += 1
    #                 non_disease_ptc_counts_end[end][0] += 1
    #             elif rel_pos > 3 and rel_pos <= 69:
    #                 non_disease_ptc_counts[1] += 1
    #                 non_disease_ptc_counts_end[end][1] += 1
    #             else:
    #                 non_disease_ptc_counts[2] += 1
    #                 non_disease_ptc_counts_end[end][2] += 1
    #
    #             if transcript_exon not in non_disease_sample_list:
    #                 if exon_length >= 138:
    #                     non_disease_nt_counts[0] += 6
    #                     non_disease_nt_counts[1] += 132
    #                     non_disease_nt_counts[2] += (exon_length - 188)
    #                 elif exon_length > 6 and exon_length < 138:
    #                     non_disease_nt_counts[0] += 6
    #                     non_disease_nt_counts[1] += (exon_length - 6)
    #                 else:
    #                     non_disease_nt_counts[0] += exon_length
    #
    #                 non_disease_sample_list.append(transcript_exon)


    disease_relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
    disease_relative_positions_list = {}
    for snp in disease_relative_positions:
        snp_id = snp[8]
        disease_relative_positions_list[snp_id] = [snp[3], int(snp[11])]

    disease_nt_counts = {}
    disease_ptc_counts = {}
    disease_exon_sample_list = []
    ranges = [5, 3]
    for i in ranges:
        disease_nt_counts[i] = [0,0,0]
        disease_ptc_counts[i] = [0,0,0]



    for ptc_id in disease_ptc_list:
        ptc = disease_ptc_list[ptc_id]
        ptc_pos = ptc[3]
        if ptc_pos not in overlaps:
            if ptc_id in disease_relative_positions_list and disease_relative_positions_list[ptc_id][0] == ptc[0]:
                exon_length = ptc[2] - ptc[1]
                transcript_exon = ptc[0]

                if transcript_exon not in disease_exon_sample_list:
                    if exon_length >= 138:
                        for i in ranges:
                            disease_nt_counts[i][0] += 3
                            disease_nt_counts[i][1] += 66
                            disease_nt_counts[i][2] += (np.divide(exon_length, 2) - 69)
                    elif exon_length > 6 and exon_length < 138:
                        for i in ranges:
                            disease_nt_counts[i][0] += 3
                            disease_nt_counts[i][1] += (np.divide(exon_length, 2) - 3)
                    else:
                        for i in ranges:
                            disease_nt_counts[i][0] += np.divide(exon_length, 2)

                disease_exon_sample_list.append(transcript_exon)

                rel_pos = disease_relative_positions_list[ptc_id][1]
                distances = [rel_pos, exon_length - rel_pos]
                min_dist = min(distances)
                if distances.index(min_dist) == 0:
                    end = 5
                else:
                    end = 3

                if rel_pos <= 3:
                    disease_ptc_counts[end][0] += 1
                elif rel_pos > 3 and rel_pos <= 69:
                    disease_ptc_counts[end][1] += 1
                else:
                    disease_ptc_counts[end][2] += 1


                # rel_pos = disease_relative_positions_list[ptc_id][1]
                # distances = [rel_pos, exon_length - rel_pos]
                # min_dist = min(distances)
                # if distances.index(min_dist) == 0:
                #     end = 5
                # else:
                #     end = 3


                # if rel_pos <= 3:
                #     disease_ptc_counts[0] += 1
                #     disease_ptc_counts_end[end][0] += 1
                # elif rel_pos > 3 and rel_pos <= 69:
                #     disease_ptc_counts[1] += 1
                #     disease_ptc_counts_end[end][1] += 1
                # else:
                #     disease_ptc_counts[2] += 1
                #     disease_ptc_counts_end[end][2] += 1
                #


    print(disease_nt_counts)


    # expected_props = []
    # for i in range(len(non_disease_ptc_counts)):
    #     expected_props.append(np.divide(non_disease_ptc_counts[i], sum(non_disease_ptc_counts)))
    #
    # ranges = ["0 - 3", "4 - 69", "69 - "]
    # fofe = {}
    # chisq = []
    # for i, interval in enumerate(ranges):
    #     expected = expected_props[i]*sum(disease_ptc_counts)
    #     observed = disease_ptc_counts[i]
    #     fofe[interval] = np.divide(observed, expected)
    #     chisq.append(np.divide((observed - expected)**2, expected))
    #
    # with open(output_file, "w") as outfile:
    #     outfile.write('non_disease\n')
    #     outfile.write('position,num_bases,num_ptc,ptc_prop\n')
    #     outfile.write('0 - 3,{0},{1},{2}\n'.format(non_disease_nt_counts[0], non_disease_ptc_counts[0], expected_props[0]))
    #     outfile.write('4 - 69,{0},{1},{2}\n'.format(non_disease_nt_counts[1], non_disease_ptc_counts[1], expected_props[1]))
    #     outfile.write('70 -,{0},{1},{2}\n'.format(non_disease_nt_counts[2], non_disease_ptc_counts[2], expected_props[2]))
    #     outfile.write('\n')
    #     outfile.write('disease\n')
    #     outfile.write(',num_bases,expected,observed\n')
    #     outfile.write(',,(non_disease_prop * num_disease_ptcs),\n')
    #     outfile.write('0 - 3,{0},{1},{2}\n'.format(disease_nt_counts[0], expected_props[0]*sum(disease_ptc_counts), disease_ptc_counts[0]))
    #     outfile.write('4 - 69,{0},{1},{2}\n'.format(disease_nt_counts[1], expected_props[1]*sum(disease_ptc_counts), disease_ptc_counts[1]))
    #     outfile.write('79 -,{0},{1},{2}\n'.format(disease_nt_counts[2], expected_props[2]*sum(disease_ptc_counts), disease_ptc_counts[2]))
    #     outfile.write('\n')
    #     outfile.write('Chisq = {0},df = {1}\n'.format(sum(chisq), len(chisq) - 1))
    #     outfile.write('\n\n')
    #     outfile.write('fo/fe\n')
    #     for interval in ranges:
    #         outfile.write('{0},{1}\n'.format(interval, fofe[interval]))
    #
    #
    #     for end in [5, 3]:
    #         outfile.write("non_disease {end}'".format(end))
    #         outfile.write('position,num_ptc,ptc_prop\n')
    #         outfile.write('0-3,{0},{1}\n'.format(non_disease_ptc_counts_end[end][0], expected_props_end[end][0]))
    #         outfile.write('4-69,{0},{1}\n'.format(non_disease_ptc_counts_end[end][1], expected_props_end[end][1]))
    #         outfile.write('70-,{0},{1}\n'.format(non_disease_ptc_counts_end[end][2], expected_props_end[end][2]))
