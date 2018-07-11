import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os
import collections
import re
import copy
import numpy as np
from scipy.stats import chisquare

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

def clinvar_ptc_locations(disease_ptcs_file, disease_snps_relative_exon_positions, ptc_file, output_file):
    '''
    Get the relative positions of disease PTCs.
    Use Chi square test to see whether disease PTCs are located nearer
    exon ends more often than expected given nullself.
    Null takes the number of base pairs in the window as a proportion of exon length.
    e.g. if 5% of base pairs are within 0-3 bp of an exon end, the null is that 5% of
    PTCs should occur in this region.
    '''
    # Get a list of 1000 genomes ptcs
    ptcs = gen.read_many_fields(ptc_file, "\t")
    ptc_list = {}
    ptc_positions = []
    for ptc in ptcs[1:]:
        ptc_id = ptc[8]
        ptc_positions.append(int(ptc[7]))
        ptc_list[ptc_id] = [ptc[3], int(ptc[1]), int(ptc[2]), int(ptc[7])]

    # Get a list of clinvar ptcs
    disease_ptcs = gen.read_many_fields(disease_ptcs_file, "\t")
    disease_ptc_list = {}
    disease_ptc_positions = []
    for ptc in disease_ptcs:
        ptc_id = ptc[8]
        disease_ptc_positions.append(int(ptc[7]))
        disease_ptc_list[ptc_id] = [ptc[3], int(ptc[1]), int(ptc[2]), int(ptc[7])]

    # get all ptcs that overlap between the two datasets
    overlaps = [i for i in ptc_positions if i in disease_ptc_positions]

    # get the relative exon positions of the clinvar ptcs
    disease_relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
    disease_relative_positions_list = {}
    for snp in disease_relative_positions:
        snp_id = snp[8]
        disease_relative_positions_list[snp_id] = [snp[3], int(snp[11])]

    disease_nt_counts = {}
    disease_ptc_counts = {}
    disease_exon_sample_list = []
    ends = [5, 3]
    for end in ends:
        disease_nt_counts[end] = [0,0,0]
        disease_ptc_counts[end] = [0,0,0]

    # for each of the clinvar ptcs not in the overlap, get the number of bp
    # in each window of the exon, the number of ptcs in that window.
    # only sample the number of bp for each exon once (some exons may have
    # more than 1 ptc and would otherwise contribute 2x)

    for ptc_id in disease_ptc_list:
        ptc = disease_ptc_list[ptc_id]
        ptc_pos = ptc[3]
        if ptc_pos not in overlaps:
            if ptc_id in disease_relative_positions_list and disease_relative_positions_list[ptc_id][0] == ptc[0]:
                exon_length = ptc[2] - ptc[1]
                transcript_exon = ptc[0]

                if transcript_exon not in disease_exon_sample_list:
                    if exon_length >= 138:
                        for i in ends:
                            disease_nt_counts[i][0] += 3
                            disease_nt_counts[i][1] += 66
                            disease_nt_counts[i][2] += (np.divide(exon_length, 2) - 69)
                    elif exon_length > 6 and exon_length < 138:
                        for i in ends:
                            disease_nt_counts[i][0] += 3
                            disease_nt_counts[i][1] += (np.divide(exon_length, 2) - 3)
                    else:
                        for i in ends:
                            disease_nt_counts[i][0] += np.divide(exon_length, 2)

                disease_exon_sample_list.append(transcript_exon)

                rel_pos = int(disease_relative_positions_list[ptc_id][1])
                distances = [rel_pos, exon_length - rel_pos]
                min_dist = min(distances)
                if distances.index(min_dist) == 0:
                    end = 5
                else:
                    end = 3

                if min_dist <= 2:
                    print(ptc_id)
                    disease_ptc_counts[end][0] += 1
                elif min_dist > 2 and min_dist <= 68:
                    disease_ptc_counts[end][1] += 1
                else:
                    disease_ptc_counts[end][2] += 1

    # write the output and chisq test to file
    intervals = ["0-3 bp", "4-69 bp", "70+ bp"]

    for end in disease_ptc_counts:
        print(end, disease_ptc_counts[end])

    total_disease_nts = []
    [total_disease_nts.extend(disease_nt_counts[end]) for end in disease_nt_counts]
    total_disease_nts = sum(total_disease_nts)
    total_disease_ptcs = []
    [total_disease_ptcs.extend(disease_ptc_counts[end]) for end in disease_ptc_counts]
    total_disease_ptcs = sum(total_disease_ptcs)

    with open(output_file, "w") as outfile:
        outfile.write('whole_exons\n')
        outfile.write('\n \n')
        outfile.write('region,nt,nt_prop,expected_disease_ptcs,observed_disease_ptcs,fo/fe\n')
        chisq = []
        expecteds = []
        nt_props = []
        for i, region in enumerate(intervals):
            region_nts = sum(disease_nt_counts[end][i] for end in disease_nt_counts)
            region_nts_prop = np.divide(region_nts, total_disease_nts)
            nt_props.append(region_nts_prop)
            expected_region_ptcs = region_nts_prop * total_disease_ptcs
            expecteds.append(expected_region_ptcs)
            region_ptcs = sum(disease_ptc_counts[end][i] for end in disease_ptc_counts)
            outfile.write("{0},{1},{2},{3},{4},{5}\n".format(region, region_nts, region_nts_prop, expected_region_ptcs, region_ptcs, np.divide(region_ptcs, expected_region_ptcs)))
            chisq.append(np.divide((region_ptcs - expected_region_ptcs)**2, expected_region_ptcs))

        observeds = [0,0,0]
        for end in disease_ptc_counts:
            for i, count in enumerate(disease_ptc_counts[end]):
                observeds[i] += count

        chisq_calc = chisquare(observeds, f_exp = expecteds)

        outfile.write("total,{0},{1},{2}\n".format(total_disease_nts, sum(nt_props), sum(expecteds)))
        outfile.write('\n\n')
        outfile.write('total_ptcs:,{0}\n'.format(total_disease_ptcs))
        outfile.write('chisq:,{0}\ndf:,{1}\npval:,{2}\n'.format(sum(chisq), len(chisq) - 1, chisq_calc.pvalue))

        for end in ends:
            outfile.write("\n \n{0}'\n".format(end))
            outfile.write('\n \n')
            outfile.write('region,nt,nt_prop,expected_disease_ptcs,observed_disease_ptcs,fo/fe\n')
            chisq = []
            expecteds = []
            nt_props = []
            for i, region in enumerate(intervals):
                region_nts = disease_nt_counts[end][i]
                region_nts_prop = np.divide(region_nts, np.divide(total_disease_nts, 2))
                nt_props.append(region_nts_prop)
                expected_region_ptcs = region_nts_prop * np.divide(total_disease_ptcs, 2)
                expecteds.append(expected_region_ptcs)
                region_ptcs = disease_ptc_counts[end][i]
                outfile.write("{0},{1},{2},{3},{4},{5}\n".format(region, region_nts, region_nts_prop, expected_region_ptcs, region_ptcs, np.divide(region_ptcs, expected_region_ptcs)))
                chisq.append(np.divide((region_ptcs - expected_region_ptcs)**2, expected_region_ptcs))

            observeds = [0,0,0]
            for i, count in enumerate(disease_ptc_counts[end]):
                observeds[i] += count
            chisq_calc = chisquare(observeds, f_exp = expecteds)

            outfile.write("total,{0},{1},{2}\n".format(sum(disease_nt_counts[end]), sum(nt_props), sum(expecteds)))
            outfile.write('\n\n')
            outfile.write('total_ptcs:,{0}\n'.format(sum(disease_ptc_counts[end])))
            outfile.write('chisq:,{0}\ndf:,{1}\npval:,{2}\n'.format(sum(chisq), len(chisq) - 1, chisq_calc.pvalue))

def get_coding_exon_nt_positions(coding_exons_list, clean_disease_ptc_list, disease_relative_positions_list, exclude_cpg=None):
    '''
    Get the index of each nt in each coding exon that isnt in the ptc list
    '''

    # get the locations of all ptcs grouped by transcript, exon, nt
    relative_ptc_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for ptc_id in clean_disease_ptc_list:
        if ptc_id in disease_relative_positions_list:
            ptc = clean_disease_ptc_list[ptc_id]
            t = ptc[4]
            e = ptc[5]
            strand = ptc[6]
            ref_allele = ptc[7]
            if strand == "-":
                ref_allele = gen.reverse_complement(ref_allele)
            rel_pos = disease_relative_positions_list[ptc_id][1]
            relative_ptc_locations[t][e][ref_allele].append(rel_pos)

    coding_exon_nt_positions = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: [])))
    for transcript in coding_exons_list:
        for exon in coding_exons_list[transcript]:
            exon_seq = list(coding_exons_list[transcript][exon])
            for i, nt in enumerate(exon_seq):
                if i not in relative_ptc_locations[transcript][exon][nt]:
                    # if we want to exclude cpg regions
                    if exclude_cpg and i < len(exon_seq)-1:
                        if nt == "G" and exon_seq[i+1] != "C" or nt == "C" and exon_seq[i+1] != "G":
                            coding_exon_nt_positions[transcript][exon][nt].append(i)
                    else:
                        coding_exon_nt_positions[transcript][exon][nt].append(i)

    return coding_exon_nt_positions


def get_real_positions(clean_ptc_list, relative_positions_list, regions):

    location_counts = {}
    ends = [5,3]
    for region in regions:
        location_counts[region] = {}
        for end in ends:
            location_counts[region][end] = 0

    for ptc_id in clean_ptc_list:
        ptc = clean_ptc_list[ptc_id]
        rel_pos_info = relative_positions_list[ptc_id]

        exon_start = int(ptc[1])
        exon_end = int(ptc[2])
        exon_length = exon_end-exon_start
        rel_pos = rel_pos_info[1]

        distances = [rel_pos, exon_length - rel_pos]
        min_dist = min(distances)
        if distances.index(min_dist) == 0:
            end = 5
        else:
            end = 3

        if min_dist <= 1:
            location_counts[regions[0]][end] += 1
        elif min_dist > 1 and min_dist <= 68:
            location_counts[regions[1]][end] += 1
        else:
            location_counts[regions[2]][end] += 1

    return location_counts

def get_possible_eses(exon_seq, rel_pos):

    possible_eses = []
    if rel_pos - 5 >= 0:
        possible_eses.append(exon_seq[rel_pos-5:rel_pos+1])
    if rel_pos - 4 >= 0 and rel_pos + 1 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-4:rel_pos+2])
    if rel_pos - 3 >= 0 and rel_pos + 2 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-3:rel_pos+3])
    if rel_pos - 2 >= 0 and rel_pos + 3 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-2:rel_pos+4])
    if rel_pos - 1 >= 0 and rel_pos + 4 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos-1:rel_pos+5])
    if rel_pos >= 0 and rel_pos + 5 <= len(exon_seq):
        possible_eses.append(exon_seq[rel_pos:rel_pos+6])

    return possible_eses

def get_real_ese_overlap(clean_ptc_list, relative_positions_list, regions, coding_exons, ese_list):

    ese_overlaps = {}
    ends = [5,3]
    for region in regions:
        ese_overlaps[region] = {}
        for end in ends:
            ese_overlaps[region][end] = 0

    for ptc_id in clean_ptc_list:
        ptc = clean_ptc_list[ptc_id]
        rel_pos_info = relative_positions_list[ptc_id]
        transcript = ptc[4]
        exon = ptc[5]
        exon_seq = coding_exons[transcript][exon]
        rel_pos = rel_pos_info[1]
        exon_length = len(exon_seq)

        distances = [rel_pos, exon_length - rel_pos]
        min_dist = min(distances)
        if distances.index(min_dist) == 0:
            end = 5
        else:
            end = 3

        possible_eses = get_possible_eses(exon_seq, rel_pos)
        ese_overlap = list(set(ese_list) & set(possible_eses))

        if min_dist <= 1 and len(ese_overlap):
            ese_overlaps[regions[0]][end] += 1
        elif min_dist > 1 and min_dist <= 68 and len(ese_overlap):
            ese_overlaps[regions[1]][end] += 1
        elif len(ese_overlap):
            ese_overlaps[regions[2]][end] += 1

    return ese_overlaps

def get_coding_exons(coding_exons_fasta):

    coding_exons = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    exon_names, exon_seqs = gen.read_fasta(coding_exons_fasta)

    for i, name in enumerate(exon_names):
        transcript = name.split('(')[0].split('.')[0]
        exon = int(name.split('(')[0].split('.')[1])
        seq = exon_seqs[i]
        coding_exons[transcript][exon] = seq

    return coding_exons


def clinvar_simulation(disease_ptcs_file, relative_exon_positions_file, ptc_file, coding_exons_fasta, simulations, output_file, ese_overlap_output_file, ese_file=None, exclude_cpg=None, clinvar=None):
    '''
    Simulation mutation locations of PTCs.
    Take the exon in which each PTC is location and randomly pick a site with the
    same nt composition.
    Locations of these matched mutations are used for null.
    '''

    # get the clean ptcs that appear only in the individual datasets
    clean_genomes_ptc_list, clean_clinvar_ptc_list = get_clean_ptc_list(ptc_file, disease_ptcs_file)
    if clinvar:
        clean_ptc_list = clean_clinvar_ptc_list
    else:
        clean_ptc_list = clean_genomes_ptc_list

    # get the relative exon positions of the clinvar snps
    relative_positions = gen.read_many_fields(relative_exon_positions_file, "\t")
    relative_positions_list = {}
    for snp in relative_positions:
        snp_id = snp[8]
        relative_positions_list[snp_id] = [snp[3], int(snp[11])]

    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t") if ese[0][0] != "#"]
    coding_exons = get_coding_exons(coding_exons_fasta)

    regions = ["0-2 bp", "3-69 bp", "70+ bp"]
    real_positions = get_real_positions(clean_ptc_list, relative_positions_list, regions)
    real_ese_overlap = get_real_ese_overlap(clean_ptc_list, relative_positions_list, regions, coding_exons, ese_list)

    simulant_list = list(range(1, simulations+1))
    processes = gen.run_in_parallel(simulant_list, ["foo", simulations, clean_ptc_list, relative_positions_list, coding_exons_fasta, ese_list, exclude_cpg], simulate_mutations)

    position_list = {}
    ese_overlap_list = {}
    for process in processes:
        result = process.get()
        position_list = {**position_list, **result[0]}
        ese_overlap_list = {**ese_overlap_list, **result[1]}

    write_to_file(real_positions, position_list, regions, output_file)
    write_to_file(real_ese_overlap, ese_overlap_list, regions, ese_overlap_output_file)

def write_to_file(real_positions, process_list, regions, output_file):

    # write all to file
    out_ends = ["all", 5, 3]
    ends = [5,3]
    header_list = ["sim_id"]
    for end in out_ends:
        for region in regions:
            header_list.append("{0}_{1}".format(end, region))

    with open(output_file, "w") as outfile:
        outfile.write("{0}\n".format(",".join(header_list)))
        outlist = ["real"]
        for region in regions:
            outlist.append(sum(real_positions[region].values()))
        for end in ends:
            for region in regions:
                outlist.append(real_positions[region][end])
        outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))

        for simulation in sorted(process_list):
            outlist = [simulation]
            for region in regions:
                outlist.append(sum(process_list[simulation][region].values()))
            for end in ends:
                for region in regions:
                    outlist.append(process_list[simulation][region][end])
            outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))



def simulate_mutations(simulant_list, simulations, clean_ptc_list, relative_positions_list, coding_exons_fasta, ese_list, exclude_cpg):
    '''
    Run the simulations
    '''

    # This needs to be done in here for parallelisation
    coding_exons = get_coding_exons(coding_exons_fasta)

    coding_exons_names, coding_exons_seqs = gen.read_fasta(coding_exons_fasta)
    coding_exons_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))
    for i, name in enumerate(coding_exons_names):
        name = name.split('(')[0]
        t = name.split('.')[0]
        e = int(name.split('.')[1])
        coding_exons_list[t][e] = coding_exons_seqs[i]

    coding_exon_nt_positions = get_coding_exon_nt_positions(coding_exons_list, clean_ptc_list, relative_positions_list, exclude_cpg)

    simulation_outputs = {}
    simulation_outputs_ese_overlaps = {}

    for simulation in simulant_list:
        np.random.seed()
        print("Running simulation {0}/{1}".format(simulation, simulations))

        location_counts = {}
        ese_overlaps = {}
        regions = ["0-2 bp", "3-69 bp", "70+ bp"]
        ends = [5,3]
        for region in regions:
            location_counts[region] = {}
            ese_overlaps[region] = {}
            for end in ends:
                location_counts[region][end] = 0
                ese_overlaps[region][end] = 0

        for i, ptc_id in enumerate(clean_ptc_list):
            if ptc_id in relative_positions_list:

                ptc = clean_ptc_list[ptc_id]
                pos = ptc[3]
                t = ptc[4]
                e = ptc[5]
                strand = ptc[6]
                ref_allele = ptc[7]
                if strand == "-":
                    gen.reverse_complement(ref_allele)
                exon_length = len(coding_exons_list[t][e])
                rel_pos_info = relative_positions_list[ptc_id]
                rel_pos = rel_pos_info[1]

                positions = coding_exon_nt_positions[t][e][ref_allele]

                if len(positions):
                    simulant_position = np.random.choice(positions, 1)[0]
                else:
                    simulant_position = rel_pos

                distances = [simulant_position, exon_length - simulant_position]
                min_dist = min(distances)
                if distances.index(min_dist) == 0:
                    end = 5
                else:
                    end = 3

                exon_seq = coding_exons[t][e]
                possible_eses = get_possible_eses(exon_seq, simulant_position)
                overlaps = list(set(ese_list) & set(possible_eses))

                if min_dist <= 1:
                    location_counts[regions[0]][end] += 1
                    if len(overlaps):
                        ese_overlaps[regions[0]][end] += 1
                elif min_dist > 1 and min_dist <= 68:
                    location_counts[regions[1]][end] += 1
                    if len(overlaps):
                        ese_overlaps[regions[1]][end] += 1
                else:
                    location_counts[regions[2]][end] += 1
                    if len(overlaps):
                        ese_overlaps[regions[2]][end] += 1


        simulation_outputs[simulation] = location_counts
        simulation_outputs_ese_overlaps[simulation] = ese_overlaps


    return simulation_outputs, simulation_outputs_ese_overlaps

def get_clean_ptc_list(ptc_file, disease_ptcs_file):
    '''
    Get a list of PTCs from the dataset
    '''
    ptcs = gen.read_many_fields(ptc_file, "\t")
    ptc_list = {}
    ptc_positions = []
    for ptc in ptcs[1:]:
        ptc_id = ptc[8]
        ptc_pos = int(ptc[7])
        ptc_positions.append(ptc_pos)
        start = int(ptc[1])
        end = int(ptc[2])
        strand = ptc[5]
        ref_allele = ptc[9]
        name = ptc[3]
        t = ptc[3].split('.')[0]
        e = int(ptc[3].split('.')[1])
        ptc_list[ptc_id] = [name, start, end, ptc_pos, t, e, strand, ref_allele]

    disease_ptcs = gen.read_many_fields(disease_ptcs_file, "\t")
    disease_ptc_list = {}
    disease_ptc_positions = []
    for ptc in disease_ptcs:
        ptc_id = ptc[8]
        ptc_pos = int(ptc[7])
        disease_ptc_positions.append(ptc_pos)
        start = int(ptc[1])
        end = int(ptc[2])
        strand = ptc[5]
        ref_allele = ptc[9]
        name = ptc[3]
        t = ptc[3].split('.')[0]
        e = int(ptc[3].split('.')[1])
        disease_ptc_list[ptc_id] = [name, start, end, ptc_pos, t, e, strand, ref_allele]

    overlaps = [i for i in ptc_positions if i in disease_ptc_positions]

    clean_genomes_ptc_list = {}
    clean_clinvar_ptc_list = {}


    for ptc_id in disease_ptc_list:
        if disease_ptc_list[ptc_id][3] not in overlaps:
            clean_clinvar_ptc_list[ptc_id] = disease_ptc_list[ptc_id]

    for ptc_id in ptc_list:
        if ptc_list[ptc_id][3] not in overlaps:
            clean_genomes_ptc_list[ptc_id] = ptc_list[ptc_id]

    return clean_genomes_ptc_list, clean_clinvar_ptc_list

def get_ptcs_in_window(ptc_list, relative_positions_list, window_start, window_end):
    '''
    Filter PTCs to only those in the given window (inclusive).
    window_start: start of the window (in base 1)
    window_end: end of the window (in base 1)
    '''

    ends = [5,3]
    ptcs_in_window = {}
    for end in ends:
        ptcs_in_window[end] = {}

    for ptc_id in ptc_list:
        rel_pos_info = relative_positions_list[ptc_id]
        ptc = ptc_list[ptc_id]
        exon_length = ptc[2] - ptc[1]

        # need to get exons longer than the window end * 2 to account for
        # both 5' and 3' ends
        if exon_length >= (window_end*2):
            rel_pos = rel_pos_info[1]
            ptc.append(rel_pos)
            if rel_pos >= window_start - 1 and rel_pos <= window_end - 1:
                ptcs_in_window[5][ptc_id] = ptc
            elif rel_pos <= exon_length - window_start and rel_pos >= exon_length - window_end:
                ptcs_in_window[3][ptc_id] = ptc

    return ptcs_in_window

def get_non_mut_nts(ptc_list, coding_exons, window_start, window_end):

    non_mut_nts = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: []))))
    relative_ptc_locations = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: []))))
    exon_lengths = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))

    for end in ptc_list:
        for ptc_id in ptc_list[end]:
            ptc = ptc_list[end][ptc_id]
            transcript = ptc[4]
            exon = ptc[5]
            exon_length = ptc[2] - ptc[1]
            rel_pos = int(ptc[8])
            exon_seq = coding_exons[transcript][exon]
            nt = exon_seq[rel_pos]

            relative_ptc_locations[end][transcript][exon][nt].append(rel_pos)
            exon_lengths[transcript][exon] = exon_length

    # now create a dictionary of all the positions that dont have a mutation
    # in each exon for each nt
    for transcript in exon_lengths:
        for exon in exon_lengths[transcript]:
            exon_nts = coding_exons[transcript][exon]
            start_nts = list(exon_nts[window_start-1:window_end])

            end_start_index = len(exon_nts) - window_end
            end_end_index = len(exon_nts) - window_start + 1
            end_nts = list(exon_nts[end_start_index:end_end_index])

            for i, nt in enumerate(start_nts, start=window_start-1):
                if i not in relative_ptc_locations[5][transcript][exon][nt]:
                    non_mut_nts[transcript][exon][5][nt].append(i)
            for i, nt in enumerate(end_nts, start=len(exon_nts)-window_end):
                if i not in relative_ptc_locations[3][transcript][exon][nt]:
                    non_mut_nts[transcript][exon][3][nt].append(i)

    return non_mut_nts

def get_real_ese_hits(ptc_list, coding_exons, ese_list):

    ends = [5,3]
    hit_counts = [0,0]
    for end in ptc_list:
        for ptc_id in ptc_list[end]:
            ptc = ptc_list[end][ptc_id]
            transcript = ptc[4]
            exon = ptc[5]
            rel_pos = ptc[8]
            exon_seq = coding_exons[transcript][exon]
            possible_eses = get_possible_eses(exon_seq, rel_pos)
            overlaps = list(set(ese_list) & set(possible_eses))
            if len(overlaps):
                hit_counts[ends.index(end)] += 1

    return hit_counts

def simulate_ese_hits(simulation_list, simulations, ptc_list, coding_exons_fasta, ese_list, window_start, window_end):

    ends = [5,3]
    hit_counts = {}

    coding_exons = get_coding_exons(coding_exons_fasta)

    for simulation in simulation_list:
        print("Simulation {0}/{1}".format(simulation+1, simulations))
        simulation_non_mut_nts = get_non_mut_nts(ptc_list, coding_exons, window_start, window_end)
        hit_count = [0,0]
        cant_count = [0,0]
        np.random.seed()
        for end in ptc_list:
            for i, ptc_id in enumerate(ptc_list[end]):
                if i:
                    ptc = ptc_list[end][ptc_id]
                    transcript = ptc[4]
                    exon = ptc[5]
                    rel_pos = ptc[8]
                    exon_seq = coding_exons[transcript][exon]
                    mut_nt = exon_seq[rel_pos]

                    if len(simulation_non_mut_nts[transcript][exon][end][mut_nt]):
                        sim_choice = np.random.choice(simulation_non_mut_nts[transcript][exon][end][mut_nt], replace=False)
                    else:
                        cant_count[ends.index(end)] += 1
                        sim_choice = rel_pos

                    possible_eses = get_possible_eses(exon_seq, sim_choice)
                    overlaps = list(set(ese_list) & set(possible_eses))
                    if len(overlaps):
                        hit_count[ends.index(end)] += 1
            hit_counts[simulation] = [hit_count, cant_count]

    return hit_counts


def ese_hit_simulation(disease_ptcs_file, relative_exon_positions_file, ptc_file, coding_exons_fasta, simulations, output_file, ese_file, window_start, window_end, exclude_cpg, clinvar=None):

    # get the ptcs that occur strictly in one dataset
    clean_genomes_ptc_list, clean_clinvar_ptc_list = get_clean_ptc_list(ptc_file, disease_ptcs_file)

    if clinvar:
        clean_ptc_list = clean_clinvar_ptc_list
    else:
        clean_ptc_list = clean_genomes_ptc_list

    # get the relative exon positions of the clinvar snps
    relative_positions = gen.read_many_fields(relative_exon_positions_file, "\t")
    relative_positions_list = {}
    for snp in relative_positions:
        snp_id = snp[8]
        relative_positions_list[snp_id] = [snp[3], int(snp[11])]

    # get the list of eses and coding exons
    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t") if ese[0][0] != "#"]
    coding_exons = get_coding_exons(coding_exons_fasta)

    # get the ptcs that are in the 4-69 bp region for each exon of exon
    # this requires exons at least 128 bp in length for comparison
    window_ptcs = get_ptcs_in_window(clean_ptc_list, relative_positions_list, window_start, window_end)
    # get indices of each nt in that region that dont contain a ptc mutation
    non_mut_nts = get_non_mut_nts(window_ptcs, coding_exons, window_start, window_end)

    # get the hit counts of eses for the real ptcs
    real_ese_hits = get_real_ese_hits(window_ptcs, coding_exons, ese_list)

    # simulate the hit counts for nt matched mutations
    simulation_list = list(range(simulations))
    # simulate_ese_hits(simulation_list, simulations, window_ptcs, coding_exons_fasta, ese_list, window_start, window_end)
    processes = gen.run_in_parallel(simulation_list, ["foo", simulations, window_ptcs, coding_exons_fasta, ese_list, window_start, window_end], simulate_ese_hits)

    simulation_outputs = {}
    for process in processes:
        simulation_hits = process.get()
        simulation_outputs = {**simulation_outputs, **simulation_hits}

    ends = [5,3]
    with open(output_file, "w") as outfile:
        outfile.write("simulation,5_ese_hit_count,5_cant_count,3_ese_hit_count,3_cant_count,both_ese_hit_count,both_cant_cant\n")
        outfile.write("real,{0},0,{1},0,{2},0\n".format(real_ese_hits[0], real_ese_hits[1], sum(real_ese_hits)))
        for simulation in sorted(simulation_outputs):
            outlist = [simulation+1]
            for end in ends:
                outlist.append(simulation_outputs[simulation][0][ends.index(end)])
                outlist.append(simulation_outputs[simulation][1][ends.index(end)])
            outlist.append(sum(simulation_outputs[simulation][0]))
            outlist.append(sum(simulation_outputs[simulation][1]))
            outfile.write("{0}\n".format(",".join(gen.stringify(outlist))))
