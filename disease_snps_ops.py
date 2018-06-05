import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os
import collections
import re
import copy
import numpy as np

def get_gene_length(cds_bed_file):
    '''
    Get info on the whole gene length
    '''
    cds_bed_lines = gen.read_many_fields(cds_bed_file, "\t")
    info = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict()))))

    for line in cds_bed_lines:
        chr = line[0]
        t = line[3].split('.')[0]
        e = int(line[3].split('.')[1])
        type = line[4]
        start = int(line[1])
        stop = int(line[2])
        strand = line[5]

        if chr not in ["chrX", "chrY"]:
            info[strand][t][e][type] = [start, stop]

    transcripts = collections.defaultdict(lambda: [])

    for strand in info:
        for t in info[strand]:
            if strand == "+":
                start_exon = min(sorted(info[strand][t]))
                end_exon = max(sorted(info[strand][t]))
                gene_start = info[strand][t][start_exon]["CDS"][0]
                if "stop_codon" in info[strand][t][end_exon]:
                    gene_end = info[strand][t][end_exon]["stop_codon"][1]
                else:
                    gene_end = info[strand][t][end_exon]["CDS"][1]
            else:
                start_exon = min(sorted(info[strand][t]))
                end_exon = max(sorted(info[strand][t]))
                gene_end = info[strand][t][start_exon]["CDS"][1]
                if "stop_codon" in info[strand][t][end_exon]:
                    gene_start = info[strand][t][end_exon]["stop_codon"][0]
                else:
                    gene_start = info[strand][t][end_exon]["CDS"][0]

            transcripts[t] = [gene_start, gene_end, gene_end-gene_start, strand]

    return(transcripts)


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

    gen.create_output_directories(output_directory)

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

def compare_ptcs(intersect_file, ptc_file, relative_exon_positions_file, exon_fasta, cds_fasta, cds_bed_file, intron_bed, output_file):
    '''
    Get the location of ptcs and whether or not they are disease causing
    '''

    introns = gen.read_many_fields(intron_bed, "\t")
    intron_list = collections.defaultdict(lambda: collections.defaultdict(lambda: [False, False]))

    for intron in introns:
        t = intron[3].split('.')[0]
        exon_3_prime = int(intron[3].split('.')[1].split('-')[0])
        exon_5_prime = int(intron[3].split('.')[1].split('-')[1])
        start = int(intron[1])
        stop = int(intron[2])
        intron_list[t][exon_3_prime][1] = stop-start
        intron_list[t][exon_5_prime][0] = stop-start

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

    gene_lengths = get_gene_length(cds_bed_file)


    ptcs = gen.read_many_fields(ptc_file, "\t")
    with open(output_file, "w") as outfile:
        outfile.write("snp_id\tsnp_pos\ttranscript\texon_no\ttotal_exons\texon_length\trel_exon_position\texon_dist\texon_end\trel_cds_position\tcds_length\tgene_length\tgene_left_length\tgene_right_length\tgene_half\tref_allele\tmut_allele\tdisease\tnmd\tflanking_intron\tother_intron\n")

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

            gene = gene_lengths[transcript]
            gene_start = gene[0]
            gene_end = gene[1]
            gene_middle = gene_start + ((gene_end - gene_start)/2)
            gene_left_length = snp_pos - gene_start
            gene_right_length = gene_end - snp_pos
            if gene_middle - snp_pos > 0:
                gene_half = 5
            else:
                gene_half = 3

            if distances.index(min_dist) == 0:
                end = 5
            else:
                end = 3

            if end == 5:
                flanking_intron = intron_list[transcript][exon][0]
                other_intron = intron_list[transcript][exon][1]
            else:
                flanking_intron = intron_list[transcript][exon][1]
                other_intron = intron_list[transcript][exon][0]

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

            outlist = gen.stringify([snp_id, snp_pos, transcript, exon, total_exons, exon_length, exon_rel_pos, min_dist, end, cds_pos, cds_length, gene[2], gene_left_length, gene_right_length, gene_half, ref_allele, mut_allele, disease, nmd, flanking_intron, other_intron])
            outfile.write("{0}\n".format("\t".join(outlist)))

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
