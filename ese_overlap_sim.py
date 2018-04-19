import bed_ops as bo
import generic as gen
import bam_ops as ba
import SNP_ops as so
import numpy as np
import generic as gen
import re
import collections
import time
import os

# in files
# ptc_file = "./results/clean_run_2/clean_run_ptc_file.txt"
syn_nonsyn_file = "./results/clean_run_2/clean_run_syn_nonsyn_file.txt"
snp_rel_pos_file = "./results/clean_run_2/clean_run_SNP_relative_exon_position.bed"
# cds_intervals_file = "./results/clean_run_2/clean_run_CDS_intervals.fasta"
# int3_ese_file = "./ese_data/CaceresHurstESEs_INT3/CaceresHurstESEs_INT3.txt"

# out files
ptc_snps_file = "./temp_data/ptc_snps_file.txt"
main_outfile = "./results/ese_overlap.csv"

stops = ["TAA", "TAG", "TGA"]


class snp_info(object):
    def __init__(self, snp, cds_list):
        self.chr = snp[0]
        self.start = snp[1]
        self.end = snp[2]
        self.transcript_id = snp[3].split('.')[0]
        self.exon_id = int(snp[3].split('.')[1])
        self.strand = snp[5]
        self.pos = snp[7]
        self.id = snp[8]
        if self.strand == "-":
            self.aa = gen.reverse_complement(snp[9])
            self.ma = gen.reverse_complement(snp[10])
        else:
            self.aa = snp[9]
            self.ma = snp[10]
        self.rel_pos = int(snp[11])
        self.seq = cds_list[self.transcript_id][self.exon_id]

    def get_seq(self, cds_list):
        return cds_list[self.transcript_id][self.exon_id]

    def interval_seq(self):
        print("{0} ** {1} ** {2}".format(self.seq[:self.rel_pos], self.seq[self.rel_pos], self.seq[1+self.rel_pos:]))

    def get_hexamer_list(self):
        hexamer_list = []
        if self.rel_pos - 5 >= 0 and self.rel_pos + 1 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos-5:self.rel_pos+1])
        if self.rel_pos - 4 >= 0 and self.rel_pos + 2 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos-4:self.rel_pos+2])
        if self.rel_pos - 3 >= 0 and self.rel_pos + 3 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos-3:self.rel_pos+3])
        if self.rel_pos - 2 >= 0 and self.rel_pos + 4 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos-2:self.rel_pos+4])
        if self.rel_pos - 1 >= 0 and self.rel_pos + 5 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos-1:self.rel_pos+5])
        if self.rel_pos + 6 <= len(self.seq):
            hexamer_list.append(self.seq[self.rel_pos:self.rel_pos+6])
        return hexamer_list


def get_ptc_snps(ptc_file, snp_rel_pos_file, ptc_snps_file):
    print('Generating file containing relative exon position of PTC snps...')
    ptcs = gen.read_many_fields(ptc_file, "\t")
    snps = gen.read_many_fields(snp_rel_pos_file, "\t")
    snp_list = collections.defaultdict(lambda: collections.defaultdict())
    for snp in snps:
        snp_list[snp[8]] = snp
    with open(ptc_snps_file, "w") as outfile:
        for ptc in ptcs[1:]:
            ptc_id = ptc[8]
            ptc_snp = snp_list[ptc_id]
            outfile.write("{0}\n".format("\t".join(ptc_snp)))


def get_cds_interval_list(cds_intervals_fasta):
    cds_interval_names, cds_interval_seqs = gen.read_fasta(cds_intervals_fasta)
    cds_interval_list = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(cds_interval_names):
        transcript_id = name.split('.')[0]
        exon = int(name.split('.')[1])
        seq = cds_interval_seqs[i]
        if seq in stops:
            exon_id = 99999
        else:
            exon_id = exon
        cds_interval_list[transcript_id][exon_id] = seq
    return cds_interval_list


def ese_counter(ese_list):
    ese_overlap_counter = {}
    for ese in ese_list:
        ese_overlap_counter[ese] = 0
    return ese_overlap_counter

def ese_overlap(ptc_file, ese_list, cds_intervals_fasta):

    ptc_snps = gen.read_many_fields(ptc_file, "\t")
    cds_interval_list = get_cds_interval_list(cds_intervals_fasta)
    ese_overlap_count = 0
    ese_overlap_counter = ese_counter(ese_list)

    for ptc in ptc_snps[1:]:
        ptc = snp_info(ptc, cds_interval_list)
        hexamer_list = ptc.get_hexamer_list()
        ptc_ese_overlaps = [hexamer for hexamer in hexamer_list if hexamer in ese_list]
        if len(ptc_ese_overlaps) > 0:
            ese_overlap_count += 1
        for ese in ptc_ese_overlaps:
            ese_overlap_counter[ese] += 1
    return(ese_overlap_count, ese_overlap_counter)

def run_simulations(simulation_list, ptc_file, syn_nonsyn_file, ese_list, cds_intervals_fasta):

    overlap_counts = []
    overlap_counters = []

    simulation_count = len(simulation_list)

    for i, simulation in enumerate(simulation_list):
        print("{0}/{1} simulation running...".format(i+1, simulation_count))
        np.random.seed()
        random_choice = np.random.random(1)[0]
        pseudo_ptcs = "./temp_data/pseudo_ptcs.{0}.txt".format(random_choice)
        pseudo_other = "./temp_data/other.{0}.txt".format(random_choice)
        so.generate_pseudo_ptc_snps(ptc_file, syn_nonsyn_file, pseudo_ptcs, pseudo_other, group_by_gene=False, without_replacement=True, match_allele_frequency=True, match_allele_frequency_window=0.05)
        overlap_count, overlap_counter = ese_overlap(pseudo_ptcs, ese_list, cds_intervals_fasta)
        overlap_counts.append(overlap_count)
        overlap_counters.append(overlap_counter)
        gen.remove_file(pseudo_ptcs)
        gen.remove_file(pseudo_other)

    return(overlap_counts, overlap_counters)



def simulations(required_simulations, ptc_file, syn_nonsyn_file, ese_list, cds_intervals_fasta):
    simulation_list = list(range(1,1+required_simulations))
    processes = gen.run_in_parallel(simulation_list, ["foo", ptc_file, syn_nonsyn_file, ese_list, cds_intervals_fasta], run_simulations)
    overlap_count_list = []
    overlap_counter_list = []
    for process in processes:
        overlap_counts, overlap_counter = process.get()
        overlap_count_list.extend(overlap_counts)
        overlap_counter_list.extend(overlap_counter)
    return overlap_count_list, overlap_counter_list


def main():

    description = "Get the number of PTCs that overlap an ESE, and simulate with allele frequency matched SNPs."
    args = gen.parse_arguments(description, ["out_prefix", "ese_file", "simulations"])
    out_prefix, ese_file, required_simulations = args.out_prefix, args.ese_file, args.simulations

    ptc_file = "{0}_ptc_file.txt".format(out_prefix)
    other_snps_file = "{0}_syn_nonsyn_file.txt".format(out_prefix)
    snp_rel_pos_file = "{0}_SNP_relative_exon_position.bed".format(out_prefix)
    cds_intervals_fasta = "{0}_CDS_intervals.fasta".format(out_prefix)

    start = time.time()

    gen.create_directory("./temp_data")

    ptc_snps_file = "{0}_ptc_SNP_relative_exon_position.txt".format(out_prefix)
    if not os.path.exists(ptc_snps_file):
        get_ptc_snps(ptc_file, snp_rel_pos_file, ptc_snps_file)

    ese_list = [ese[0] for ese in gen.read_many_fields(ese_file, "\t") if "#" not in ese[0]]
    cds_interval_list = get_cds_interval_list(cds_intervals_fasta)

    real_overlap_count, real_ese_overlap_counter = ese_overlap(ptc_snps_file, ese_list, cds_intervals_fasta)
    simulation_overlap_counts, simulation_overlap_counters = simulations(int(required_simulations), ptc_file, other_snps_file, ese_list, cds_intervals_fasta)

    output_file = "{0}_ese_overlap_simulation.csv".format(out_prefix)
    with open(output_file, "w") as outfile:
        outfile.write('id,snps_overlapping_eses,{0}\n'.format(",".join(sorted(ese_list))))
        outfile.write('real,{0},{1}\n'.format(real_overlap_count, ",".join([str(real_ese_overlap_counter[ese]) for ese in sorted(real_ese_overlap_counter)])))

        for i, count in enumerate(simulation_overlap_counts):
            counter = simulation_overlap_counters[i]
            outfile.write('{0},{1},{2}\n'.format(i, count, ",".join([str(counter[ese]) for ese in sorted(counter)])))

    gen.get_time(start)

if __name__ == "__main__":
    main()
