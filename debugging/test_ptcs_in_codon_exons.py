import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re

codon_exon_fasta = "./results/clean_run/clean_run_coding_exons.fasta"
ptc_file = "./results/clean_run/clean_run_ptc_file.txt"
snp_file = "./results/clean_run/clean_run_SNP_relative_exon_position.bed"
cds_fasta = "./results/clean_run/clean_run_CDS.fasta"

exon_names, exon_seqs = gen.read_fasta(codon_exon_fasta)
cds_names, cds_seqs = gen.read_fasta(cds_fasta)
ptcs = gen.read_many_fields(ptc_file, "\t")
snps = gen.read_many_fields(snp_file, "\t")


class info(object):

    def __init__(self, ptc):
        self.exon_id = ptc[3]
        self.cds_id = ptc[3].split('.')[0]
        self.strand = ptc[5]
        self.snp_pos = ptc[7]
        self.snp_id = ptc[8]
        self.snp_aa = ptc[9]
        self.snp_ma = ptc[10]
        self.snp_rel_pos = int(ptc[11])
        # self.snp_cds_codon = ptc[13].split('$')[0].split('=')[-1]
        # self.snp_mut_codon = ptc[13].split('$')[1].split('=')[-1]

# put the exons into a dictionary sorted by name
exon_list = {}
for i, exon in enumerate(exon_names):
    exon_list[exon] = exon_seqs[i]

# get a dict of all snps by snp id
snp_list = {}
for snp in snps:
    snp_list[info(snp).snp_id] = snp

# get all the cds
cds_list = {}
for i, name in enumerate(cds_names):
    cds_list[name] = cds_seqs[i]

# get the snp info for the ptc
ptc_list = []
for ptc in ptcs[1:]:
    ptc = info(ptc)
    if ptc.snp_id in snp_list:
        ptc_list.append(snp_list[ptc.snp_id])

ptc_cds_positions = {}
for ptc in ptcs[1:]:
    ptc_cds_positions[info(ptc).snp_id] = ptc

def print_entry(ptc, ptc_entry, exon_seq, cds_seq):
    print('SNP ID: ', ptc.snp_id)
    print('Exon ID: ', ptc.exon_id)
    print('CDS ID: ', ptc.cds_id)
    print('\n')
    print('Exon seq')
    print(exon_seq)
    print('\n')
    print('SNP relative to exon index: ', ptc.snp_rel_pos)
    print('AA: ', ptc.snp_aa)
    print('MA: ', ptc.snp_ma)
    print(real_exon_parts)
    print(mut_exon_parts)
    print('\n')
    print('CDS seq')
    print(cds_seq[exon_index:], '**', cds_seq[exon_index:exon_index + len(exon_seq) + 1], '**', cds_seq[exon_index + len(exon_seq) + 1:])
    print('\n')
    print('SNP relative to CDS index: ', ptc_entry.snp_rel_pos)
    print('AA: ', ptc_entry.snp_aa)
    print('MA: ', ptc_entry.snp_ma)
    print(cds_seq[:ptc_entry.snp_rel_pos], '**', cds_seq[ptc_entry.snp_rel_pos + 1], '**', cds_seq[ptc_entry.snp_rel_pos + 1:])
    print(cds_seq[:ptc_entry.snp_rel_pos], '**', ptc_entry.snp_ma, '**', cds_seq[ptc_entry.snp_rel_pos + 1:])


aa_fail = 0
fail_list = []


# for each of the PTCs, get the seq and check ptc
for i, ptc in enumerate(ptc_list[1:10]):
    ptc = info(ptc)
    exon_seq = exon_list[ptc.exon_id]

    if ptc.strand == "-":
        exon_seq = list(exon_seq)
        exon_seq = [gen.reverse_complement(nt) for nt in exon_seq]
        exon_seq = "".join(exon_seq)
        exon_seq = exon_seq[::-1]
        ptc.snp_aa = gen.reverse_complement(ptc.snp_aa)
        ptc.snp_ma = gen.reverse_complement(ptc.snp_ma)

    else:
        if ptc.snp_aa != exon_seq[ptc.snp_rel_pos]:
            aa_fail += 1

        mut_exon = "{0}{1}{2}".format(exon_seq[:ptc.snp_rel_pos], ptc.snp_ma, exon_seq[1+ptc.snp_rel_pos:])
        real_exon_parts = "{0} ** {1} ** {2}".format(exon_seq[:ptc.snp_rel_pos], exon_seq[ptc.snp_rel_pos], exon_seq[1+ptc.snp_rel_pos:])
        mut_exon_parts = "{0} ** {1} ** {2}".format(exon_seq[:ptc.snp_rel_pos], ptc.snp_ma, exon_seq[1+ptc.snp_rel_pos:])

        codon1 = mut_exon[ptc.snp_rel_pos:ptc.snp_rel_pos+3]
        codon2 = mut_exon[ptc.snp_rel_pos-1:ptc.snp_rel_pos+2]
        codon3 = mut_exon[ptc.snp_rel_pos-2:ptc.snp_rel_pos+1]

        stops = ["TAA", "TAG", "TGA"]
        if codon1 not in stops and codon2 not in stops and codon3 not in stops:
            fail_list.append(ptc.snp_id)


        cds_seq = cds_list[ptc.cds_id]
        exon_index = cds_list[ptc.cds_id].index(exon_seq)
        cds_pos = exon_index + ptc.snp_rel_pos

        ptc_entry = ptc_cds_positions[ptc.snp_id]
        ptc_entry = info(ptc_entry)

        if ptc.snp_id == "rs142458192":
            print_entry(ptc, ptc_entry, exon_seq, cds_seq)

print(fail_list)

print(aa_fail)
print(len(fail_list))
