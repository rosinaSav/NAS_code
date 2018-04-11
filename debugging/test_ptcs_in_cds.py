import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re

cds_fasta = "./results/clean_run/clean_run_CDS.fasta"
ptc_file = "./results/clean_run/clean_run_ptc_file.txt"

cds_names, cds_seqs = gen.read_fasta(cds_fasta)
ptcs = gen.read_many_fields(ptc_file, "\t")

# put the cdss into a dictionary sorted by name
cds_list = {}
for i, cds in enumerate(cds_names):
    cds_list[cds] = cds_seqs[i]

class PTC(object):

    def __init__(self, ptc):
        self.exon_id = ptc[3]
        self.cds_id = ptc[3].split('.')[0]
        self.strand = ptc[5]
        self.snp_pos = ptc[7]
        self.snp_id = ptc[8]
        self.snp_aa = ptc[9]
        self.snp_ma = ptc[10]
        self.snp_rel_pos = int(ptc[11])
        self.snp_cds_codon = ptc[13].split('$')[0].split('=')[-1]
        self.snp_mut_codon = ptc[13].split('$')[1].split('=')[-1]

aa_fail = 0
cds_codon_fail = 0
mut_codon_fail = 0
ptc_fail = 0

# for each of the PTCs, get the seq and check ptc
for i, ptc in enumerate(ptcs[1:]):
    ptc = PTC(ptc)
    cds_seq = cds_list[ptc.cds_id]

    # revserse complement the ancestral allele and snp allele if on the - strand
    if ptc.strand == "-":
        ptc.snp_aa = gen.reverse_complement(ptc.snp_aa)
        ptc.snp_ma = gen.reverse_complement(ptc.snp_ma)

    mut_cds = "{0}{1}{2}".format(cds_seq[:ptc.snp_rel_pos], ptc.snp_ma, cds_seq[1+ptc.snp_rel_pos:])

    real_cds_parts = "{0} {1} {2}".format(cds_seq[:ptc.snp_rel_pos], cds_seq[ptc.snp_rel_pos], cds_seq[1+ptc.snp_rel_pos:])
    mut_cds_parts = "{0} {1} {2}".format(cds_seq[:ptc.snp_rel_pos], ptc.snp_ma, cds_seq[1+ptc.snp_rel_pos:])


    # check that the aa match
    if ptc.snp_aa != cds_seq[ptc.snp_rel_pos]:
        aa_fail += 1

    snp_index = len(cds_seq[:ptc.snp_rel_pos]) % 3
    if snp_index == 0:
        cds_codon = cds_seq[ptc.snp_rel_pos:ptc.snp_rel_pos+3]
        mut_codon = mut_cds[ptc.snp_rel_pos:ptc.snp_rel_pos+3]
    elif snp_index == 1:
        cds_codon = cds_seq[ptc.snp_rel_pos-1:ptc.snp_rel_pos+2]
        mut_codon = mut_cds[ptc.snp_rel_pos-1:ptc.snp_rel_pos+2]
    elif snp_index == 2:
        cds_codon = cds_seq[ptc.snp_rel_pos-2:ptc.snp_rel_pos+1]
        mut_codon = mut_cds[ptc.snp_rel_pos-2:ptc.snp_rel_pos+1]

    # check that the called cds codon matches
    if ptc.snp_cds_codon != cds_codon:
        cds_codon_fail += 1

    # check that the called mut codon matches
    if ptc.snp_mut_codon != mut_codon:
        mut_codon_fail += 1

    # get all the codons for the mutated cds
    mut_codons = re.findall('.{3}', mut_cds[:-3])
    if "TAA" not in mut_codons and "TAG" not in mut_codons and "TGA" not in mut_codons:
        ptc_fail += 1

    if ptc.snp_id == "rs12751100":
        print(real_cds_parts)


print(aa_fail)
print(cds_codon_fail)
print(mut_codon_fail)
print(ptc_fail)
