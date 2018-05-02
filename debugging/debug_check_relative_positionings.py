'''
Check whether the relative positions are being calculated correctly
'''

import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re

cds_bed_file = "./results/clean_run/clean_run_CDS.bed"
cds_fasta_file = "./results/clean_run/clean_run_CDS_intervals.fasta"
snp_cds_relative_positions_file = "./results/clean_run/clean_run_SNP_relative_exon_position.bed"

cds_fasta_names, cds_fasta_seqs = gen.read_fasta(cds_fasta_file)
cds_bed = gen.read_many_fields(cds_bed_file, "\t")
snp_cds_relative_positions = gen.read_many_fields(snp_cds_relative_positions_file, "\t")

class Snp(object):
    def __init__(self, ptc):
        self.exon_id = ptc[3]
        self.cds_id = ptc[3].split('.')[0]
        self.strand = ptc[5]
        self.snp_pos = ptc[7]
        self.snp_id = ptc[8]
        self.snp_aa = ptc[9]
        self.snp_ma = ptc[10]
        self.snp_rel_pos = int(ptc[11])

cds_fasta_list = {}
for i, name in enumerate(cds_fasta_names):
    name = name.split('.')[0] + '.' + name.split('.')[1]
    cds_fasta_list[name] = cds_fasta_seqs[i]

# get a dictionary of snps listed by id
snp_list = {}
for snp in snp_cds_relative_positions:
    snp_list[Snp(snp).snp_id] = snp

# set the query snp
test_snp_id = "rs12751100"
query_snp = snp_list[test_snp_id]

# get the parts of the snp
query_snp = Snp(query_snp)
print(query_snp.cds_id)

cds_bed_list = {}
for cds in cds_bed:
    id_split = cds[3].split('.')
    cds_id = id_split[0]
    if cds_id == query_snp.cds_id:
        if cds[4] == "stop_codon":
            exon_id = 99999
        else:
            exon_id = int(id_split[1])
        cds_bed_list[exon_id] = cds

exon_id = int(query_snp.exon_id.split('.')[-1])
cds_rel_pos = 0

for i in range(1, exon_id):
    part = cds_bed_list[i]
    cds_exon_length = int(part[2]) - int(part[1])
    cds_rel_pos += cds_exon_length
    id_split = part[3].split('.')
    cds_exon_id = "{0}.{1}".format(id_split[0], id_split[1])
    print(cds_exon_id)
    print(cds_fasta_list[cds_exon_id])
    print(len(cds_fasta_list[cds_exon_id]))

cds_rel_pos += query_snp.snp_rel_pos

# print(cds_rel_pos)

# for cds_exon in sorted(cds_bed_list):
#     print(cds_exon)
