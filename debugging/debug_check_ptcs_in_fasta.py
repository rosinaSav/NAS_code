'''
Check whether SNPs that are being labelled PTC snps are actually generating a PTC
'''


import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re

ptc_file = "./results/clean_run/clean_run_ptc_file.txt"
snp_file = "./results/clean_run/clean_run_SNP_file.txt"
snp_relative_exon_position_file = "./results/clean_run/clean_run_SNP_relative_exon_position.bed"
snp_cds_file = "./results/clean_run/clean_run_CDS_SNP_intersect.bed"
coding_exons_file = "./results/clean_run/clean_run_coding_exons.bed"
coding_exons_fasta = "./results/clean_run/clean_run_coding_exons.fasta"
coding_exons_file_changed = "./results/clean_run/clean_run_coding_exons_formatted.bed"


genome_file = "./source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
out = "./temp_data/ptc_fasta.fasta"
ptc_bed = "./temp_data/ptc_bed.bed"


cds_fasta = "./results/clean_run/clean_run_CDS.fasta"

ptcs = gen.read_many_fields(ptc_file, "\t")
names, cdss = gen.read_fasta(cds_fasta)

cds_list = {}
for i, name in enumerate(names):
    cds_list[name] = cdss[i]

stops = ["TAA", "TAG", "TGA"]

without_ptcs = []

print(len(ptcs))

for ptc in ptcs[1:]:

    ptc_count = 0

    exon_id = ptc[3]
    cds_id = exon_id.split('.')[0]

    snp_aa = ptc[9]
    snp_ma = ptc[10]
    rel_pos = int(ptc[11])

    strand = ptc[5]

    cds_seq = cds_list[cds_id]
    #
    # if strand == "-":
    #     snp_ma = gen.reverse_complement(snp_ma)

    mutation_seq = cds_seq[:rel_pos] + snp_ma + cds_seq[1+rel_pos:]

    codons = re.findall('.{3}', mutation_seq[:-3])
    for stop in stops:
        if stop in codons:
            ptc_count += 1

    if ptc_count == 0:
        without_ptcs.append(ptc[8])
        print(snp_aa, snp_ma, strand)
        print(cds_seq[:rel_pos], cds_seq[rel_pos], cds_seq[1+rel_pos:])
        print(cds_seq[:rel_pos], snp_ma, cds_seq[1+rel_pos:])

with open('debugging/error_ids.txt', "w") as out:
    for id in without_ptcs:
        out.write("{0}\n".format(id))

print(without_ptcs)
print(len(without_ptcs))
