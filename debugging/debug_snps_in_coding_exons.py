'''
Check whether the ancestral alleles of all SNPS match.
Check how many SNPs could potentially be called a PTC.
'''


import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen

ptc_file = "./results/clean_run/clean_run_ptc_file.txt"
snp_relative_exon_position_file = "./results/clean_run/clean_run_SNP_relative_exon_position.bed"
snp_cds_file = "./results/clean_run/clean_run_CDS_SNP_intersect.bed"
coding_exons_file = "./results/clean_run/clean_run_coding_exons.bed"
coding_exons_fasta = "./results/clean_run/clean_run_coding_exons.fasta"
coding_exons_file_changed = "./results/clean_run/clean_run_coding_exons_formatted.bed"
genome_file = "./source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
out = "./temp_data/ptc_fasta.fasta"
ptc_bed = "./temp_data/ptc_bed.bed"

# gen.extract_head_of_file(snp_relative_exon_position_file, 8000)

comps = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

names, exons = gen.read_fasta(coding_exons_fasta)
exon_list = {}
for i, name in enumerate(names):
    exon_list[name] = exons[i]

snps = gen.read_many_fields(snp_relative_exon_position_file, "\t")

not_count = 0
aa_fail = 0

# print(snps[0])

stops = ["TAA", "TAG", "TGA"]
poss_ptcs = 0

for i, snp in enumerate(snps[:]):
    exon_id = snp[3]
    snp_id = snp[8]
    aa = snp[9]
    ma = snp[10]
    rel_pos = int(snp[11])
    strand = snp[5]

    if len(aa) == 1 and exon_id in exon_list:

        exon_seq = exon_list[exon_id]

        if strand == "-":
            nts = list(exon_seq)
            seq_nts = [comps[nt] for nt in nts]
            exon_seq = "".join(seq_nts)

        seq_nt = exon_seq[rel_pos]

        if aa != seq_nt:
            aa_fail += 1

        new_seq = exon_seq[:rel_pos] + ma + exon_seq[1+rel_pos:]
        #
        # print(exon_seq[:rel_pos], exon_seq[rel_pos], exon_seq[1+rel_pos:])
        # print(exon_seq[:rel_pos], ma, exon_seq[1+rel_pos:])

        codon1 = new_seq[rel_pos-2:rel_pos+1]
        codon2 = new_seq[rel_pos-1:rel_pos+2]
        codon3 = new_seq[rel_pos:rel_pos+3]

        if codon1 in stops or codon2 in stops or codon3 in stops:
            poss_ptcs += 1

print(aa_fail)
print(poss_ptcs)
