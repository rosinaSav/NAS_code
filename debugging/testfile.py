import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen

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

ptcs = gen.read_many_fields(ptc_file, "\t")

comps = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}

def complement(seq):

    list_seq = list(seq)
    comp = [comps[i] for i in list_seq]
    return("".join(comp))


with open(ptc_bed, "w") as file:

    for ptc in ptcs[1:]:
        outs = []
        outs.append(ptc[6])
        outs.append(str(int(ptc[7])-3))
        outs.append(str(int(ptc[7])+2))
        outs.append(ptc[8])
        outs.append('.')
        outs.append(ptc[5])
        outs.append(ptc[9])
        outs.append(ptc[10])
        file.write("{0}\n".format("\t".join(outs)))
bo.fasta_from_intervals(ptc_bed, out, genome_file, names=True)

nptc = gen.read_many_fields(ptc_bed, "\t")
gnames, gseqs = gen.read_fasta(out)

g = {}
for name in gnames:
    g[name] = gseqs[gnames.index(name)]

count = 0
for ptc in nptc:
    id = ptc[3]
    seq = g[id]
    aa = ptc[6]
    ma = ptc[7]
    strand = ptc[5]

    if strand == "-":
        seq = complement(seq)
        ma = complement(ma)

    mut_seq = list(seq)
    mut_seq[2] = ma
    mut_seq = "".join(mut_seq)

    if strand == "-":
        mut_seq = mut_seq[::-1]

    if "TAA" not in mut_seq and "TAG" not in mut_seq and "TGA" not in mut_seq:
        count += 1

print(count)
