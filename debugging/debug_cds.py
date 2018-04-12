import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re
import collections
import numpy as np

cds_intervals_file = "./results/clean_run/clean_run_CDS_intervals.fasta"
cds_file = "./results/clean_run/clean_run_CDS.fasta"
cds_bed_file = "./results/clean_run/clean_run_CDS.bed"

interval_names, interval_seqs = gen.read_fasta(cds_intervals_file)
cds_names, cds_fasta = gen.read_fasta(cds_file)

cds_list = {}
for i, name in enumerate(cds_names):
    cds_list[name] = cds_fasta[i]


query_name = "ENST00000325425"

required_names = []
required_seqs = []

for i, name in enumerate(interval_names):
    if query_name in name:
        required_names.append(name)
        required_seqs.append(interval_seqs[i])

exon_list = {}
for i, name in enumerate(required_names):
    splits = name.split(".")
    if len(required_seqs[i]) == 3:
        exon = 99999
    else:
        exon = int(splits[1])
    exon_list[exon] = required_seqs[i]

seq = ""
for i in sorted(exon_list):
    print(i)
    seq += exon_list[i]


print('\nGenerated seq:')
print(seq)
print('\nCDS fasta seq:')
print(cds_list[query_name])

test_fasta = "./temp_data/test_outfile.fasta"
test_names, test_seqs = gen.read_fasta(test_fasta)
print('\ntest fasta seq:')
print(test_seqs[0])

if seq == cds_list[query_name] and seq == test_seqs[0] and cds_list[query_name] == test_seqs[0]:
    print("same")

#
# cds_bed_file = "./results/clean_run/clean_run_CDS.bed"
# cds_intervals_fasta_file = "./results/clean_run/clean_run_CDS_intervals_not_filtered.bed"
# outfile = "temp_data/test_outfile.fasta"
# genome_fasta = "./source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
#
# bo.extract_cds_from_bed(cds_bed_file, outfile, genome_fasta, cds_intervals_fasta_file)

rand = ""
nts = ["A", "C", "G", "T"]
for i in range(0, 100):
    rand += np.random.choice(nts, 1, replace=True)[0]
print(rand)
