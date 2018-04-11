import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re
import collections

cds_intervals_file = "./results/clean_run/clean_run_CDS_intervals.fasta"
cds_file = "./results/clean_run/clean_run_CDS.fasta"
cds_bed_file = "./results/clean_run/clean_run_CDS.bed"

cds_interval_locations = gen.read_many_fields(cds_bed_file, "\t")
cds_names, cds_seqs = gen.read_fasta(cds_file)
cds_interval_names, cds_interval_seqs = gen.read_fasta(cds_intervals_file)

# get all full coding sequences
cds_list = {}
for i, name in enumerate(cds_names):
    cds_list[name] = cds_seqs[i]

# get the strands
strands = {}
for cds in cds_interval_locations:
    cds_id = cds[3].split('.')[0]
    if cds_id not in strands:
        strands[cds_id] = cds[-1]

# tag the stop codons
for i, name in enumerate(cds_interval_names):
    if len(cds_interval_seqs[i]) == 3:
        cds_interval_names[i] = name + '.stop_codon'


# id_list = ["ENST00000381567", "ENST00000381130"]
id_list = ["ENST00000325425"]
interval_list = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(cds_interval_names):
    splits = name.split('.')
    cds_id = splits[0]

    if cds_id in id_list:
        if splits[-1] == "stop_codon":
            exon_id = 99999
        else:
            exon_id = int(splits[1])

        interval_list[cds_id][exon_id] = cds_interval_seqs[i]

cds_seqs = {}

for cds_id in id_list:
    strand = strands[cds_id]
    exons = sorted(interval_list[cds_id])
    cds_seq = ""
    for exon in exons:
        print("{0}.{1}".format(cds_id, exon))
        print(interval_list[cds_id][exon])
        cds_seq += interval_list[cds_id][exon]
    cds_seqs[cds_id] = cds_seq

for cds_id in cds_seqs:
    print(cds_id, len(cds_seqs[cds_id]))
    # print(cds_seqs[cds_id])
