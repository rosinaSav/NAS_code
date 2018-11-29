import generic as gen
import disease_snps_ops as dso
import bed_ops as bo
import collections
import re
import os
import numpy as np
import pandas as pd
import random

ptcs = "results/clean_run_2/clean_run_ptc_SNP_relative_exon_position.txt"
cds_seqs_fasta = "results/clean_run_2/clean_run_CDS.fasta"
coding_exons_fasta = "results/clean_run_2/clean_run_CDS_intervals.fasta"
exons_bed = "results/clean_run_2/clean_run_exons.bed"
genome_fasta = "../source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa"

queries = [
    "ENST00000255409.8",
    "ENST00000456763.12",
    "ENST00000265316.3",
    "ENST00000272065.5",
]

changes = {
    "A": "T", "T": "A", "G": "C", "C": "G"
}

cds_names, cds_seqs = gen.read_fasta(cds_seqs_fasta)
exon_names, exon_seqs = gen.read_fasta(coding_exons_fasta)
ptcs = gen.read_many_fields(ptcs, "\t")
all_exons_list = gen.read_many_fields(exons_bed, "\t")

cds = {name: cds_seqs[i] for i, name in enumerate(cds_names)}

exons = collections.defaultdict(lambda: collections.defaultdict())
stops = collections.defaultdict(lambda: collections.defaultdict())
for i, name in enumerate(exon_names):
    id = name.split(".")[0]
    exon_no = int(name.split(".")[1].split("(")[0])
    if len(exon_seqs[i]) > 3:
        exons[id][exon_no] = exon_seqs[i]
    else:
        stops[id][exon_no] = exon_seqs[i]

ptcs = [i for i in ptcs if i[3] in queries]
ptc_list = {i[3]: i for i in ptcs}

all_exons = collections.defaultdict(lambda: collections.defaultdict())
for i, exon in enumerate(all_exons_list):
    id = exon[3].split(".")[0]
    exon_no = int(exon[3].split(".")[1])
    all_exons[id][exon_no] = exon

output_file = "results/experimental/candidate_seq_info.csv"

with open(output_file, "w") as outfile:
    outfile.write("transcript_id,exon_number,mutation,mutation_position_in_exon,codon_change,upstream_start_index,,5'_exon_seq,exon_seq,3'_exon_seq,3'_is_final_exon,,5'_intron_seq,3'_intron_seq,,full_cds_seq\n")
    for q in queries:
        id = q.split(".")[0]
        exon_no = int(q.split(".")[1])
        exon = exons[id][exon_no]
        ptc = ptc_list[q]

        start = int(ptc[1])
        stop = int(ptc[2])
        strand = ptc[5]
        aa = ptc[9]
        ma = ptc[10]
        pos = int(ptc[11])

        if strand == "-":
            aa = changes[aa]
            ma = changes[ma]

        upstream_exon = exons[id][exon_no-1]
        downstream_exon = exons[id][exon_no+1]

        if id in stops:
            if exon_no+1 in stops[id]:
                downstream_exon += stops[id][exon_no+1]

        cds_seq = cds[id]

        upstream_info = all_exons[id][exon_no-1]
        downstream_info = all_exons[id][exon_no+1]

        if strand == "-":
            up_start = stop
            up_stop = int(upstream_info[1])
            upstream_intron_info = [ptc[0], up_start, up_stop, "{0}.{1}-{2}".format(id, exon_no-1, exon_no), ".", strand]
            down_start = int(downstream_info[2])
            down_stop = start
            downstream_intron_info = [ptc[0], down_start, down_stop, "{0}.{1}-{2}".format(id, exon_no, exon_no+1), ".", strand]
        else:
            up_start = int(upstream_info[2])
            up_stop = start
            upstream_intron_info = [ptc[0], up_start, up_stop, "{0}.{1}-{2}".format(id, exon_no-1, exon_no), ".", strand]
            down_start = stop
            down_stop = int(downstream_info[1])
            downstream_intron_info = [ptc[0], down_start, down_stop, "{0}.{1}-{2}".format(id, exon_no, exon_no+1), ".", strand]

        ri = random.random()
        temp_file = "temp_files/{0}.bed".format(ri)
        with open(temp_file, "w") as temp_output:
            temp_output.write("{0}\n".format("\t".join(gen.stringify(upstream_intron_info))))
            temp_output.write("{0}\n".format("\t".join(gen.stringify(downstream_intron_info))))

        temp_seqs = "temp_files/{0}.fasta".format(ri)
        bo.fasta_from_intervals(temp_file, temp_seqs, genome_fasta, names = True)

        intron_seqs = gen.read_fasta(temp_seqs)[1]
        upstream_intron_seq = intron_seqs[0].lower()
        downstream_intron_seq = intron_seqs[1].lower()

        gen.remove_file(temp_file)
        gen.remove_file(temp_seqs)


        exon_start_index = cds_seq.index(exon) % 3
        if exon_start_index == 1:
            start_pos = 2
        elif exon_start_index == 2:
            start_pos = 1
        else:
            start_pos = 0

        cds_before = cds_seq[:start_pos]
        cds_exon = cds_seq[start_pos:start_pos+len(exon)]
        cds_after = cds_seq[start_pos+len(exon):]

        codons = []
        for i in range(0, len(cds_before), 3):
            codons.append(cds_before[i:i+3])

        query_codon = ""
        exon_codons = []
        for i in range(start_pos, len(exon), 3):
            codon = exon[i:i+3]
            exon_codons.append(codon)
            if pos in list(range(i, i+3)):
                query_codon = codon

        mut_exon = list(exon)
        mut_exon[pos] = ma
        mut_exon = "".join(mut_exon)

        mut_codon = ""
        mut_exon_codons = []
        for i in range(start_pos, len(exon), 3):
            codon = mut_exon[i:i+3]
            mut_exon_codons.append(codon)
            if pos in list(range(i, i+3)):
                mut_codon = codon

        upstream_start_index = cds_seq.index(upstream_exon) % 3

        final_exon = ""
        if exon_no + 1 == len(exons[id]):
            final_exon = "*"



        outfile.write("{0},{1},{2} > {3},{4},{5} > {6},{7},,{8},{9},{10},{11},,{12},{13},,{14}\n".format(id, exon_no, aa, ma, pos+1, query_codon, mut_codon, upstream_start_index, upstream_exon, exon, downstream_exon, final_exon, upstream_intron_seq, downstream_intron_seq, cds_seq))
