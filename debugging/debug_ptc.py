'''
Try and work out why 1 ptc is not correctly called
'''


import bed_ops as bo
import generic as gen
import bam_ops as ba
import numpy as np
import generic as gen
import re
import collections


ptc_file = "./results/clean_run_2/clean_run_ptc_file.txt"
snp_file = "./results/clean_run_2/clean_run_SNP_relative_exon_position.bed"
coding_exons_file = "./results/clean_run_2/clean_run_coding_exons.fasta"
cds_intervals_file = "./results/clean_run_2/clean_run_CDS_intervals.fasta"
cds_fasta_file = "./results/clean_run_2/clean_run_CDS.fasta"
snps_file = "./results/clean_run_2/clean_run_SNP_file.txt"

test_snp_ids = ["rs58993422", "rs2161752"]

def get_sequences():
    ptcs = gen.read_many_fields(ptc_file, "\t")
    snps = gen.read_many_fields(snp_file, "\t")
    exon_names, exon_seqs = gen.read_fasta(coding_exons_file)
    cds_interval_names, cds_interval_seqs = gen.read_fasta(cds_intervals_file)
    cds_names, cds_seqs = gen.read_fasta(cds_fasta_file)


    test_seqs = ["ENST00000400186"]
    test_exons = [24, 27]


    tests = {
        24: ["G", "A", 95, 3572],
        27: ["C", "T", 145, 4098],
    }

    stops = ["TAA", "TAG", "TGA"]

    cds_list = {}
    for i, name in enumerate(cds_names):
        cds_list[name] = cds_seqs[i]

    intervals = collections.defaultdict(lambda: collections.defaultdict())
    for i, name in enumerate(cds_interval_names):
        name = name.split('.')
        if name[0] in test_seqs:
            if cds_interval_seqs[i] in stops:
                exon_id = 99999
            else:
                exon_id = int(name[1])
            intervals[name[0]][int(exon_id)] = cds_interval_seqs[i]

    exon_bits = collections.defaultdict(lambda: collections.defaultdict())

    for gene in intervals:
        exons = sorted(intervals[gene])
        split_seq = ''
        seq = ''
        for exon in exons:
            if exon in test_exons:
                split_seq += " **{0}** ".format(exon)
            else:
                split_seq += " *{0}* ".format(exon)

            seq += intervals[gene][exon]
            split_seq += intervals[gene][exon]
            exon_bits[gene][int(exon)] = intervals[gene][exon]

            if exon in test_exons:
                split_seq += " **{0}** ".format(exon)
            else:
                split_seq += " *{0}* ".format(exon)

        # print(split_seq)
        #
        # print(seq[:4098], "**", seq[4098], "**", seq[4099:4104])
        # print(seq[:3572], "**", seq[3572], "**", seq[3572:3578])


        # list_seq = list(seq)
        # for i, nt in enumerate(list_seq):
        #     if nt != test_seq[i]:
        #         print(i, nt)


    print("\n")

    ptc_list = {}
    for ptc in ptcs[1:]:
        ptc_list[ptc[8]] = ptc


    ptc_ids = ["rs58993422", "rs2161752"]

    snp_list = {}
    for snp in snps:
        snp_list[snp[8]] = snp

    exon_list = {}
    for i, name in enumerate(exon_names):
        exon_list[name.split("(")[0]] = exon_seqs[i]


    aa_fail = []
    ptc_fail = []


    for ptc in ptc_ids:
        snp = snp_list[ptc]
        ptc_entry = ptc_list[ptc]

        exon_id = snp[3]
        strand = snp[5]
        aa = snp[9]
        ma = snp[10]
        rel_pos = int(snp[11])

        seq = exon_list[exon_id]

        transcript = exon_id.split('.')[0]
        exon = exon_id.split('.')[1]

        # print(transcript, exon)
        bit = exon_bits[transcript][int(exon)]

        if strand == "+":
            if aa != seq[rel_pos]:
                aa_fail.append(1)
            else:
                aa_fail.append(0)
                mut_seq = seq[:rel_pos] + ma + seq[rel_pos + 1:]
                # print(seq[:rel_pos], "**", ma, seq[rel_pos + 1:])
                codon1 = mut_seq[rel_pos-2:rel_pos+1]
                codon2 = mut_seq[rel_pos-1:rel_pos+2]
                codon3 = mut_seq[rel_pos:rel_pos+3]

            if codon1 not in stops and codon2 not in stops and codon3 not in stops:
                ptc_fail.append(1)
                print(ptc, exon_id, strand, aa, ma, rel_pos, ptc_entry[11], len(mut_seq), len(bit))
                # print(seq[:rel_pos], "**", ma, "**", seq[rel_pos + 1:])
            else:
                ptc_fail.append(0)

        else:
            # seq = list(seq)
            # seq = [gen.reverse_complement(nt) for nt in seq]
            # seq = "".join(seq)
            aa = gen.reverse_complement(aa)
            ma = gen.reverse_complement(ma)
            if aa != seq[rel_pos]:
                aa_fail.append(1)
                # ma = gen.reverse_complement(ma)


            else:
                aa_fail.append(0)
                # ma = gen.reverse_complement(ma)
                mut_seq = seq[:rel_pos] + ma + seq[rel_pos + 1:]
                # print(seq[:rel_pos], "**", ma, seq[rel_pos + 1:])
                codon1 = mut_seq[rel_pos-2:rel_pos+1]
                codon2 = mut_seq[rel_pos-1:rel_pos+2]
                codon3 = mut_seq[rel_pos:rel_pos+3]

                if codon1 not in stops and codon2 not in stops and codon3 not in stops:
                    ptc_fail.append(1)
                    print(ptc, exon_id, strand, aa, ma, rel_pos, ptc_entry[11], len(mut_seq), len(bit))
                    print(seq[:rel_pos], "**", ma, "**", seq[rel_pos + 1:])
                else:
                    ptc_fail.append(0)

    print(sum(aa_fail), sum(aa_fail)/len(aa_fail))
    print(sum(ptc_fail), sum(ptc_fail)/len(ptc_fail))

    index = 0
    total = 0

    full_seq = ''

    for gene in exon_bits:
        for bit in sorted(exon_bits[gene]):
            exon_bit = exon_bits[gene][bit]
            # print("{0}.{1} {2}".format(gene, bit, len(exon_bit)))
            cds_bit = cds_list[gene][index:index+len(exon_bits[gene][bit])]
            full_seq += cds_bit

            total += len(exon_bit)

            if bit in tests:
                seq_index = tests[bit][2]
                print(seq_index)
                print(len(exon_bit))
                print(index + tests[bit][2], tests[bit][3])
                exon_mut = "{0} ** {1} ** {2}".format(exon_bit[:seq_index], exon_bit[seq_index], exon_bit[seq_index+1:])
                cds_mut = "{0} ** {1} ** {2}".format(cds_bit[:seq_index], cds_bit[seq_index], cds_bit[seq_index+1:])
                aa = tests[bit][0]
                ma = tests[bit][1]

                if index < tests[bit][3] < total:
                    query_index = tests[bit][3]
                    query_seq = "{0} ** {1} ** {2}".format(full_seq[:query_index], full_seq[query_index], full_seq[query_index+1:])
                else:
                    query_seq = ''

            # else:
            #     exon_mut = exon_bit
            #     cds_mut = cds_bit
            #     aa = ''
            #     ma = ''
            #     query_seq = ''


                print("{0}.{1} {2} {3}".format(gene, bit, aa, ma))
                print(exon_mut, len(exon_bit), total)
                print(cds_mut, len(cds_bit), total)
                print(query_seq)
            # print(exon_bit == cds_bit)
            index += len(exon_bit)


def get_snps():
    snps = gen.read_many_fields(snps_file, "\t")
    snp_list = {}
    for snp in snps[1:]:
        if snp[8] in test_snp_ids:
            snp_list[snp[8]] = snp


    for snp in snp_list:
        single_snp = snp_list[snp]
        print(single_snp[:15])



if __name__ == "__main__":
    # get_sequences()
    get_snps()
