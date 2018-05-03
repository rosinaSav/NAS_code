import generic as gen
import bed_ops as bo
import collections
import re
import os

description = "Get PTC information"
args = gen.parse_arguments(description, ["fresh_run", "create_files"], flags = [0,1])
fresh_run, create_files =  args.fresh_run, args.create_files

ptc_file = './results/clean_run_2/clean_run_ptc_file.txt'
gtf_file = '../source_data/Homo_sapiens.GRCh37.87.gtf'
exon_fasta_file = './results/clean_run_2/clean_run_CDS_intervals.fasta';
genes_file = './temp_data/gene_list.bed'
cds_file = './results/clean_run_2/clean_run_CDS.bed'
cds_fasta_file = './results/clean_run_2/clean_run_CDS.fasta'
exons_file = './results/clean_run_2/clean_run_exons.bed'
introns_file = './temp_data/intron_list.bed'
transcript_file = "./temp_data/transcript_list.bed"
cds_rel_pos_file = './results/clean_run_2/clean_run_ptc_SNP_relative_exon_position.txt'

ptc_gene_file = "./temp_data/ptc_gene_list.bed"
ptc_transcript_file = "./temp_data/ptc_transcript_list.bed"

introns_file = "../source_data/introns.bed"

def get_feature(query_feature, gtf_file, file, fresh_run):
    if os.path.exists(file) and not fresh_run:
        features = gen.read_many_fields(file, "\t")
    else:
        features = [feature.split("\t") for feature in gen.run_process(["grep", "\t{0}\t".format(query_feature), gtf_file]).split('\n')[:-1]]
        with open(file, "w") as out:
            for feature in features:
                out.write("{0}\n".format("\t".join(feature)))
    return(features)

class get_info(object):
    def __init__(self, ptc):
        self.chr = ptc[0]
        self.transcript = ptc[3].split('.')[0]
        self.exon = int(ptc[3].split('.')[1])
        self.strand = ptc[5]
        self.pos = ptc[7]
        self.id = ptc[8]
        self.aa = ptc[9]
        self.ma = ptc[10]
        self.rel_cds_pos = ptc[11]


def get_ptc_list(ptc_file):
    ptc_list = {}
    for ptc in gen.read_many_fields(ptc_file, "\t")[1:]:
        ptc_info = get_info(ptc)
        ptc_list[ptc_info.id] = ptc

    return ptc_list

def get_cds_list(cds_fasta_file):
    cds_list = {}
    cds_names, cds_seqs = gen.read_fasta(cds_fasta_file)
    for i, name in enumerate(cds_names):
        cds_list[name] = cds_seqs[i]
    return cds_list


# exon_seq_list = collections.defaultdict(lambda: collections.defaultdict())
# exon_names, exon_seqs = gen.read_fasta(exon_fasta_file)
# for i, name in enumerate(exon_names):
#     t = name.split('.')[0]
#     e = int(name.split('.')[1])
#     if exon_seqs[i] in ["TAG", "TAA", "TGA"]:
#         e = 999
#     exon_seq_list[t][e] = exon_seqs[i]

# ptcs = gen.read_many_fields(ptc_file, "\t")
# features = gen.read_many_fields(gtf_file, "\t")
# exon_names, exon_seqs = gen.read_fasta(exon_fasta_file)


# # get a list of cds
# cds_list = {}
# for i, name in enumerate(cds_names):
#     cds_list[name] = cds_seqs[i]
#
# # retrive a list of ptc transcripts
# transcript_list = []
# for ptc in ptcs[1:]:
#     transcript_list.append(ptc[3].split('.')[0])
#
# # get a list of gene names for the transcripts
# gene_names = collections.defaultdict()
# for name in exon_names:
#     transcript_id = name.split('.')[0]
#     gene_id = name.split('(')[0].split('.')[-1]
#     gene_names[transcript_id] = gene_id



    # # get a file with all the transcripts
    # transcripts = get_feature('transcript', gtf_file, transcript_file, fresh_run)
    # with open(ptc_transcript_file, "w") as outfile:
    #     for transcript in transcripts:
    #         if len(transcript) > 0:
    #             info = transcript[8]
    #             transcript_id = re.findall('transcript_id\s"([^";]*)', info)[0]
    #             if transcript_id in transcript_list:
    #                 outfile.write("{0}\n".format("\t".join(transcript)))
    #
    #


def get_exon_list(cds_file):
    cds_entries = gen.read_many_fields(cds_file, "\t")
    exon_list = collections.defaultdict(lambda: collections.defaultdict())
    for exon in cds_entries:
        transcript_id = exon[3].split('.')[0]
        exon_id = int(exon[3].split('.')[1])
        type = exon[4]
        strand = exon[5]
        start = int(exon[1])
        stop = int(exon[2])

        if type == "stop_codon":
            exon_id = 999
        exon_list[transcript_id][exon_id] = [type, strand, start, stop]
    return exon_list


def get_intron_list(exon_list):
    intron_list = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict())))

    for i, transcript in enumerate(exon_list):
        if i:
            listed_exons = sorted(exon_list[transcript])
            for j, exon in enumerate(sorted(exon_list[transcript])):
                info = exon_list[transcript][exon]
                type = info[0]
                strand = info[1]

                focal_exon_start = info[2]
                focal_exon_stop = info[3]


                # the first exon
                if exon + 1 in listed_exons and exon - 1 not in listed_exons:
                    exon_type = "first_exon"
                    next_exon_start = exon_list[transcript][exon + 1][2]
                    next_exon_stop = exon_list[transcript][exon + 1][3]

                    if strand == "-":
                        intron_start = next_exon_stop
                        intron_end = focal_exon_start
                    else:
                        intron_start = focal_exon_stop
                        intron_end = next_exon_start

                    downstream_intron = [intron_start, intron_end]
                    upstream_intron = None

                # the last exon
                elif exon + 1 not in listed_exons and exon - 1 in listed_exons:
                    exon_type = "last_exon"
                    prev_exon_start = exon_list[transcript][exon - 1][2]
                    prev_exon_stop = exon_list[transcript][exon - 1][3]

                    if strand == "-":
                        intron_start = focal_exon_stop
                        intron_end = prev_exon_start
                    else:
                        intron_start = prev_exon_stop
                        intron_end = focal_exon_start

                    downstream_intron = None
                    upstream_intron = [intron_start, intron_end]

                elif exon + 1 in listed_exons and exon - 1 in listed_exons:
                    exon_type = "internal_exon"
                    next_exon_start = exon_list[transcript][exon + 1][2]
                    next_exon_stop = exon_list[transcript][exon + 1][3]
                    prev_exon_start = exon_list[transcript][exon - 1][2]
                    prev_exon_stop = exon_list[transcript][exon - 1][3]

                    if strand == "-":
                        upstream_intron_start = focal_exon_stop
                        upstream_intron_end = prev_exon_start
                        downstream_intron_start = next_exon_stop
                        downstream_intron_end = focal_exon_start
                    else:
                        upstream_intron_start = prev_exon_stop
                        upstream_intron_end = focal_exon_start
                        downstream_intron_start = focal_exon_stop
                        downstream_intron_end = next_exon_start

                    upstream_intron = [upstream_intron_start, upstream_intron_end]
                    downstream_intron = [downstream_intron_start, downstream_intron_end]
                else:
                    exon_type = "single_exon"
                    upstream_intron = None
                    downstream_intron = None

                # if exon_type != "first_exon" and exon_type != "single_exon":
                #     print("U Intron {0}".format(j+1), upstream_intron, upstream_intron[1] - upstream_intron[0])
                # print("Exon {0} {1}".format(exon, strand), [focal_exon_start, focal_exon_stop], focal_exon_stop-focal_exon_start)
                # print("{0}...{1}".format(exon_seq_list[transcript][exon][:6], exon_seq_list[transcript][exon][-6:]))
                # if exon_type != "last_exon" and exon_type != "single_exon" and exon_seq_list[transcript][exon] not in ["TAG", "TAA", "TGA"]:
                #     print("D Intron {0}".format(j), downstream_intron, downstream_intron[1] - downstream_intron[0])
                #
                # print("\n")

                intron_list[transcript][exon]["type"] = exon_type
                intron_list[transcript][exon]["upstream"] = upstream_intron
                intron_list[transcript][exon]["downstream"] = downstream_intron

    return intron_list

def get_intron_stats(intron_list, exon):
    internal_intron_count = len(intron_list)
    type = intron_list[exon]["type"]

    return internal_intron_count, type

def get_exon_rel_pos_list(cds_rel_pos_file):
    exon_rel_pos_list = {}
    ptcs = gen.read_many_fields(cds_rel_pos_file, "\t")
    for ptc in ptcs:
        t = ptc[3].split('.')[0]
        e = int(ptc[3].split('.')[1])
        id = ptc[8]
        pos = ptc[7]
        aa = ptc[9]
        ma = ptc[10]
        rel_pos = ptc[11]
        exon_rel_pos_list[id] = [t, e, pos, aa, ma, rel_pos]
    return exon_rel_pos_list

def get_gene_list(gene_list_file):
    gene_list = {}
    genes = gen.read_many_fields(gene_list_file, "\t")
    for gene in genes:
        print(gene)

def create_gene_file(gtf_file, genes_file, fresh_run):

    # get the genes from the gtf file
    gene_list = {}
    genes = get_feature('gene', gtf_file, genes_file, fresh_run)
    for gene in genes:
        # ensure the gene information exists
        if len(gene) > 0:
            # gene_info = gene[8]
            # gene_id = re.findall('gene_id\s"([^";]*)', gene_info)[0]
            # gene_list[gene_id] = gene
            with open(ptc_gene_file, "w") as outfile:
            # for transcript in transcript_list:
            #     transcript_gene_id = gene_names[transcript]
            #     gene_info = gene_list[transcript_gene_id]
                outfile.write("{0}\n".format("\t".join(gene)))


if create_files:
    create_gene_file(gtf_file, genes_file, fresh_run)

exon_list = get_exon_list(cds_file)
intron_list = get_intron_list(exon_list)
cds_list = get_cds_list(cds_fasta_file)
exon_rel_pos_list = get_exon_rel_pos_list(cds_rel_pos_file)
gene_list = get_gene_list(ptc_gene_file)



ptc_list = get_ptc_list(ptc_file)
for ptc in ptc_list:
    ptc = get_info(ptc_list[ptc])

    id = ptc.id
    transcript = ptc.transcript
    exon_number = ptc.exon
    cds_rel_pos = ptc.rel_cds_pos
    cds_seq = cds_list[ptc.transcript]
    cds_length = len(cds_seq) - 3   # ignore stop codon
    internal_intron_count, type = get_intron_stats(intron_list[ptc.transcript], ptc.exon)

    exon_rel_pos_info = exon_rel_pos_list[ptc.id]

    if transcript == exon_rel_pos_info[0] and exon_number == exon_rel_pos_info[1] and ptc.aa == exon_rel_pos_info[3] and ptc.ma == exon_rel_pos_info[4]:
        rel_exon_position = exon_rel_pos_info[5]

    exon_length = exon_list[transcript][exon_number][3] - exon_list[transcript][exon_number][2]

    if exon_length % 3 == 0:
        skip_frame = "in_frame"
    else:
        skip_frame = "out_of_frame"



    # print(transcript, exon_number, exon_length)
    # print(internal_intron_number, type)

# for transcript in intron_list:
#     print(transcript)
#     for exon in sorted(exon_list[transcript]):
#         if exon != 999:
#             exon_info = exon_list[transcript][exon]
#             upstream_intron = intron_list[transcript][exon]["upstream"]
#             downstream_intron = intron_list[transcript][exon]["downstream"]
#             type = intron_list[transcript][exon]["type"]
#
#             if(upstream_intron):
#                 print("U", upstream_intron[0], upstream_intron[1])
#             print(exon, exon_info[2], exon_info[3], type)
#             if downstream_intron:
#                 print("D", downstream_intron[0], downstream_intron[1])
#             print('\n')
