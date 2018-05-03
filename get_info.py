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
exons_file = './results/clean_run_2/clean_run_exons.bed'
introns_file = './temp_data/intron_list.bed'
transcript_file = "./temp_data/transcript_list.bed"

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


if create_files:

    # get the genes from the gtf file
    gene_list = {}
    genes = get_feature('gene', gtf_file, genes_file, fresh_run)
    for gene in genes:
        # ensure the gene information exists
        if len(gene) > 0:
            gene_info = gene[8]
            gene_id = re.findall('gene_id\s"([^";]*)', gene_info)[0]
            gene_list[gene_id] = gene
    with open(ptc_gene_file, "w") as outfile:
        for transcript in transcript_list:
            transcript_gene_id = gene_names[transcript]
            gene_info = gene_list[transcript_gene_id]
            outfile.write("{0}\n".format("\t".join(gene_info)))

    # get a file with all the transcripts
    transcripts = get_feature('transcript', gtf_file, transcript_file, fresh_run)
    with open(ptc_transcript_file, "w") as outfile:
        for transcript in transcripts:
            if len(transcript) > 0:
                info = transcript[8]
                transcript_id = re.findall('transcript_id\s"([^";]*)', info)[0]
                if transcript_id in transcript_list:
                    outfile.write("{0}\n".format("\t".join(transcript)))


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

for i, transcript in enumerate(exon_list):
    if i == 1070:
        listed_exons = sorted(exon_list[transcript])
        print(transcript)
        for i, exon in enumerate(sorted(exon_list[transcript])):
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

                upstream_intron = [intron_start, intron_end]
                downstream_intron = [0, 0]

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

                upstream_intron = [0, 0]
                downstream_intron = [intron_start, intron_end]

            elif exon + 1 in listed_exons and exon - 1 in listed_exons:
                exon_type = "internal_exon"
                next_exon_start = exon_list[transcript][exon + 1][2]
                next_exon_stop = exon_list[transcript][exon + 1][3]
                prev_exon_start = exon_list[transcript][exon - 1][2]
                prev_exon_stop = exon_list[transcript][exon - 1][3]

                if strand == "-":
                    upstream_intron_start = next_exon_stop
                    upstream_intron_end = focal_exon_start
                    downstream_intron_start = next_exon_stop
                    downstream_intron_end = focal_exon_start
                else:
                    upstream_intron_start = focal_exon_stop
                    upstream_intron_end = next_exon_start
                    downstream_intron_start = focal_exon_stop
                    downstream_intron_end = next_exon_start

                upstream_intron = [upstream_intron_start, upstream_intron_end]
                downstream_intron = [downstream_intron_start, downstream_intron_end]


            print("Exon {0}".format(exon), [focal_exon_start, focal_exon_stop], focal_exon_stop-focal_exon_start)
            print("Intron {0}".format(i+1), upstream_intron, upstream_intron[1] - upstream_intron[0])
