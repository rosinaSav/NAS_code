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
cds_file = './results/clean_run_2/clean_run_CDS.fasta'
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


ptcs = gen.read_many_fields(ptc_file, "\t")
features = gen.read_many_fields(gtf_file, "\t")
exon_names, exon_seqs = gen.read_fasta(exon_fasta_file)
cds_names, cds_seqs = gen.read_fasta(cds_file)

# get a list of cds
cds_list = {}
for i, name in enumerate(cds_names):
    cds_list[name] = cds_seqs[i]

# retrive a list of ptc transcripts
transcript_list = []
for ptc in ptcs[1:]:
    transcript_list.append(ptc[3].split('.')[0])

# get a list of gene names for the transcripts
gene_names = collections.defaultdict()
for name in exon_names:
    transcript_id = name.split('.')[0]
    gene_id = name.split('(')[0].split('.')[-1]
    gene_names[transcript_id] = gene_id


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

introns = gen.read_many_fields(introns_file, "\t")
for intron in introns:
    print(intron)
