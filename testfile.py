import generic as gen
import re
import bed_ops as bedops
gen.extract_head_of_file('./source_data/Homo_sapiens.GRCh37.87.gtf', 5000)
# gen.extract_head_of_file('./source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa', 500)


# def extract_exons(gtf, bed):
#     '''Given a GTF file, extract exon coordinates and write them to .bed.
#     EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf", "../source_data/Homo_sapiens.GRCh37.87_exons.bed")'''
#     gen.remove_file('./source_data/Homo_sapiens.GRCh37.87.extracted.500.bed')
#     #extract exons from GTF
#     exons = gen.run_process(["grep", "\texon\t", gtf])
#     #filter down to only protein-coding ones
#     exons = gen.run_process(["grep", "gene_biotype \"protein_coding\""], input_to_pipe = exons)
#     #split lines
#     exons = [i.split("\t") for i in exons.split("\n")]
#     #format as .bed. Switch to base 0.
#     exons = [["chr{0}".format(i[0]), int(i[3]) - 1, i[4], i[8], ".", i[6]] for i in exons if len(i) >= 8]
#     #pre-compile regex
#     trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
#     exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
#     #extract transcript IDs and exon numbers
#     for pos, exon in enumerate(exons):
#         to_parse = exon[3]
#         trans = re.search(trans_regex, to_parse).group(0)
#         exon_no = re.search(exon_no_regex, to_parse).group(0)
#         exons[pos][3] = "{0}.{1}".format(trans, exon_no)
#     #write to bed
#     with open(bed, "w") as file:
#         for exon in exons:
#             file.write("{0}\n".format("\t".join([str(i) for i in exon])))


def extract_exons(gtf, bed):
    '''Given a GTF file, extract exon coordinates and write them to .bed.
    EX.: extract_exons("../source_data/Homo_sapiens.GRCh37.87.gtf", "../source_data/Homo_sapiens.GRCh37.87_exons.bed")'''
    gen.remove_file('./source_data/Homo_sapiens.GRCh37.87.extracted.500.bed')
    #extract exons from GTF
    exons = gen.run_process(["grep", "\tCDS\t", gtf])
    #filter down to only protein-coding ones
    exons = gen.run_process(["grep", "gene_biotype \"protein_coding\""], input_to_pipe = exons)
    # exons = gen.run_process(["grep", "transcript_id \"ENST00000437963\""], input_to_pipe = exons)
    #split lines
    exons = [i.split("\t") for i in exons.split("\n")]
    #format as .bed. Switch to base 0.
    exons = [["{0}".format(i[0]), int(i[3]) - 1, i[4], i[8], ".", i[6]] for i in exons if len(i) >= 8]
    #pre-compile regex
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
    #extract transcript IDs and exon numbers
    for pos, exon in enumerate(exons):
        to_parse = exon[3]
        trans = re.search(trans_regex, to_parse).group(0)
        exon_no = re.search(exon_no_regex, to_parse).group(0)
        exons[pos][3] = "{0}.{1}".format(trans, exon_no)
    #write to bed
    with open(bed, "w") as file:
        for exon in exons:
            file.write("{0}\n".format("\t".join([str(i) for i in exon])))



def bed_to_fasta():
    genome_fasta = './source_data/Genomes/hg37/Homo_sapiens.GRCh37.dna.primary_assembly.fa'
    bed_file = './source_data/Homo_sapiens.GRCh37.87.extracted.5000.bed'
    fasta_file = './source_data/test_output.fasta'
    bedops.fasta_from_intervals(bed_file, fasta_file, genome_fasta, force_strand = True, names = True)

extract_exons('./source_data/Homo_sapiens.GRCh37.87.extracted.5000.gtf', './source_data/Homo_sapiens.GRCh37.87.extracted.5000.bed')
bed_to_fasta()
