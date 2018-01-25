import generic as gen
import re

def extract_exons(gtf, bed):
    '''Given a GTF file, extract exon coordinates and write them to .bed.'''
    exons = gen.run_process(["grep", "\texon\t", gtf])
    exons = gen.run_process(["grep", "gene_biotype \"protein_coding\""], input_to_pipe = exons)
    exons = [i.split("\t") for i in exons.split("\n")]
    exons = [["chr{0}".format(i[0]), int(i[3]) - 1, i[4], i[8], ".", i[6]] for i in exons if len(i) >= 8]
    trans_regex = re.compile("(?<=transcript_id \")ENST[0-9]*")
    exon_no_regex = re.compile("(?<=exon_number \")[0-9]*")
    for pos, exon in enumerate(exons):
        to_parse = exon[3]
        trans = re.search(trans_regex, to_parse).group(0)
        exon_no = re.search(exon_no_regex, to_parse).group(0)
        exons[pos][3] = "{0}.{1}".format(trans, exon_no)
    with open(bed, "w") as file:
        for exon in exons:
            file.write("{0}\n".format("\t".join([str(i) for i in exon])))
        

