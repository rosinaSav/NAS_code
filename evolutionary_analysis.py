import generic as gen
import os
import re
import collections
import random
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline
from io import StringIO

def blast_seqs(in_fasta_file, blast_dir, in_query_file, output_file):
    '''
    Blast two fasta files against each other
    '''

    print('Blasting sequences against each other...')
    # create the blast dir if it doesnt already exist
    gen.create_output_directories(blast_dir)
    # make the blast database
    db_name = "{0}/{1}_blast_db".format(blast_dir, in_fasta_file.split('/')[-1])
    gen.run_process(["makeblastdb", "-in", in_fasta_file, "-out", db_name, "-dbtype", "nucl"])

    # now run the second file against the first
    gen.run_process(["blastn", "-task", "blastn", "-db", db_name, "-query", in_query_file, "-outfmt", "10", "-out", output_file, "-num_threads", str(int((os.cpu_count()/2)-1)), "-evalue", "1e-04"])

def pair_sequences(blast_output_file, run_pairing, output_file):
    '''
    Pair the sequences from the local blast.
    '''

    blast_results = gen.read_many_fields(blast_output_file, ',')
    transcript_results = collections.defaultdict(lambda: [])
    for result in blast_results:
        transcript_results[result[1]].append(result)

    matches = {}
    used_matches = []
    for transcript in transcript_results:
        blasted_results = transcript_results[transcript]
        similarities = []
        for result in blasted_results:
            similarities.append(float(result[2]))
        best_match = blasted_results[similarities.index(max(similarities))]
        matches[transcript] = best_match
        used_matches.append(best_match[0])

    duplicates = [item for item in used_matches if used_matches.count(item) > 1]

    if run_pairing or not os.path.isfile(output_file):
        with open(output_file, "w") as outfile:
            for pair in matches:
                outfile.write("{0}\n".format("\t".join(matches[pair])))

    sequence_pairs = {}
    matches = gen.read_many_fields(output_file, "\t")
    for match in matches:
        sequence_pairs[match[1]] = match[0]

    return sequence_pairs, duplicates

def problematic_fasta(fasta):
    '''
    Read in a problematic fasta where there are line breaks.
    '''

    with open(fasta, "r") as file:
        lines = file.readlines()
        names = [line.strip('>').strip('\n') for line in lines if line.startswith(">")]
        seqs = []
        seq = ""
        for line in lines:
            if line.startswith('>'):
                if len(seq):
                    seqs.append(seq)
                seq = ""
            else:
                seq += line.strip("\n")

    return names, seqs


def align_pairs(sequence_pairs, cds_fasta_file, macaque_cds_fasta):

    cds_names, cds_seqs = gen.read_fasta(cds_fasta_file)
    macaque_names, macaque_seqs = problematic_fasta(macaque_cds_fasta)
    name_regex = re.compile('([^\s]+)')
    macaque_names = [re.search(name_regex, name).group(0) for name in macaque_names]

    for i, pair in enumerate(sequence_pairs):

        if i < 5:
            # get the sequences and convert to protein sequences
            cds_seq = cds_seqs[cds_names.index(pair)]
            macaque_seq = macaque_seqs[macaque_names.index(sequence_pairs[pair])]

            cds_seq = Seq(cds_seq, IUPAC.unambiguous_dna)
            macaque_seq = Seq(macaque_seq, IUPAC.unambiguous_dna)

            cds_aa = cds_seq.translate()
            macaque_aa = macaque_seq.translate()

            # write to temporary file for muscle
            temp_file = "./temp_data/pair_fasta{0}.fasta".format(random.random())
            with open(temp_file, "w") as outfile:
                outfile.write('>{0}\n{1}\n>{2}\n{3}\n'.format(pair, cds_aa, sequence_pairs[pair], macaque_aa))

            # set the path for the muscle output
            output_path = "./temp_data/muscle_output.fasta"
            # align using muscle
            align = MuscleCommandline("./binaries/muscle3.8.31_i86darwin64", input=temp_file, out=output_path)
            stdout, stderr = align()
            # get the string result
            with open(output_path, "r") as mout:
                str = "".join(mout)
            gen.remove_file(temp_file)
            gen.remove_file(output_path)
            muscle_string = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", str)
            aligned_prot_sequences = re.findall("^[A-Z\-]+(?=\n)",muscle_string, re.MULTILINE)
            print(aligned_prot_sequences)

def main():



    description = "Check whether genes associated with the PTCs are faster evolving."
    args = gen.parse_arguments(description, ["run_blast", "run_pairing"], flags = [0, 1])
    run_blast, run_pairing = args.run_blast, args.run_pairing

    cds_fasta_file = "./results/clean_run_2/clean_run_CDS.extracted.2000.fasta"
    macaque_cds_fasta = "./source_data/Macaca_mulatta.MMUL_1.75.cds.all.fa"
    blast_output_file = "./temp_data/blast.txt"
    blast_dir = "blast"

    if run_blast or not os.path.isfile(blast_output_file):
        blast_seqs(cds_fasta_file, blast_dir, macaque_cds_fasta, blast_output_file)

    paired_sequences = "./temp_data/paired_sequences.txts"
    sequence_pairs, duplicates = pair_sequences(blast_output_file, run_pairing, paired_sequences)

    align_pairs(sequence_pairs, cds_fasta_file, macaque_cds_fasta)

if __name__ == "__main__":
    main()
