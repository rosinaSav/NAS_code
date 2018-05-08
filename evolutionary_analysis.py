'''
Authors: Liam Abrahams & Rosina Savisaar.
Evolutionary conservation analysis.
'''

import generic as gen
import os
import re
import collections
import random
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from io import StringIO


codon_dict = {'GGA': 'G', 'GGC': 'G', 'TTA': 'L', 'CTG': 'L', 'TTC': 'F', 'GAA': 'E', 'CTC': 'L', 'CGC': 'R', 'CAT': 'H', 'GTG': 'V', 'GGG': 'G',
 'GTC': 'V', 'CTA': 'L', 'GCC': 'A', 'TCG': 'S', 'CAC': 'H', 'AAG': 'K', 'AAT': 'N', 'CCG': 'P', 'TGT': 'C', 'CCC': 'P', 'CAA': 'Q',
 'CGT': 'R', 'GCG': 'A', 'ATA': 'I', 'ACA': 'T', 'TGC': 'C', 'CTT': 'L', 'ACT': 'T', 'CGA': 'R', 'TAC': 'Y', 'AGA': 'R', 'ACC': 'T',
 'GTA': 'V', 'ACG': 'T', 'TGG': 'W', 'AAC': 'N', 'CAG': 'Q', 'AGC': 'S', 'TCC': 'S', 'GCA': 'A', 'AGG': 'R', 'ATG': 'M', 'GAC': 'D',
 'TCA': 'S', 'TAT': 'Y', 'ATT': 'I', 'CCA': 'P', 'ATC': 'I', 'AAA': 'K', 'TTG': 'L', 'CCT': 'P', 'GGT': 'G', 'AGT': 'S', 'GCT': 'A',
 'TTT': 'F', 'TCT': 'S', 'GTT': 'V', 'GAG': 'E', 'CGG': 'R', 'GAT': 'D'}


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

def read_problematic_fasta(fasta):
    '''
    Read in a problematic fasta where there are line breaks.
    '''

    with open(fasta, "r") as file:
        lines = file.read()
        # names = [line.strip('>').strip('\n') for line in lines if line.startswith(">")]
        seqs = []
        names = []
        groups = lines.split('>')

        for group in groups:
            parts = group.split('\n')
            seq = "".join(parts[1:])
            names.append(parts[0])
            seqs.append(seq.strip('\n'))

    if len(names) != len(seqs):
        print('\nERROR: the number of names is not equal to the number of sequences...\n')
        raise Exception
    return names, seqs


def lists_to_dict(index_list, out_list):

    new_dict = {}
    for i, item in enumerate(index_list):
        new_dict[item] = out_list[i]
    return(new_dict)

def align_pairs(sequence_pairs, cds_fasta_file, macaque_cds_fasta):

    cds_names, cds_seqs = gen.read_fasta(cds_fasta_file)
    macaque_names, macaque_seqs = gen.read_fasta(macaque_cds_fasta)
    name_regex = re.compile('([^\s]+)')

    cds_seqs = lists_to_dict(cds_names, cds_seqs)
    macaque_seqs = lists_to_dict(macaque_names, macaque_seqs)


    for i, pair in enumerate(sequence_pairs):

        if i < 1:

            # get the sequences and convert to protein sequences
            cds_seq = cds_seqs[pair]
            if sequence_pairs[pair] in macaque_seqs:
                macaque_seq = macaque_seqs[sequence_pairs[pair]]

                cds_seq = Seq(cds_seq, IUPAC.unambiguous_dna)
                macaque_seq = Seq(macaque_seq, IUPAC.unambiguous_dna)

                cds_aa = cds_seq.translate()
                macaque_aa = macaque_seq.translate()

                # write to temporary file for muscle
                temp_file = "./temp_data/pair_fasta{0}.fasta".format(random.random())
                with open(temp_file, "w") as outfile:
                    outfile.write('>{0}\n{1}\n>{2}\n{3}\n'.format(pair, cds_aa, sequence_pairs[pair], macaque_aa))

                # align using muscle
                output_path = "./temp_data/muscle_output.fasta"
                alignment = MuscleCommandline("./binaries/muscle3.8.31_i86darwin64", input=temp_file, out=output_path)
                stdout, stderr = alignment()
                # get the string result
                with open(output_path, "r") as muscle_out:
                    str = "".join(muscle_out)
                gen.remove_file(temp_file)
                gen.remove_file(output_path)
                muscle_string = re.sub("([A-Z\-])\n([A-Z\-])","\\1\\2", str)
                aligned_prot_sequence = re.findall("^[A-Z\-]+(?=\n)",muscle_string, re.MULTILINE)

                # cds_nt = Seq(convert_protein_to_nt(cds_seq, aligned_prot_sequence[0])[:-3], IUPAC.unambiguous_dna)
                m_nt = Seq(convert_protein_to_nt(macaque_seq, aligned_prot_sequence[0]), IUPAC.unambiguous_dna)
                if m_nt[-3:] in ["TAA", "TAG", "TGA"]:
                    m_nt = m_nt[:-3]

                print( len(m_nt))
                # print(len(cds_nt), len(m_nt))
                # write to phylip file
                # multiple_alignment = MultipleSeqAlignment([SeqRecord(cds_nt, id=pair), SeqRecord(m_nt, id=sequence_pairs[pair])])
                # print(multiple_alignment)


def convert_protein_to_nt(nt_seq, aligned_prot_seq):

    print(len(aligned_prot_seq))
    print(len(nt_seq))
    nts = list(nt_seq)
    protein_sequence = list(aligned_prot_seq)

    alignment = []

    codon_counter = 0
    for i in protein_sequence:
        if i == "-":
            alignment.extend(["-", "-", "-"])
        else:
            chosen_nts = nts[codon_counter:codon_counter + 3]
            codon = "".join(chosen_nts)
            alignment.extend(chosen_nts)
        codon_counter += 3
    return "".join(alignment)


def filter_fasta_sequences(in_file, out_file, check_non_actg=None):

    names, seqs = read_problematic_fasta(in_file)
    pass_names, pass_seqs = [], []
    with open(out_file, "w") as out:
        for i, name in enumerate(names):
            seq = seqs[i]
            seq_pass = True
            if check_non_actg:
                for i in seq:
                    if i not in ["A", "C", "G", "T"]:
                        seq_pass = False
            if seq_pass and len(name):
                pass_names.append(name.split(' ')[0])
                pass_seqs.append(seq)

    gen.write_to_fasta(pass_names, pass_seqs, out_file)

def main():

    description = "Check whether genes associated with the PTCs are faster evolving."
    args = gen.parse_arguments(description, ["filter_sequences", "run_blast", "run_pairing"], flags = [0,1,2])
    filter_sequences, run_blast, run_pairing = args.filter_sequences, args.run_blast, args.run_pairing

    cds_fasta_file = "./results/clean_run_2/clean_run_CDS.extracted.2000.fasta"
    macaque_cds_fasta = "./source_data/Macaca_mulatta.MMUL_1.75.cds.all.fa"
    blast_output_file = "./temp_data/blast.txt"
    blast_dir = "blast"

    # filter the sequences to remove any not containing actg
    filtered_macaque_cds_fasta = "./source_data/Macaca_mulatta.MMUL_1.75.cds.all.filtered.fa"
    if filter_sequences or not os.path.isfile(filtered_macaque_cds_fasta):
        filter_fasta_sequences(macaque_cds_fasta, filtered_macaque_cds_fasta, check_non_actg=True)

    # run the blast
    if run_blast or not os.path.isfile(blast_output_file):
        blast_seqs(cds_fasta_file, blast_dir, filtered_macaque_cds_fasta, blast_output_file)

    # pair the sequences
    paired_sequences = "./temp_data/paired_sequences.txts"
    sequence_pairs, duplicates = pair_sequences(blast_output_file, run_pairing, paired_sequences)

    # align the pairs
    align_pairs(sequence_pairs, cds_fasta_file, filtered_macaque_cds_fasta)

if __name__ == "__main__":
    main()
