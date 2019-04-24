import generic as gen
import bed_ops as bo
import collections
import re
import os
import numpy as np
import pandas as pd
import random
from Bio.Seq import Seq
import difflib

file = "./results/experimental/minigenes.csv"

entries = [i for i in gen.read_many_fields(file, ",") if "ENST" in i[1]]

def print_result(result_text, result = False):
    if result:
        print("{0}... Pass".format(result_text))
    else:
        print("{0}... FAIL".format(result_text))


class Entry(object):
    def __init__(self, entry):
        self.construct_id = entry[0]
        self.transcript_id = entry[1]
        self.mutation_position = int(entry[2])
        self.codon_change = entry[3]
        self.remove_5 = int(entry[4])
        self.remove_3 = int(entry[5])
        self.upstream_exon = entry[6]
        self.upstream_intron = entry[7]
        self.exon = entry[8]
        self.downstream_intron = entry[9]
        self.downstream_exon = entry[10]
        self.downstream_is_final = entry[11]
        self.stop_codon = entry[12]
        self.exons_seq = entry[13]
        self.exons_length = int(entry[14])
        self.minigene_seq = entry[15]
        self.minigene_length = int(entry[16])
        self.protein_seq = entry[17]

    def test_generate_exon_seq(self):
        """
        Generate the sequence, make sure it is the one in the file
        """
        upstream_exon = self.upstream_exon
        if self.remove_5 != 0:
            upstream_exon = self.upstream_exon[self.remove_5:]
        downstream_exon = self.downstream_exon
        if self.remove_3 != 0:
            downstream_exon = self.downstream_exon[:-self.remove_3]
        stop = ""
        if self.downstream_is_final != "*":
            stop = self.stop_codon
        generated_seq = "ATG" + upstream_exon + self.exon + downstream_exon + stop

        print_result("Generated sequences match", self.exons_seq == generated_seq)

    def test_generate_minigene_seq(self):
        """
        Test if the minigene sequence is correct
        """

        upstream_exon = self.upstream_exon
        if self.remove_5 != 0:
            upstream_exon = self.upstream_exon[self.remove_5:]
        downstream_exon = self.downstream_exon
        if self.remove_3 != 0:
            downstream_exon = self.downstream_exon[:-self.remove_3]
        stop = ""
        if self.downstream_is_final != "*":
            stop = self.stop_codon

        minigene = "ATG" + upstream_exon + self.upstream_intron + self.exon + self.downstream_intron + downstream_exon + stop

        print_result("Minigene sequences match", self.minigene_seq == minigene)

        differences = []
        if self.minigene_seq != minigene:
            minigene_nts = list(self.minigene_seq)
            made_nts = list(minigene)
            for i, nt in enumerate(minigene_nts):
                if nt != made_nts[i]:
                    differences.append([i, "{0}|{1}".format(nt, made_nts[i])])
        if len(differences):
            print(differences)


    def test_minigene_length(self):
        """
        Test that the minigene lengths match
        """

        print_result("Minigene sequence lengths match", len(self.minigene_seq) == self.minigene_length)


    def test_exon_seq_length(self):
        """
        Test the length of the exon sequence matches
        """

        print_result("Exon sequence length matches", len(self.exons_seq) == self.exons_length)
        print_result("Exons sequence is a multiple of 3", len(self.exons_seq) % 3 == 0)


    def test_translation(self):
        """
        Test the translated sequence matches
        """

        translation = Seq(self.exons_seq).translate()
        print_result("Exons translation", self.protein_seq == translation)


    def test_mutation_position(self):
        """
        Test the position and codon of the mutation
        """

        upstream_exon = self.upstream_exon
        if self.remove_5 != 0:
            upstream_exon = self.upstream_exon[self.remove_5:]
        downstream_exon = self.downstream_exon
        if self.remove_3 != 0:
            downstream_exon = self.downstream_exon[:-self.remove_3]

        mutation_index = self.mutation_position - 1
        if "wt" in self.construct_id:
            codon = self.codon_change.split(">")[0].strip(" ")
        else:
            codon = self.codon_change.split(">")[1].strip(" ")

        seq_mutation_index = len(upstream_exon) + self.mutation_position -1
        if seq_mutation_index % 3 == 1:
            query_codon = self.exons_seq[seq_mutation_index-1:seq_mutation_index+2]
        elif seq_mutation_index % 3 == 2:
            query_codon = self.exons_seq[seq_mutation_index-2:seq_mutation_index+1]
        elif seq_mutation_index % 3 == 0:
            query_codon = self.exons_seq[seq_mutation_index-3:seq_mutation_index]

        print(codon, query_codon)
        print_result("Sample codon is the same as in the minigne", codon == query_codon)

for entry in entries:
    entry = Entry(entry)
    print("\n", entry.construct_id, entry.transcript_id)
    entry.test_generate_exon_seq()
    entry.test_exon_seq_length()
    entry.test_generate_minigene_seq()
    entry.test_minigene_length()
    entry.test_translation()
    entry.test_mutation_position()
