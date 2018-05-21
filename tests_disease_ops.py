from disease_ops import *
import generic as gen
import unittest

class Test_disease_ops(unittest.TestCase):

    def test_refactor_files(self):
        input_dir = "test_data/disease_ops/test_refactor_files/input_files"
        expected = "test_data/disease_ops/test_refactor_files/expected.txt"
        observed = "test_data/disease_ops/test_refactor_files/observed_1.txt"
        gen.remove_file(observed)
        output_dir = "test_data/disease_ops/test_refactor_files"
        output_prefix = "observed"
        refactor_files(input_dir, output_dir, output_prefix)
        gen.remove_file("{0}/{1}_readme.txt".format(output_dir, output_prefix))
        observed = gen.read_many_fields(observed, "\t")
        expected = gen.read_many_fields(expected, "\t")
        self.assertEqual(observed, expected)

    def test_check_line(self):
        outlist = [
            "Chromosome", "Start_position", "dbSNP_RS", "Strand", "Reference_Allele", "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
            "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification", "Variant_Type",
            "Refseq_mRNA_Id", "Refseq_prot_Id"
        ]
        input_file = "test_data/disease_ops/test_check_line/input.txt"
        input_lines = gen.read_many_fields(input_file, "\t")
        expected = [True, True, False, False, False, False, False, False]
        observed = []
        header = input_lines[0]
        for j, line in enumerate(input_lines[1:]):
            line_dict = {}
            for i, entry in enumerate(line):
                line_dict[header[i]] = entry
            print(j)
            line_pass, line_out = check_line(line_dict, outlist)
            observed.append(line_pass)
        self.assertEqual(observed, expected)