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

    def test_junction_raw_counts_to_bed(self):
        rawfile = "counts_file.txt"
        expected = "test_data/disease_ops/test_junction_raw_counts_to_bed/expected.bed"
        expected_samples = "test_data/disease_ops/test_junction_raw_counts_to_bed/expected_samples.bed"
        observed = "test_data/disease_ops/test_junction_raw_counts_to_bed/observed.bed"
        observed_samples = "test_data/disease_ops/test_junction_raw_counts_to_bed/observed_samples.bed"
        dir = "test_data/disease_ops/test_junction_raw_counts_to_bed"
        gen.remove_file(observed)
        junction_raw_counts_to_bed(rawfile, dir, observed, observed_samples)
        observed = gen.read_many_fields(observed, "\t")
        expected = gen.read_many_fields(expected, "\t")
        self.assertEqual(observed, expected)
        observed_samples = gen.read_many_fields(observed_samples, "\t")
        expected_samples = gen.read_many_fields(expected_samples, "\t")
        self.assertEqual(observed_samples, expected_samples)

    def test_raw_counts_to_samples(self):
        sample_file = "test_data/disease_ops/test_raw_counts_to_samples/samples.bed"
        intersect_file = "test_data/disease_ops/test_raw_counts_to_samples/IN.intersect.bed"
        output_dir = "test_data/disease_ops/test_raw_counts_to_samples/IN"
        expected1 = "test_data/disease_ops/test_raw_counts_to_samples/expected1.bed"
        expected2 = "test_data/disease_ops/test_raw_counts_to_samples/expected2.bed"
        observed1 = "test_data/disease_ops/test_raw_counts_to_samples/IN/observed1.bed"
        observed2 = "test_data/disease_ops/test_raw_counts_to_samples/IN/observed2.bed"
        gen.remove_file(observed1)
        gen.remove_file(observed2)
        raw_counts_to_samples(intersect_file, sample_file, output_dir)
        observed1 = gen.read_many_fields(observed1, "\t")
        observed2 = gen.read_many_fields(observed2, "\t")
        expected1 = gen.read_many_fields(expected1, "\t")
        expected2 = gen.read_many_fields(expected2, "\t")
        self.assertEqual(observed1, expected1)
        self.assertEqual(observed2, expected2)

    def test_get_ptcs_in_window(self):
        ptc_file = "test_data/disease_ops/test_get_ptcs_in_window/ptc_list.txt"
        relative_positions_file = "test_data/disease_ops/test_get_ptcs_in_window/relative_positions.txt"
        ptcs = gen.read_many_fields(ptc_file, "\t")
        relative_positions = gen.read_many_fields(relative_positions_file, "\t")
        ptc_list = {}
        relative_positions_list = {}
        for i, ptc in enumerate(ptcs):
            ptc_list[i] = ptc
            relative_positions_list[i] = relative_positions[i]
        expected_file = "test_data/disease_ops/test_get_ptcs_in_window/expected.txt"
        expected_list = gen.read_many_fields(expected_file, "\t")
        expected = {}
        ends = [5,3]
        for i in ends:
            expected[i] = {}
        for entry in expected_list:
            if entry[8] in ends:
                expected[entry[8]][entry[9]] = entry[:8]
        observed = get_ptcs_in_window(ptc_list, relative_positions_list, 3, 69)
        self.assertEqual(observed, expected)
