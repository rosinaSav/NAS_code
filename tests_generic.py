from generic import *
import unittest

class Test_generic(unittest.TestCase):

    def test_find_families_ensembl(self):
        ensembl_file = "test_data/generic/test_find_families/ensembl_file.txt"
        names = "ENST4, ENST1, ENST9, ENST2, ENST8, ENST3, ENST7, ENST5, ENST6"
        observed = "test_data/generic/test_find_families/observed.txt"
        expected = "test_data/generic/test_find_families/expected.txt"
        remove_file(observed)
        find_families_ensembl(ensembl_file, names, observed)
        with open(expected) as file:
            expected = "".join(file)
        with open(observed) as file:
            observed = "".join(file)
        self.assertEqual(expected, observed)
        
    def test_get_extension(self):
        file = "test.jpg"
        length = 3
        allowed = ["png", "gif", "jpg"]
        expected = "jpg"
        observed = get_extension(file, length, valid_list = allowed)
        self.assertEqual(expected, observed)

    def test_line_count(self):
        file = "test_data/generic/test_line_count/test_bed_second_go_expected.txt"
        expected = 3
        observed = line_count(file)
        self.assertEqual(expected, observed)

    def test_line_count2(self):
        file = "test_data/single-exon_wo_retrocopies_UCSC_ESE_positions_1000_test.txt"
        expected = 646
        observed = line_count(file)
        self.assertEqual(expected, observed)

    def test_list_to_dict(self):
        input_list = [[1, 5], [3, 8], [2, 0]]
        expected = {}
        expected[1] = 5
        expected[3] = 8
        expected[2] = 0
        observed = list_to_dict(input_list, 0, 1)
        self.assertEqual(observed, expected)

    def test_list_to_dict2(self):
        input_list = [[1, 4, 5], [3, 2, 8], [2, 9, 0]]
        expected = {}
        expected[1] = 5
        expected[3] = 8
        expected[2] = 0
        observed = list_to_dict(input_list, 0, 2)
        self.assertEqual(observed, expected)

    def test_list_to_dict_as_list(self):
        input_list = [[1, 4, 5], [3, 2, 8], [2, 9, 0], [3, 9, 15]]
        expected = {}
        expected[1] = [5]
        expected[3] = [8, 15]
        expected[2] = [0]
        observed = list_to_dict(input_list, 0, 2, as_list = True)
        self.assertEqual(observed, expected)

    def test_list_to_dict_as_list_uniquify(self):
        input_list = [[1, 4, 5], [3, 2, 8], [2, 9, 0], [3, 9, 15], [3, 60, 15]]
        expected = {}
        expected[1] = [5]
        expected[3] = [8, 15]
        expected[2] = [0]
        observed = list_to_dict(input_list, 0, 2, as_list = True, uniquify = True)
        self.assertEqual(observed, expected)

    def test_read_fasta(self):
        expected = [[],[]]
        expected[0] = ["1:17-25(+)", "2:0-12(+)", "1:21-30(-)"]
        expected[1] = ["CATAGACA", "GTCCCCCCCCAA", "AAATATGTC"]
        expected = tuple(expected)
        observed = read_fasta("test_data/generic/test_read_fasta/test_converted_fasta.fasta")
        self.assertEqual(expected, observed)
