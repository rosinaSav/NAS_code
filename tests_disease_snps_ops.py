from disease_snps_ops import *
import generic as gen
import unittest

class Test_disease_snps_ops(unittest.TestCase):

    def test_get_ptcs_in_window(self):
        ptc_file = "test_data/disease_snps_ops/test_get_ptcs_in_window/ptc_list.txt"
        relative_positions_file = "test_data/disease_snps_ops/test_get_ptcs_in_window/relative_positions.txt"
        ptcs = gen.read_many_fields(ptc_file, "\t")
        relative_positions = gen.read_many_fields(relative_positions_file, "\t")
        ptc_list = {}
        relative_positions_list = {}
        for i, ptc in enumerate(ptcs):
            ptc[1], ptc[2], ptc[3] = int(ptc[1]), int(ptc[2]), int(ptc[3])
            ptc_list[i] = ptc
            rel_pos = relative_positions[i]
            rel_pos[1] = int(rel_pos[1])
            relative_positions_list[i] = rel_pos
        expected_file = "test_data/disease_snps_ops/test_get_ptcs_in_window/expected.txt"
        expected_list = gen.read_many_fields(expected_file, "\t")
        expected = {}
        ends = [5,3]
        for i in ends:
            expected[i] = {}
        for entry in expected_list:
            if entry[9] != '.' and int(entry[9]) in ends:
                required_entry = entry[:9]
                required_entry[1], required_entry[2], required_entry[3], required_entry[8] = int(required_entry[1]), int(required_entry[2]), int(required_entry[3]), int(required_entry[8])
                expected[int(entry[9])][int(entry[10])] = required_entry
        observed = get_ptcs_in_window(ptc_list, relative_positions_list, 4, 69)
        self.assertEqual(observed, expected)
