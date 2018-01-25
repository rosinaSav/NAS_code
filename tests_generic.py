from generic import *
import unittest

class Test_generic(unittest.TestCase):

    def test_get_extension(self):
        file = "test.jpg"
        length = 3
        allowed = ["png", "gif", "jpg"]
        expected = "jpg"
        observed = get_extension(file, length, valid_list = allowed)
        self.assertEqual(expected, observed)

    def test_line_count(self):
        file = "test_data/test_bed_second_go_expected.txt"
        expected = 3
        observed = line_count(file)
        self.assertEqual(expected, observed)

    def test_line_count2(self):
        file = "test_data/single-exon_wo_retrocopies_UCSC_ESE_positions_1000_test.txt"
        expected = 646
        observed = line_count(file)
        self.assertEqual(expected, observed)
