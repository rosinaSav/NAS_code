from bed_ops import *
import generic as gen
import unittest

class Test_bed_ops(unittest.TestCase):

    def test_extract_exons(self):
        gtf = "test_data/test_extract_exons.gtf"
        observed = "test_data/test_extract_exons_observed.bed"
        gen.remove_file(observed)
        extract_exons(gtf, observed)
        expected = gen.read_many_fields("test_data/test_extract_exons_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)
        