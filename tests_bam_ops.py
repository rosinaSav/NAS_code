from bam_ops import *
import generic as gen
import os
import unittest

class Test_bam_ops(unittest.TestCase):

    def test_group_flags(self):
        input_bed = "test_data/test_tabix.bed"
        observed = "test_data/test_group_flags_observed.bed"
        gen.remove_file(observed)
        flag_start = 3
        group_flags(input_bed, observed, flag_start)
        expected = gen.read_many_fields("test_data/test_group_flags_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)
        
