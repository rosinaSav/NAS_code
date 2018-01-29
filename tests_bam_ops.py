from bam_ops import *
import generic as gen
import os
import unittest

class Test_bam_ops(unittest.TestCase):

    def test_extract_exons(self):
        gtf = "test_data/test_extract_exons.gtf"
        observed = "test_data/test_extract_exons_observed.bed"
        gen.remove_file(observed)
        extract_exons(gtf, observed)
        expected = gen.read_many_fields("test_data/test_extract_exons_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_group_flags(self):
        input_bed = "test_data/test_tabix.bed"
        observed = "test_data/test_group_flags_observed.bed"
        gen.remove_file(observed)
        flag_start = 3
        group_flags(input_bed, observed, flag_start)
        expected = gen.read_many_fields("test_data/test_group_flags_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_default(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_default_expected.bed"
        observed_file = "test_data/test_intersect_bed_default_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

##    def test_intersect_bed_default_bedops(self):
##        A_file = "test_data/test_intersect_bed_A_file.bed"
##        B_file = "test_data/test_intersect_bed_B_file.bed"
##        expected_file = "test_data/test_intersect_bed_default_bedops_expected.bed"
##        observed_file = "test_data/test_intersect_bed_default_bedops_observed.bed"
##        gen.remove_file(observed_file)
##        intersect_bed(A_file, B_file, output_file = observed_file, use_bedops = True)
##        expected = gen.read_many_fields(expected_file, "\t")
##        observed = gen.read_many_fields(observed_file, "\t")
##        print("\n")
##        for i in expected:
##            print(i)
##        print("\n")
##        for i in observed:
##            print(i)
##        self.assertEqual(expected, observed)

    def test_intersect_bed_no_dups(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_default_bedops_expected.bed"
        observed_file = "test_data/test_intersect_bed_no_dups_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_hit_count_expected.bed"
        observed_file = "test_data/test_intersect_bed_hit_count_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count_unsorted(self):
        A_file = "test_data/test_intersect_bed_A_file_unsorted.bed"
        B_file = "test_data/test_intersect_bed_B_file_unsorted.bed"
        expected_file = "test_data/test_intersect_bed_hit_count_unsorted_expected.bed"
        observed_file = "test_data/test_intersect_bed_hit_count_unsorted_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_overlap(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_overlap_expected.bed"
        observed_file = "test_data/test_intersect_bed_overlap_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, overlap = 0.5)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_force_strand_expected.bed"
        observed_file = "test_data/test_intersect_bed_force_strand_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand_hit_count(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_force_strand_hit_count_expected.bed"
        observed_file = "test_data/test_intersect_bed_force_strand_hit_count_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True, hit_count = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_write_both(self):
        A_file = "test_data/test_intersect_bed_A_file.bed"
        B_file = "test_data/test_intersect_bed_B_file.bed"
        expected_file = "test_data/test_intersect_bed_write_both_expected.bed"
        observed_file = "test_data/test_intersect_bed_write_both_observed.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, write_both = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

        
