from bam_ops import *
import generic as gen
import os
import unittest

class Test_bam_ops(unittest.TestCase):

    def test_group_flags(self):
        input_bed = "test_data/bam_ops/test_group_flags/test_tabix.bed"
        observed = "test_data/bam_ops/test_group_flags/observed_test_group_flags.bed"
        gen.remove_file(observed)
        flag_start = 3
        group_flags(input_bed, observed, flag_start)
        expected = gen.read_many_fields("test_data/bam_ops/test_group_flags/expected_test_group_flags.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bam(self):
        bam_file = "test_data/bam_ops/test_intersect_bam/test_input_bam.bam"
        bed_file = "test_data/bam_ops/test_intersect_bam/test_input_bed.bed"
        observed_bam_output = "test_data/bam_ops/test_intersect_bam/observed_bam_intersect.bam"
        observed_bed_output = "test_data/bam_ops/test_intersect_bam/observed_bam_intersect.bed"
        expected_bed_output = "test_data/bam_ops/test_intersect_bam/expected_intersect_bed.bed"
        intersect_bed(bam_file, bed_file, output_file=observed_bam_output, intersect_bam=True)
        expected = gen.read_many_fields(expected_bed_output, "\t")
        #convert bam to bed to check correct output
        #use samtools to extract in the same format as bed
        samtools_args = ["samtools", "view", observed_bam_output]
        gen.run_process(samtools_args, file_for_output=observed_bed_output)
        observed = gen.read_many_fields(observed_bed_output, "\t")
        self.assertEqual(observed, expected)

    def test_intersect_bed_default(self):
        A_file = "test_data/bam_ops/test_intersect_bed_default/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_default/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_default/expected_test_intersect_bed_default.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_default/observed_test_intersect_bed_default.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_default_bedops(self):
        A_file = "test_data/bam_ops/test_intersect_bed_default_bedops/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_default_bedops/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_default_bedops/expected_test_intersect_bed_default_bedops.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_default_bedops/observed_test_intersect_bed_default_bedops.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, use_bedops = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_no_dups(self):
        A_file = "test_data/bam_ops/test_intersect_bed_no_dups/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_no_dups/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_no_dups/expected_test_intersect_bed_default_bedops.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_no_dups/observed_test_intersect_bed_no_dups.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count(self):
        A_file = "test_data/bam_ops/test_intersect_bed_hit_count/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_hit_count/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_hit_count/expected_test_intersect_bed_hit_count.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_hit_count/observed_test_intersect_bed_hit_count.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_hit_count_unsorted(self):
        A_file = "test_data/bam_ops/test_intersect_bed_hit_count_unsorted/test_intersect_bed_A_file_unsorted.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_hit_count_unsorted/test_intersect_bed_B_file_unsorted.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_hit_count_unsorted/expected_test_intersect_bed_hit_count_unsorted.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_hit_count_unsorted/observed_test_intersect_bed_hit_count_unsorted.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, hit_count = True, no_dups = False)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_overlap(self):
        A_file = "test_data/bam_ops/test_intersect_bed_overlap/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_overlap/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_overlap/expected_test_intersect_bed_overlap.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_overlap/observed_test_intersect_bed_overlap.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, overlap = 0.5)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand(self):
        A_file = "test_data/bam_ops/test_intersect_bed_force_strand/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_force_strand/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_force_strand/expected_test_intersect_bed_force_strand.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_force_strand/observed_test_intersect_bed_force_strand.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_force_strand_hit_count(self):
        A_file = "test_data/bam_ops/test_intersect_bed_force_strand_hit_count/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_force_strand_hit_count/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_force_strand_hit_count/expected_test_intersect_bed_force_strand_hit_count.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_force_strand_hit_count/observed_test_intersect_bed_force_strand_hit_count.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, force_strand = True, hit_count = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_write_both(self):
        A_file = "test_data/bam_ops/test_intersect_bed_write_both/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_write_both/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_write_both/expected_test_intersect_bed_write_both.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_write_both/observed_test_intersect_bed_write_both.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, write_both = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_sort_bed(self):
        infile = "test_data/bam_ops/test_sort_bed/test_intersect_bed_A_file_unsorted.bed"
        expected_file = "test_data/bam_ops/test_sort_bed/expected_test_intersect_bed_A_file.bed"
        observed_file = "test_data/bam_ops/test_sort_bed/observed_test_sort_bed.bed"
        gen.remove_file(observed_file)
        sort_bed(infile, observed_file)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_intersect(self):
        A_file = "test_data/bam_ops/test_intersect_bed_intersect/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_intersect/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_intersect/expected_test_intersect_bed_intersect.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_intersect/observed_test_intersect_bed_intersect.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, intersect = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)

    def test_intersect_bed_intersect_bedops(self):
        A_file = "test_data/bam_ops/test_intersect_bed_intersect_bedops/test_intersect_bed_A_file.bed"
        B_file = "test_data/bam_ops/test_intersect_bed_intersect_bedops/test_intersect_bed_B_file.bed"
        expected_file = "test_data/bam_ops/test_intersect_bed_intersect_bedops/expected_test_intersect_bed_intersect_bedops.bed"
        observed_file = "test_data/bam_ops/test_intersect_bed_intersect_bedops/observed_test_intersect_bed_intersect_bedops.bed"
        gen.remove_file(observed_file)
        intersect_bed(A_file, B_file, output_file = observed_file, no_dups = False, use_bedops = True, intersect = True)
        expected = gen.read_many_fields(expected_file, "\t")
        observed = gen.read_many_fields(observed_file, "\t")
        self.assertEqual(expected, observed)
