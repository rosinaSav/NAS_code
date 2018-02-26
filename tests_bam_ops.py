from bam_ops import *
import generic as gen
import os
import unittest

class Test_bam_ops(unittest.TestCase):

    def test_bam_flag_filter_improper_paired_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_improper_paired_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_improper_paired_reads/expected_flag_filtered_improper_paired_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_improper_paired_reads/observed_flag_improper_paired_reads_filtered_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_improper_paired_reads/observed_flag_improper_paired_reads_filtered.sam"
        bam_flag_filter(input_bam, observed, get_paired_reads=True, get_improper_paired_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_mapped_proper_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_mapped_proper_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_mapped_proper_reads/expected_flag_filtered_mapped_proper_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_mapped_proper_reads/observed_flag_filtered_mapped_proper_reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_mapped_proper_reads/observed_flag_filtered_mapped_proper_reads.sam"
        bam_flag_filter(input_bam, observed, get_mapped_reads=True, get_paired_reads=True, get_proper_paired_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_mapped_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_mapped_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_mapped_reads/expected_flag_filtered_mapped_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_mapped_reads/observed_flag_filtered_mapped_reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_mapped_reads/observed_flag_filtered_mapped_reads.sam"
        bam_flag_filter(input_bam, observed, get_mapped_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_mate_mapped_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_mate_mapped_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_mate_mapped_reads/expected_flag_filtered_mate_mapped_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_mate_mapped_reads/observed_flag_filtered_mate_mapped_reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_mate_mapped_reads/observed_flag_filtered_mate_mapped_reads.sam"
        bam_flag_filter(input_bam, observed, get_paired_reads=True, get_mate_mapped_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_mate_unmapped_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_mate_unmapped_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_mate_unmapped_reads/expected_flag_filtered_mate_unmapped_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_mate_unmapped_reads/observed_flag_filtered_mate_unmapped_reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_mate_unmapped_reads/observed_flag_filtered_mate_unmapped_reads.sam"
        bam_flag_filter(input_bam, observed, get_paired_reads=True, get_mate_unmapped_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_paired_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_paired_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_paired_reads/expected_flag_filtered_paired_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_paired_reads/observed_flag_paired_reads_filtered_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_paired_reads/observed_flag_paired_reads_filtered.sam"
        bam_flag_filter(input_bam, observed, get_paired_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_proper_paired_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_proper_paired_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_proper_paired_reads/expected_flag_filtered_proper_paired_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_proper_paired_reads/observed_flag_proper_paired_reads_filtered_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_proper_paired_reads/observed_flag_proper_paired_reads_filtered.sam"
        bam_flag_filter(input_bam, observed, get_paired_reads=True, get_proper_paired_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_unmapped_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_unmapped_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_unmapped_reads/expected_flag_filtered_unmapped_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_unmapped_reads/observed_flag_filtered_unmapped_reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_unmapped_reads/observed_flag_filtered_unmapped_reads.sam"
        bam_flag_filter(input_bam, observed, get_unmapped_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_unpaired_mapped_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_unpaired_mapped_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_unpaired_mapped_reads/expected_flag_filtered_unpaired_mapped_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_unpaired_mapped_reads/observed_flag_filtered_unpaired__reads_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_unpaired_mapped_reads/observed_flag_filtered_unpaired_mapped_reads.sam"
        bam_flag_filter(input_bam, observed, get_unpaired_reads= True, get_mapped_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_flag_filter_unpaired_reads(self):
        input_bam = "test_data/bam_ops/test_bam_flag_filter_unpaired_reads/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_flag_filter_unpaired_reads/expected_flag_filtered_unpaired_reads.sam"
        observed = "test_data/bam_ops/test_bam_flag_filter_unpaired_reads/observed_flag_unpaired_reads_filtered_bam.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_flag_filter_unpaired_reads/observed_flag_unpaired_reads_filtered.sam"
        bam_flag_filter(input_bam, observed, get_unpaired_reads=True)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_nm_filter(self):
        input_bam = "test_data/bam_ops/test_bam_nm_filter/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_nm_filter/expected_bam_nm_filter.sam"
        observed = "test_data/bam_ops/test_bam_nm_filter/observed_bam_nm_filter.bam"
        observed_sam_file = "test_data/bam_ops/test_bam_nm_filter/observed_bam_nm_filter.sam"
        bam_nm_filter(input_bam, observed, nm_less_equal_to=6)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_file)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_file, "\t")
        self.assertEqual(expected, observed)

    def test_bam_xt_filter(self):
        input_bam = "test_data/bam_ops/test_bam_xt_filter/input_bam.bam"
        expected = "test_data/bam_ops/test_bam_xt_filter/expected_bam_xt_filter.sam"
        observed = "test_data/bam_ops/test_bam_xt_filter/observed_bam_xt_filter.bam"
        observed_sam_file = "test_data/bam_ops/test_bam_xt_filter/observed_bam_xt_filter.sam"
        bam_xt_filter(input_bam, observed, xt_filter="U")
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_file)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_file, "\t")
        self.assertEqual(expected, observed)


    def test_bam_quality_filter(self):
        input_bam = "test_data/bam_ops/test_bam_quality_filter/test_bam.bam"
        expected = "test_data/bam_ops/test_bam_quality_filter/expected_bam_quality_filter.sam"
        observed = "test_data/bam_ops/test_bam_quality_filter/observed_bam_quality_filter.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_quality_filter/observed_bam_quality_filter.sam"
        expected = gen.read_many_fields(expected, "\t")
        bam_quality_filter(input_bam, observed, quality_greater_than_equal_to=100, quality_less_than_equal_to=200)
        #convert bam to sam to check correct output
        #use samtools to extract in the same format as sam
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_quality_filter_greater_than(self):
        input_bam = "test_data/bam_ops/test_bam_quality_filter_greater_than/test_bam.bam"
        expected = "test_data/bam_ops/test_bam_quality_filter_greater_than/expected_bam_quality_filter_greater_than.sam"
        observed = "test_data/bam_ops/test_bam_quality_filter_greater_than/observed_bam_quality_filter_greater_than.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_quality_filter_greater_than/observed_bam_quality_filter_greater_than.sam"
        expected = gen.read_many_fields(expected, "\t")
        bam_quality_filter(input_bam, observed, quality_greater_than_equal_to=200)
        #convert bam to sam to check correct output
        #use samtools to extract in the same format as sam
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_bam_quality_filter_less_than(self):
        input_bam = "test_data/bam_ops/test_bam_quality_filter_less_than/test_bam.bam"
        expected = "test_data/bam_ops/test_bam_quality_filter_less_than/expected_bam_quality_filter_less_than.sam"
        observed = "test_data/bam_ops/test_bam_quality_filter_less_than/observed_bam_quality_filter_less_than.bam"
        observed_sam_output = "test_data/bam_ops/test_bam_quality_filter_less_than/observed_bam_quality_filter_less_than.sam"
        expected = gen.read_many_fields(expected, "\t")
        bam_quality_filter(input_bam, observed, quality_less_than_equal_to=250)
        #convert bam to sam to check correct output
        #use samtools to extract in the same format as sam
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        observed = gen.read_many_fields(observed_sam_output, "\t")
        self.assertEqual(expected, observed)

    def test_compare_PSI(self):
        SNPs = "test_data/bam_ops/test_compare_PSI/SNPs.bed"
        bam_folder = "test_data/bam_ops/test_compare_PSI/bam_folder"
        expected = gen.read_many_fields("test_data/bam_ops/test_compare_PSI/expected.txt", "\t")
        observed = "test_data/bam_ops/test_compare_PSI/observed.txt"
        gen.remove_file(observed)
        compare_PSI(SNPs, bam_folder, observed)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_count_junction_reads(self):
        sam = "test_data/bam_ops/test_count_junction_reads/reads.sam"
        junctions = {"chr1": {51: {54: {"exon": ["ENST1.1", "ENST1.2"], "type": ["incl", "incl"]},
                             59: {"exon": ["ENST1.2"], "type": ["skip"]}},
                    56: {59: {"exon": ["ENST1.2", "ENST1.3"], "type": ["incl", "incl"]},
                         67: {"exon": ["ENST1.3"], "type": ["skip"]}},
                    63: {67: {"exon": ["ENST1.3", "ENST1.4"], "type": ["incl", "incl"]}},
                    76: {79: {"exon": ["ENST1.5", "ENST1.6"], "type": ["incl", "incl"]}}},
                    "chr2": {34: {37: {"exon": ["ENST2.1", "ENST2.2"], "type": ["incl", "incl"]}},
                    46: {50: {"exon": ["ENST2.3", "ENST2.4"], "type": ["incl", "incl"]},
                         57: {"exon": ["ENST2.4"], "type": ["skip"]}},
                    54: {57: {"exon": ["ENST2.4", "ENST2.5"], "type": ["incl", "incl"]},
                         62: {"exon": ["ENST2.5"], "type": ["skip"]}},
                    59: {62: {"exon": ["ENST2.5", "ENST2.6"], "type": ["incl", "incl"]}}}}
        expected = gen.read_many_fields("test_data/bam_ops/test_count_junction_reads/expected.txt", "\t")
        observed = "test_data/bam_ops/test_count_junction_reads/observed.txt"
        gen.remove_file(observed)
        count_junction_reads(sam, junctions, observed, 1000)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

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

    def test_map_intron_from_cigar(self):
        cigars = [["75M", 55], ["21M1I3M1I4M3I42M", 845], ["10M2N65M", 48], ["4I18M16N53M", 40], ["10M2N5M3N3M3N4M", 48], ["18M3I16N57M", 40]]
        expected = [None, None, [[56, 59]], [[56, 73]], [[56, 59], [63, 67], [69, 73]], [[56, 73]]]
        observed = [map_intron_from_cigar(i[0], i[1]) for i in cigars]
        self.assertEqual(expected, observed)

    def test_read_exon_junctions(self):
        junctions_file = "test_data/bam_ops/test_read_exon_junctions/junctions.bed"
        expected = {"chr1": {51: {54: {"exon": ["ENST1.1", "ENST1.2"], "type": ["incl", "incl"]},
                             59: {"exon": ["ENST1.2"], "type": ["skip"]}},
                    56: {59: {"exon": ["ENST1.2", "ENST1.3"], "type": ["incl", "incl"]},
                         67: {"exon": ["ENST1.3"], "type": ["skip"]}},
                    63: {67: {"exon": ["ENST1.3", "ENST1.4"], "type": ["incl", "incl"]}},
                    76: {79: {"exon": ["ENST1.5", "ENST1.6"], "type": ["incl", "incl"]}}},
                    "chr2": {34: {37: {"exon": ["ENST2.1", "ENST2.2"], "type": ["incl", "incl"]}},
                    46: {50: {"exon": ["ENST2.3", "ENST2.4"], "type": ["incl", "incl"]},
                         57: {"exon": ["ENST2.4"], "type": ["skip"]}},
                    54: {57: {"exon": ["ENST2.4", "ENST2.5"], "type": ["incl", "incl"]},
                         62: {"exon": ["ENST2.5"], "type": ["skip"]}},
                    59: {62: {"exon": ["ENST2.5", "ENST2.6"], "type": ["incl", "incl"]}}}}
        observed = read_exon_junctions(junctions_file)
        self.assertEqual(expected, observed)

    def test_merge_bams(self):
        input_bam1 = "test_data/bam_ops/test_merge_bams/input1.bam"
        input_bam2 = "test_data/bam_ops/test_merge_bams/input2.bam"
        input_list = [input_bam1, input_bam2]
        expected = "test_data/bam_ops/test_merge_bams/expected_merge_bams.sam"
        observed = "test_data/bam_ops/test_merge_bams/observed_merge_bams.bam"
        observed_sam_output = "test_data/bam_ops/test_merge_bams/observed_merge_bams.sam"
        merge_bams(input_list, observed)
        #convert bam to sam to check correct output
        samtools_args = ["samtools", "view", observed]
        gen.run_process(samtools_args, file_for_output=observed_sam_output)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed_sam_output, "\t")
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
