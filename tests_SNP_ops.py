import generic as gen
import re
import numpy as np
from SNP_ops import *
import unittest

class Test_SNP_ops(unittest.TestCase):

    def test_filter_by_snp_type(self):
        input_snps = "test_data/snp_ops/test_filter_by_snp_type/input_snps.bed"
        expected = "test_data/snp_ops/test_filter_by_snp_type/expected_snps.bed"
        observed = "test_data/snp_ops/test_filter_by_snp_type/observed_snps.bed"
        gen.remove_file(observed)
        filter_by_snp_type(input_snps, observed, "non")
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_filter_motif_SNPs(self):
        motifs = "test_data/snp_ops/test_filter_motif_snps/ESEs.txt"
        fasta = "test_data/snp_ops/test_filter_motif_snps/CDS.fasta"
        SNPs = "test_data/snp_ops/test_filter_motif_snps/snps.bed"
        expected = "test_data/snp_ops/test_filter_motif_snps/expected.txt"
        observed = "test_data/snp_ops/test_filter_motif_snps/observed.txt"
        gen.remove_file(observed)
        filter_motif_SNPs(fasta, SNPs, motifs, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_filter_motif_SNPs_complement(self):
        motifs = "test_data/snp_ops/test_filter_motif_snps_complement/ESEs.txt"
        fasta = "test_data/snp_ops/test_filter_motif_snps_complement/CDS.fasta"
        SNPs = "test_data/snp_ops/test_filter_motif_snps_complement/snps.bed"
        expected = "test_data/snp_ops/test_filter_motif_snps_complement/expected.txt"
        observed = "test_data/snp_ops/test_filter_motif_snps_complement/observed.txt"
        gen.remove_file(observed)
        filter_motif_SNPs(fasta, SNPs, motifs, observed, complement = True)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_generate_pseudo_monomorphic_ptcs(self):
        fasta = "test_data/snp_ops/test_generate_pseudo_monomorphic_ptcs/input_nt_indicies.fasta"
        ptc_file = "test_data/snp_ops/test_generate_pseudo_monomorphic_ptcs/input_ptc_file.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_monomorphic_ptcs/observed_monomorphoc_ptc_simulants.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_monomorphic_ptcs/expected_monomorphoc_ptc_simulants.bed"
        gen.remove_file(observed)
        generate_pesudo_monomorphic_ptcs(ptc_file, fasta, observed, seed=5)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_generate_pseudo_ptc_snps(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps/expected_pseudo_ptc_snps.bed"
        expected_remaining = "test_data/snp_ops/test_generate_pseudo_ptc_snps/expected_remaining_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps/observed_pseudo_ptc_snps.bed"
        observed_remaining = "test_data/snp_ops/test_generate_pseudo_ptc_snps/observed_remaining_snps.bed"
        gen.remove_file(observed)
        gen.remove_file(observed_remaining)
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, observed_remaining, seed=5)
        expected = gen.read_many_fields(expected, "\t")
        expected_remaining = gen.read_many_fields(expected_remaining, "\t")
        observed = gen.read_many_fields(observed, "\t")
        observed_remaining = gen.read_many_fields(observed_remaining, "\t")
        self.assertEqual(observed, expected)
        self.assertEqual(observed_remaining, expected_remaining)

    def test_generate_pseudo_ptc_snps_group_by_gene(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/expected_pseudo_ptc_snps.bed"
        expected_remainder = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/expected_remaining_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/observed_pseudo_ptc_snps.bed"
        observed_remainder = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/observed_remaining_snps.bed"
        gen.remove_file(observed)
        gen.remove_file(observed_remainder)
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, observed_remainder, group_by_gene=True, seed=10)
        expected = gen.read_many_fields(expected, "\t")
        expected_remainder = gen.read_many_fields(expected_remainder, "\t")
        observed = gen.read_many_fields(observed, "\t")
        observed_remainder = gen.read_many_fields(observed_remainder, "\t")
        self.assertEqual(observed, expected)
        self.assertEqual(observed_remainder, expected_remainder)

    def test_generate_pseudo_ptc_snps_match_allele_frequency(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/input_nonsyn_snps.bed"
        expected1 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/expected_pseudo_ptc_snps1.bed"
        expected_remaining1 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/expected_remaining_snps1.bed"
        expected2 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/expected_pseudo_ptc_snps2.bed"
        expected_remaining2 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/expected_remaining_snps2.bed"
        observed1 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/observed_pseudo_ptc_snps1.bed"
        observed_remaining1 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/observed_remaining_snps1.bed"
        observed2 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/observed_pseudo_ptc_snps2.bed"
        observed_remaining2 = "test_data/snp_ops/test_generate_pseudo_ptc_snps_match_allele_frequency/observed_remaining_snps2.bed"
        gen.remove_file(observed1)
        gen.remove_file(observed_remaining1)
        gen.remove_file(observed2)
        gen.remove_file(observed_remaining2)
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed1, observed_remaining1, match_allele_frequency=True, match_allele_frequency_window=0.05, seed=3)
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed2, observed_remaining2, match_allele_frequency=True, match_allele_frequency_window=0.075, seed=3)
        expected1 = gen.read_many_fields(expected1, "\t")
        expected_remaining1 = gen.read_many_fields(expected_remaining1, "\t")
        expected2 = gen.read_many_fields(expected2, "\t")
        expected_remaining2 = gen.read_many_fields(expected_remaining2, "\t")
        observed1 = gen.read_many_fields(observed1, "\t")
        observed_remaining1 = gen.read_many_fields(observed_remaining1, "\t")
        observed2 = gen.read_many_fields(observed2, "\t")
        observed_remaining2 = gen.read_many_fields(observed_remaining2, "\t")
        self.assertEqual(observed1, expected1)
        self.assertEqual(observed_remaining1, expected_remaining1)
        self.assertEqual(observed2, expected2)
        self.assertEqual(observed_remaining2, expected_remaining2)

    def test_generate_pseudo_ptc_snps_without_replacement(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/expected_pseudo_ptc_snps.bed"
        expected_remaining = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/expected_remaining_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/observed_pseudo_ptc_snps.bed"
        observed_remaining = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/observed_remaining_snps.bed"
        gen.remove_file(observed)
        gen.remove_file(observed_remaining)
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, observed_remaining, seed=5, without_replacement=True)
        expected = gen.read_many_fields(expected, "\t")
        expected_remaining = gen.read_many_fields(expected_remaining, "\t")
        observed = gen.read_many_fields(observed, "\t")
        observed_remaining = gen.read_many_fields(observed_remaining, "\t")
        self.assertEqual(observed, expected)
        self.assertEqual(observed_remaining, expected_remaining)

    def test_get_allele_frequency(self):
        input_snps = "test_data/snp_ops/test_get_allele_frequency/input_snps.bed"
        input_snps = gen.read_many_fields(input_snps, "\t")
        expected = [np.divide(2, 6), np.divide(3, 6), np.divide(4, 6), np.divide(3, 6), np.divide(24, 5008), np.divide(15, 5008)]
        observed = []
        for snp in input_snps[1:]:
            observed.append(get_allele_frequency(snp))
        self.assertEqual(observed, expected)

    def test_get_codon_start(self):
        cases = [(10, 14), (11, 14), (9, 14), (0, 14), (13, 14)]
        expected = [9, 9, 9, 0, 12]
        observed = [get_codon_start(i[0], i[1]) for i in cases]
        self.assertEqual(expected, observed)

    def test_get_codon_start_shift(self):
        cases = [(9, 15), (0, 15), (4, 15), (8, 15), (13, 15)]
        expected = [8, "error", 4, 7, "error"]
        observed = [get_codon_start(i[0], i[1], shift = True) for i in cases]
        self.assertEqual(expected, observed)

    def test_get_snp_relative_cds_position(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position/test_snp_relative_exon_position.bed", "\t")
        full_bed_file = "test_data/snp_ops/test_get_snp_relative_cds_position/test_full_bed.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, full_bed_file)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_get_snp_relative_cds_position_minus_strand(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand/test_snp_relative_exon_position.bed", "\t")
        bed_file = "test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand/full_bed.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, bed_file)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_get_snp_relative_cds_position_plus_strand(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand/test_snp_relative_exon_position.bed", "\t")
        bed_file = "test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand/full_bed.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, bed_file)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_get_snp_relative_cds_position_minus_strand_split(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand_split/test_snp_relative_exon_position.bed", "\t")
        bed_file = "test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand_split/full_bed.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand_split/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position_minus_strand_split/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, bed_file)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_get_snp_relative_cds_position_plus_strand_split(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand_split/test_snp_relative_exon_position.bed", "\t")
        bed_file = "test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand_split/full_bed.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand_split/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position_plus_strand_split/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, bed_file)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_get_snp_relative_exon_position(self):
        intersect_file = "test_data/snp_ops/test_get_snp_relative_exon_position/test_cds_bed_snps_vcf_intersect.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_exon_position/expected_test_snp_relative_exon_position.bed", "\t")
        observed = get_snp_relative_exon_position(intersect_file)
        self.assertEqual(observed, expected)

    def test_get_snp_relative_exon_position_minus_strand(self):
        intersect_file = "test_data/snp_ops/test_get_snp_relative_exon_position_minus_strand/test_cds_bed_snps_vcf_intersect.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_exon_position_minus_strand/expected_test_snp_relative_exon_position.bed", "\t")
        observed = get_snp_relative_exon_position(intersect_file)
        self.assertEqual(observed, expected)

    def test_get_snp_change_status(self):
        snp_relative_cds_position = "test_data/snp_ops/test_get_snp_change_status/test_snp_relative_cds_position.bed"
        cds_fasta = "test_data/snp_ops/test_get_snp_change_status/test_get_snp_change_status_cds_fasta.fasta"
        expected_ptc_snps = gen.read_many_fields("test_data/snp_ops/test_get_snp_change_status/expected_get_snp_ptcs.bed", "\t")
        expected_other_snps = gen.read_many_fields("test_data/snp_ops/test_get_snp_change_status/expected_get_snp_others.bed", "\t")
        observed_ptc_snps = "test_data/snp_ops/test_get_snp_change_status/observed_get_snp_ptc_status.bed"
        observed_other_snps = "test_data/snp_ops/test_get_snp_change_status/observed_get_snp_other_status.bed"
        gen.remove_file(observed_ptc_snps)
        gen.remove_file(observed_other_snps)
        get_snp_change_status(snp_relative_cds_position, cds_fasta, observed_ptc_snps, observed_other_snps)
        observed_ptc_snps = gen.read_many_fields(observed_ptc_snps, "\t")
        observed_other_snps = gen.read_many_fields(observed_other_snps, "\t")
        self.assertEqual(observed_ptc_snps, expected_ptc_snps)
        self.assertEqual(observed_other_snps, expected_other_snps)

    def test_get_snp_type(self):
        cds_list = gen.read_many_fields("test_data/snp_ops/test_get_snp_type/test_cdss.bed", "\t")
        snp_info = gen.read_many_fields("test_data/snp_ops/test_get_snp_type/test_snp_cds_info.bed", "\t")
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_type/expected_snp_types.bed", "\t")
        observed = []
        for i, snp in enumerate(snp_info):
            cds_codon, snp_codon, mutation_type = get_snp_type(cds_list[i][0], snp)
            observed.append([cds_codon, snp_codon, mutation_type])
        self.assertEqual(observed, expected)

    def test_get_snps_in_cds(self):
        bed = "test_data/snp_ops/test_get_snps_in_cds/exons.bed"
        vcf_folder = "test_data/snp_ops/test_get_snps_in_cds/per_sample_vcfs"
        names = ["HG1", "HG3"]
        sample_file = "test_data/snp_ops/test_get_snps_in_cds/sample_file.txt"
        panel_file = "test_data/snp_ops/test_get_snps_in_cds/panel_file.txt"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snps_in_cds/expected.bed", "\t")
        observed = "test_data/snp_ops/test_get_snps_in_cds/observed.bed"
        gen.remove_file(observed)
        gen.remove_file(sample_file)
        get_snps_in_cds(bed, bed, vcf_folder, panel_file, names, sample_file, observed, out_prefix = "test_data/snp_ops/test_get_snps_in_cds/test")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_merge_and_header(self):
        file1 = "test_data/snp_ops/test_merge_and_header/file1.txt"
        file2 = "test_data/snp_ops/test_merge_and_header/file2.txt"
        expected = "test_data/snp_ops/test_merge_and_header/expected.txt"
        observed = "test_data/snp_ops/test_merge_and_header/observed.txt"
        gen.remove_file(observed)
        merge_and_header(file1, file2, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_tabix(self):
        bed_file = "test_data/snp_ops/test_tabix/test_tabix_bed.txt"
        expected = gen.read_many_fields("test_data/snp_ops/test_tabix/expected_test_tabix.txt", "\t")
        observed = "test_data/snp_ops/observed_test_tabix.bed"
        gen.remove_file(observed)
        vcf = "../source_data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.gz"
        tabix(bed_file, observed, vcf)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(sorted(observed), sorted(expected))

    def test_tabix_samples(self):
        bed_file = "test_data/snp_ops/test_tabix_samples/test_tabix.bed"
        with open("test_data/snp_ops/test_tabix_samples/expected_test_tabix_samples.txt") as file:
            expected = "".join(file)
        observed = "test_data/snp_ops/test_tabix_samples/observed_test_tabix_samples.txt"
        gen.remove_file(observed)
        gen.remove_file(observed + ".gz")
        vcf_folder = "../source_data/per_sample_vcfs"
        panel_file = "../source_data/integrated_call_samples_v3.20130502.ALL.panel"
        tabix_samples(bed_file, observed + ".gz", panel_file, vcf_folder, superpop = "AFR")
        gen.run_process(["bgzip", "-d", observed + ".gz"])
        with open(observed) as file:
            observed = "".join(file)
        expected = re.sub("0\.[0-9]*\.vcf", "N.vcf", expected)
        observed = re.sub("0\.[0-9]*\.vcf", "N.vcf", observed)
        expected = re.sub("source_[0-9]*\.[0-9]*", "source_N", expected)
        observed = re.sub("source_[0-9]*\.[0-9]*", "source_N", observed)
        self.assertEqual(observed, expected)

    def test_tabix_samples2(self):
        bed_file = "test_data/snp_ops/test_tabix_samples2/test_tabix.bed"
        with open("test_data/snp_ops/test_tabix_samples2/expected_test_tabix_samples2.txt") as file:
            expected = "".join(file)
        observed = "test_data/snp_ops/test_tabix_samples2/observed_test_tabix_samples2.txt"
        gen.remove_file(observed)
        gen.remove_file(observed + ".gz")
        vcf_folder = "../source_data/per_sample_vcfs"
        panel_file = "../source_data/integrated_call_samples_v3.20130502.ALL.panel"
        tabix_samples(bed_file, observed + ".gz", panel_file, vcf_folder, samples = ["NA18917", "NA19024"])
        gen.run_process(["bgzip", "-d", observed + ".gz"])
        with open(observed) as file:
            observed = "".join(file)
        expected = re.sub("0\.[0-9]*\.vcf", "N.vcf", expected)
        observed = re.sub("0\.[0-9]*\.vcf", "N.vcf", observed)
        expected = re.sub("source_[0-9]*\.[0-9]*", "source_N", expected)
        observed = re.sub("source_[0-9]*\.[0-9]*", "source_N", observed)
        self.assertEqual(observed, expected)
