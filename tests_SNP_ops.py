import generic as gen
import re
from SNP_ops import *
import unittest

class Test_SNP_ops(unittest.TestCase):

    def test_filter_by_snp_type(self):
        input_snps = "test_data/snp_ops/test_filter_by_snp_type/input_snps.bed"
        expected = "test_data/snp_ops/test_filter_by_snp_type/expected_snps.bed"
        observed = "test_data/snp_ops/test_filter_by_snp_type/observed_snps.bed"
        filter_by_snp_type(input_snps, observed, "non")
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_generate_pseudo_ptc_snps_group_by_gene(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/expected_pseudo_ptc_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps_group_by_gene/observed_pseudo_ptc_snps.bed"
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, group_by_gene=True, seed=10)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_generate_pseudo_ptc_snps(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps/expected_pseudo_ptc_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps/observed_pseudo_ptc_snps.bed"
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, seed=5)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

    def test_generate_pseudo_ptc_snps_without_replacement(self):
        input_ptc_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/input_ptc_snps.bed"
        input_nonsyn_snps = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/input_nonsyn_snps.bed"
        expected = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/expected_pseudo_ptc_snps.bed"
        observed = "test_data/snp_ops/test_generate_pseudo_ptc_snps_without_replacement/observed_pseudo_ptc_snps.bed"
        generate_pseudo_ptc_snps(input_ptc_snps, input_nonsyn_snps, observed, seed=5, without_replacement=True)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed, expected)

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
