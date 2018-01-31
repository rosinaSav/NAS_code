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

    def test_extract_exon_junctions(self):
        exons = "test_data/test_extract_exon_junctions.bed"
        observed = "test_data/test_extract_exon_junctions_observed.bed"
        gen.remove_file(observed)
        extract_exon_junctions(exons, observed)
        expected = gen.read_many_fields("test_data/test_extract_exon_junctions_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_extract_exon_junctions_window(self):
        exons = "test_data/test_extract_exon_junctions.bed"
        observed = "test_data/test_extract_exon_window_junctions_observed.bed"
        gen.remove_file(observed)
        extract_exon_junctions(exons, observed, 30)
        expected = gen.read_many_fields("test_data/test_extract_exon_window_junctions_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_fasta_from_intervals(self):
        gen.remove_file("test_data/observed_fasta_from_intervals.fasta")
        bed_file = "test_data/test_bed_for_fasta_conversion.bed"
        expected = gen.read_fasta("test_data/test_converted_fasta.fasta")
        observed_file_name = "test_data/observed_converted_fasta.fasta"
        fasta_from_intervals(bed_file, observed_file_name, "source_data/Genomes/test_genome/test_genome.fa")
        observed = gen.read_fasta(observed_file_name)
        self.assertEqual(observed,expected)

    def test_extract_features_cdss(self):
        gtf_file = "test_data/test_extract_features.gtf"
        observed = "test_data/test_extract_features_cdss_observed.bed"
        gen.remove_file(observed)
        extract_features(gtf_file, observed, ['CDS'])
        expected = gen.read_many_fields("test_data/test_extract_features_cdss_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed,expected)

    def test_extract_features_cdss_stops(self):
        gtf_file = "test_data/test_extract_features.gtf"
        observed = "test_data/test_extract_features_cdss_stops_observed.bed"
        gen.remove_file(observed)
        extract_features(gtf_file, observed, ['CDS', 'stop_codon'])
        expected = gen.read_many_fields("test_data/test_extract_features_cdss_stops_expected.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed,expected)

    def test_extract_fasta_temp(self):
        # gen.remove_file("test_data/observed_fasta_from_intervals.fasta")
        bed_file = "./test_data/test_extract_cds_bed.bed"
        expected = gen.read_fasta("./test_data/test_extract_cds_temp_fasta_expected.fasta")
        observed_file_name = "./temp_files/temp_fasta_files/test_extract_cds_temp_fasta_observed.fasta"
        extract_fasta_temp(bed_file, observed_file_name, "./test_data/test_extract_cds_genome_fasta.fa")
        observed = gen.read_fasta(observed_file_name)
        self.assertEqual(observed,expected)
