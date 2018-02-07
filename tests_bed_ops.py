from bed_ops import *
import generic as gen
import unittest

class Test_bed_ops(unittest.TestCase):

    def test_extract_cds(self):
        gtf_file = "./test_data/bed_ops/test_extract_cds/test_extract_cds.gtf"
        genome_file = "./test_data/bed_ops/test_extract_cds/test_extract_cds_genome.fa"
        observed = "./test_data/bed_ops/test_extract_cds/observed_test_extract_cds_fasta.fasta"
        gen.remove_file(observed)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds/expected_test_extract_cds_fasta.fasta")
        extract_cds(gtf_file, observed, genome_file, random_directory=True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_quality_control(self):
        gtf_file = "./test_data/bed_ops/test_extract_cds_quality_control/test_extract_cds_quality_control.gtf"
        genome_file = "./test_data/bed_ops/test_extract_cds_quality_control/test_extract_cds_quality_control_genome.fa"
        observed = "./test_data/bed_ops/test_extract_cds_quality_control/observed_test_extract_cds_quality_control_fasta.fasta"
        gen.remove_file(observed)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_quality_control/expected_test_extract_cds_quality_control_fasta.fasta")
        extract_cds(gtf_file, observed, genome_file, random_directory=True, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_from_bed(self):
        bed_file = "./test_data/bed_ops/test_extract_cds_from_bed/test_extract_cds_from_bed_bed.bed"
        observed = "./test_data/bed_ops/test_extract_cds_from_bed/observed_test_extract_cds_from_bed_fasta.fasta"
        gen.remove_file(observed)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_from_bed/expected_test_extract_cds_from_bed_fasta.fasta")
        extract_cds_from_bed(bed_file, observed, "./test_data/bed_ops/test_extract_cds_from_bed/test_extract_cds_from_bed_genome.fa", random_directory=True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_from_bed_quality_control(self):
        bed_file = "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/test_extract_cds_from_bed_quality_control_bed.bed"
        observed = "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/observed_test_extract_cds_from_bed_quality_control_fasta.fasta"
        gen.remove_file(observed)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_from_bed_quality_control/expected_test_extract_cds_from_bed_quality_control_fasta.fasta")
        extract_cds_from_bed(bed_file, observed, "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/test_extract_cds_from_bed_quality_control_genome.fa", random_directory=True, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_exons(self):
        gtf = "test_data/bed_ops/test_extract_exons/test_extract_exons.gtf"
        observed = "test_data/bed_ops/test_extract_exons/observed_test_extract_exons.bed"
        gen.remove_file(observed)
        extract_exons(gtf, observed)
        expected = gen.read_many_fields("test_data/bed_ops/test_extract_exons/expected_test_extract_exons.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_extract_exon_junctions(self):
        exons = "test_data/bed_ops/test_extract_exon_junctions/test_extract_exon_junctions.bed"
        observed = "test_data/bed_ops/test_extract_exon_junctions/observed_test_extract_exon_junctions.bed"
        gen.remove_file(observed)
        extract_exon_junctions(exons, observed)
        expected = gen.read_many_fields("test_data/bed_ops/test_extract_exon_junctions/expected_test_extract_exon_junctions.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_extract_exon_junctions_window(self):
        exons = "test_data/bed_ops/test_extract_exon_junctions_window/test_extract_exon_junctions.bed"
        observed = "test_data/bed_ops/test_extract_exon_junctions_window/observed_test_extract_exon_window_junctions.bed"
        gen.remove_file(observed)
        extract_exon_junctions(exons, observed, 30)
        expected = gen.read_many_fields("test_data/bed_ops/test_extract_exon_junctions_window/expected_test_extract_exon_window_junctions.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_extract_features_cdss(self):
        gtf_file = "test_data/bed_ops/test_extract_features_cdss/test_extract_features.gtf"
        observed = "test_data/bed_ops/test_extract_features_cdss/observed_test_extract_features_cdss.bed"
        gen.remove_file(observed)
        extract_features(gtf_file, observed, ['CDS'])
        expected = gen.read_many_fields("test_data/bed_ops/test_extract_features_cdss/expected_test_extract_features_cdss.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed,expected)

    def test_extract_features_cdss_stops(self):
        gtf_file = "test_data/bed_ops/test_extract_features_cdss_stops/test_extract_features.gtf"
        observed = "test_data/bed_ops/test_extract_features_cdss_stops/observed_test_extract_features_cdss_stops.bed"
        gen.remove_file(observed)
        extract_features(gtf_file, observed, ['CDS', 'stop_codon'])
        expected = gen.read_many_fields("test_data/bed_ops/test_extract_features_cdss_stops/expected_test_extract_features_cdss_stops.bed", "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(observed,expected)

    def test_fasta_from_intervals(self):
        observed = "test_data/bed_ops/test_fasta_from_intervals/observed_converted_fasta.fasta"
        gen.remove_file(observed)
        bed_file = "test_data/bed_ops/test_fasta_from_intervals/test_bed_for_fasta_conversion.bed"
        expected = gen.read_fasta("test_data/bed_ops/test_fasta_from_intervals/expected_converted_fasta.fasta")
        fasta_from_intervals(bed_file, observed, "./source_data/Genomes/test_genome/test_genome.fa")
        observed = gen.read_fasta(observed)
        self.assertEqual(observed,expected)

    def test_fasta_sequence_quality_control(self):
        fasta_parts = gen.read_fasta('./test_data/bed_ops/test_fasta_sequence_quality_control/test_fasta_sequence_quality_control.fasta')
        fasta_parts_names = fasta_parts[0]
        fasta_parts_seqs = fasta_parts[1]
        names, seqs = check_sequence_quality(fasta_parts_names, fasta_parts_seqs, check_acgt=True, check_stop=True, check_start=True, check_length=True, check_inframe_stop=True)
        observed = (names, seqs)
        expected = gen.read_fasta('./test_data/bed_ops/test_fasta_sequence_quality_control/expected_test_fasta_sequence_quality_control.fasta')
        self.assertEqual(observed, expected)

    def test_filter_bed_from_fasta(self):
        bed = "test_data/bed_ops/test_filter_bed_from_fasta/test_filter_bed_from_fasta.bed"
        fasta = "test_data/bed_ops/test_filter_bed_from_fasta/test_filter_bed_from_fasta.fasta"
        observed = "test_data/bed_ops/test_filter_bed_from_fasta/observed_test_filter_bed_from_fasta.bed"
        gen.remove_file(observed)
        expected = "test_data/bed_ops/test_filter_bed_from_fasta/expected_test_filter_bed_from_fasta.bed"
        filter_bed_from_fasta(bed, fasta, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)
