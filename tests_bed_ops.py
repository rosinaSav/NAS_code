from bed_ops import *
import generic as gen
import unittest

class Test_bed_ops(unittest.TestCase):

    def test_change_bed_names(self):
        full_bed = "test_data/bed_ops/test_change_bed_names/input_full.bed"
        short_bed = "test_data/bed_ops/test_change_bed_names/input_short.bed"
        expected_full = "test_data/bed_ops/test_change_bed_names/expected_full.bed"
        expected_short = "test_data/bed_ops/test_change_bed_names/expected_short.bed"
        observed_full = "test_data/bed_ops/test_change_bed_names/observed_full.bed"
        observed_short = "test_data/bed_ops/test_change_bed_names/observed_short.bed"
        gen.remove_file(observed_full)
        gen.remove_file(observed_short)
        change_bed_names(full_bed, observed_full, full_names=True, header=False)
        change_bed_names(short_bed, observed_short, full_names=False, header=False)
        expected_full = gen.read_many_fields(expected_full, "\t")
        expected_short = gen.read_many_fields(expected_short, "\t")
        observed_full = gen.read_many_fields(observed_full, "\t")
        observed_short = gen.read_many_fields(observed_short, "\t")
        self.assertEqual(expected_full, observed_full)
        self.assertEqual(expected_short, observed_short)

    def test_check_coding(self):
        exon_file = "test_data/bed_ops/test_check_coding/exons.bed"
        CDS_file = "test_data/bed_ops/test_check_coding/CDSs.bed"
        expected = "test_data/bed_ops/test_check_coding/expected_check_coding.bed"
        observed = "test_data/bed_ops/test_check_coding/observed_check_coding.bed"
        gen.remove_file(observed)
        check_coding(exon_file, CDS_file, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_extract_cds(self):
        gtf_file = "./test_data/bed_ops/test_extract_cds/test_extract_cds.gtf"
        genome_file = "./test_data/bed_ops/test_extract_cds/test_extract_cds_genome.fa"
        bed_output = "./test_data/bed_ops/test_extract_cds/test_extract_cds.bed"
        observed = "./test_data/bed_ops/test_extract_cds/observed_test_extract_cds_fasta.fasta"
        gen.remove_file(observed)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds/expected_test_extract_cds_fasta.fasta")
        extract_cds(gtf_file, bed_output, observed, genome_file, full_chr_name = True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_clean_chrom(self):
        gtf_file = "./test_data/bed_ops/test_extract_cds_clean_chrom/test_extract_cds.gtf"
        genome_file = "./test_data/bed_ops/test_extract_cds_clean_chrom/test_extract_cds_genome.fa"
        bed_output = "./test_data/bed_ops/test_extract_cds_clean_chrom/test_extract_cds.bed"
        observed = "./test_data/bed_ops/test_extract_cds_clean_chrom/observed_test_extract_cds_fasta.fasta"
        intervals = "./test_data/bed_ops/test_extract_cds_clean_chrom/observed_intervals.fasta"
        gen.remove_file(observed)
        gen.remove_file(intervals)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_clean_chrom/expected_test_extract_cds_fasta.fasta")
        extract_cds(gtf_file, bed_output, observed, genome_file, intervals, clean_chrom_only = True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_quality_control(self):
        gtf_file = "./test_data/bed_ops/test_extract_cds_quality_control/test_extract_cds_quality_control.gtf"
        genome_file = "./test_data/bed_ops/test_extract_cds_quality_control/test_extract_cds_quality_control_genome.fa"
        bed_output = "./test_data/bed_ops/test_extract_cds_quality_control/observed_test_extract_cds_quality_control_cds.bed"
        observed = "./test_data/bed_ops/test_extract_cds_quality_control/observed_test_extract_cds_quality_control_fasta.fasta"
        intervals = "./test_data/bed_ops/test_extract_cds_quality_control/observed_intervals.fasta"
        gen.remove_file(observed)
        gen.remove_file(intervals)
        gen.remove_file(bed_output)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_quality_control/expected_test_extract_cds_quality_control_fasta.fasta")
        extract_cds(gtf_file, bed_output, observed, genome_file, intervals, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_from_bed(self):
        bed_file = "./test_data/bed_ops/test_extract_cds_from_bed/test_extract_cds_from_bed_bed.bed"
        observed = "./test_data/bed_ops/test_extract_cds_from_bed/observed_test_extract_cds_from_bed_fasta.fasta"
        intervals = "./test_data/bed_ops/test_extract_cds_from_bed/observed_intervals.fasta"
        gen.remove_file(observed)
        gen.remove_file(intervals)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_from_bed/expected_test_extract_cds_from_bed_fasta.fasta")
        extract_cds_from_bed(bed_file, observed, "./test_data/bed_ops/test_extract_cds_from_bed/test_extract_cds_from_bed_genome.fa", intervals)
        observed = gen.read_fasta(observed)
        self.assertEqual(observed, expected)

    def test_extract_cds_from_bed_quality_control(self):
        bed_file = "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/test_extract_cds_from_bed_quality_control_bed.bed"
        observed = "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/observed_test_extract_cds_from_bed_quality_control_fasta.fasta"
        intervals = "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/observed_intervals.fasta"
        gen.remove_file(observed)
        gen.remove_file(intervals)
        expected = gen.read_fasta("./test_data/bed_ops/test_extract_cds_from_bed_quality_control/expected_test_extract_cds_from_bed_quality_control_fasta.fasta")
        extract_cds_from_bed(bed_file, observed, "./test_data/bed_ops/test_extract_cds_from_bed_quality_control/test_extract_cds_from_bed_quality_control_genome.fa", intervals, check_acgt=True, check_start=True, check_length=True, check_stop=True, check_inframe_stop=True)
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

    def test_extract_nt_indices(self):
        fasta_file = "test_data/bed_ops/test_extract_nt_indices/test_fasta_input.fasta"
        observed_a = "test_data/bed_ops/test_extract_nt_indices/observed_indices_a.fasta"
        observed_c = "test_data/bed_ops/test_extract_nt_indices/observed_indices_c.fasta"
        observed_g = "test_data/bed_ops/test_extract_nt_indices/observed_indices_g.fasta"
        observed_t = "test_data/bed_ops/test_extract_nt_indices/observed_indices_t.fasta"
        expected_a = "test_data/bed_ops/test_extract_nt_indices/expected_indices_a.fasta"
        expected_c = "test_data/bed_ops/test_extract_nt_indices/expected_indices_c.fasta"
        expected_g = "test_data/bed_ops/test_extract_nt_indices/expected_indices_g.fasta"
        expected_t = "test_data/bed_ops/test_extract_nt_indices/expected_indices_t.fasta"
        gen.remove_file(observed_a)
        gen.remove_file(observed_c)
        gen.remove_file(observed_g)
        gen.remove_file(observed_t)
        observed_files = {
            "A": observed_a,
            "C": observed_c,
            "G": observed_g,
            "T": observed_t,
        }
        extract_nt_indices(fasta_file, observed_files)
        expected_a = gen.read_fasta(expected_a)
        expected_c = gen.read_fasta(expected_c)
        expected_g = gen.read_fasta(expected_g)
        expected_t = gen.read_fasta(expected_t)
        observed_a = gen.read_fasta(observed_a)
        observed_c = gen.read_fasta(observed_c)
        observed_g = gen.read_fasta(observed_g)
        observed_t = gen.read_fasta(observed_t)
        self.assertEqual(observed_a,expected_a)
        self.assertEqual(observed_c,expected_c)
        self.assertEqual(observed_g,expected_g)
        self.assertEqual(observed_t,expected_t)

    def test_fasta_from_intervals(self):
        observed = "test_data/bed_ops/test_fasta_from_intervals/observed_converted_fasta.fasta"
        gen.remove_file(observed)
        bed_file = "test_data/bed_ops/test_fasta_from_intervals/test_bed_for_fasta_conversion.bed"
        expected = gen.read_fasta("test_data/bed_ops/test_fasta_from_intervals/expected_converted_fasta.fasta")
        fasta_from_intervals(bed_file, observed, "test_data/bed_ops/test_fasta_from_intervals/test_genome.fa")
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

    def test_filter_bed_from_fasta_families(self):
        bed = "test_data/bed_ops/test_filter_bed_from_fasta_families/test_filter_bed_from_fasta.bed"
        fasta = "test_data/bed_ops/test_filter_bed_from_fasta_families/test_filter_bed_from_fasta.fasta"
        families_file = "test_data/bed_ops/test_filter_bed_from_fasta_families/families.txt"
        observed = "test_data/bed_ops/test_filter_bed_from_fasta_families/observed_test_filter_bed_from_fasta.bed"
        gen.remove_file(observed)
        expected = "test_data/bed_ops/test_filter_bed_from_fasta_families/expected.txt"
        filter_bed_from_fasta(bed, fasta, observed, families_file = families_file)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_filter_fasta_intervals_from_fasta(self):
        fasta = "test_data/bed_ops/test_filter_fasta_intervals_from_fasta/test_filter_fasta_intervals_from_fasta_fasta.fasta"
        fasta_intervals = "test_data/bed_ops/test_filter_fasta_intervals_from_fasta/test_filter_fasta_intervals_from_fasta_intervals.fasta"
        observed = "test_data/bed_ops/test_filter_fasta_intervals_from_fasta/observed_filtered_intervals.fasta"
        gen.remove_file(observed)
        expected = "test_data/bed_ops/test_filter_fasta_intervals_from_fasta/expected_filtered_intervals.fasta"
        filter_fasta_intervals_from_fasta(fasta_intervals, fasta, observed)
        expected = gen.read_fasta(expected)
        observed = gen.read_fasta(observed)
        self.assertEqual(expected, observed)

    def test_filter_exon_junctions(self):
        exon_junctions_file = "test_data/bed_ops/test_filter_exon_junctions/exon_junctions.bed"
        exons_file = "test_data/bed_ops/test_filter_exon_junctions/exons.bed"
        expected = "test_data/bed_ops/test_filter_exon_junctions/expected_filter_exon_junctions.bed"
        observed = "test_data/bed_ops/test_filter_exon_junctions/observed_filter_exon_junctions.bed"
        gen.remove_file(observed)
        filter_exon_junctions(exon_junctions_file, exons_file, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_get_descriptions(self):
        gtf = "test_data/bed_ops/test_get_descriptions/descriptions.gtf"
        names = ["ENST100", "ENST7", "ENST0003", "ENST5"]
        expected = "test_data/bed_ops/test_get_descriptions/expected_get_descriptions.txt"
        observed = "test_data/bed_ops/test_get_descriptions/observed_get_descriptions.txt"
        gen.remove_file(observed)
        get_descriptions(names, gtf, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_link_genes_and_transcripts(self):
        bed = "test_data/bed_ops/test_link_genes_and_transcripts/input.bed"
        expected = {"ENSG000001": ["ENST0003246"], "ENSG00000223972": ["ENST00032132", "ENST00032323"]}
        observed = link_genes_and_transcripts(bed)
        self.assertEqual(expected, observed)

    def test_remove_overlaps(self):
        in_bed = "test_data/bed_ops/test_remove_overlaps/in.bed"
        expected = "test_data/bed_ops/test_remove_overlaps/expected.bed"
        observed = "test_data/bed_ops/test_remove_overlaps/observed.bed"
        gen.remove_file(observed)
        remove_overlaps(in_bed, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_remove_overlaps2(self):
        in_bed = "test_data/bed_ops/test_remove_overlaps2/in.bed"
        expected = "test_data/bed_ops/test_remove_overlaps2/expected.bed"
        observed = "test_data/bed_ops/test_remove_overlaps2/observed.bed"
        gen.remove_file(observed)
        remove_overlaps(in_bed, observed)
        expected = gen.read_many_fields(expected, "\t")
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(expected, observed)

    def test_uniquify_trans(self):
        gene_to_trans = {"ENSG1": ["ENST6"], "ENSG2": ["ENST2", "ENST3"], "ENSG3": ["ENST66", "ENST8", "ENST5"]}
        names = ["ENST66", "ENST3", "ENST6", "ENST8", "ENST2", "ENST5"]
        seqs = ["ATGGA", "ATGT", "ATTGGATAGT", "ATTT", "ATGATGATGAT", "ATGGGTTTTTT"]
        expected = (["ENST6", "ENST2", "ENST5"], ["ATTGGATAGT", "ATGATGATGAT", "ATGGGTTTTTT"])
        observed = uniquify_trans(names, seqs, gene_to_trans)
        self.assertEqual(expected, observed)
