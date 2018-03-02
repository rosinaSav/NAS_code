import unittest
import simulate_eses_ops as se
import generic as gen

class Test_simulate_eses_ops(unittest.TestCase):

    def test_generate_motifs(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        dinucleotides = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs/dinucleotides.txt", ",")
        dinucleotides = [i[0] for i in dinucleotides]
        seed_list = [5]
        expected = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs/expected_simulated_motifs.txt", ",")
        expected = [i[0] for i in expected]
        observed = se.generate_motifs([0], motif_set, dinucleotides, seed=seed_list)[0]
        self.assertEqual(expected, observed)

    def test_generate_motifs_diff_lengths(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs_diff_lengths/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        dinucleotides = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs_diff_lengths/dinucleotides.txt", ",")
        dinucleotides = [i[0] for i in dinucleotides]
        seed_list = [5]
        expected = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs_diff_lengths/expected_simulated_motifs.txt", ",")
        expected = [i[0] for i in expected]
        observed = se.generate_motifs([0], motif_set, dinucleotides, seed=seed_list)[0]
        self.assertEqual(expected, observed)

    def test_generate_motifs_sets(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_generate_motifs_sets/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        seed_list = range(3)
        expected = "test_data/simulate_eses/test_generate_motifs_sets/expected_simulated_motifs.txt"
        observed = "test_data/simulate_eses/test_generate_motifs_sets/observed_simulated_motifs.txt"
        se.generate_motifs_sets(motif_set, 3, output_file=observed, seed_list=seed_list)
        expected = gen.read_many_fields(expected, "|")
        observed = gen.read_many_fields(observed, "|")
        self.assertEqual(expected, observed)

    def test_get_dinucleotides_reg(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_get_dinucleotides_reg/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        expected = gen.read_many_fields("test_data/simulate_eses/test_get_dinucleotides_reg/expected_dinucleotides.txt", ",")
        expected = [i[0] for i in expected]
        observed = se.get_dinucleotides(motif_set)
        self.assertEqual(expected, observed)

    def test_get_dinucleotides_contact(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_get_dinucleotides_concat/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        expected = gen.read_many_fields("test_data/simulate_eses/test_get_dinucleotides_concat/expected_dinucleotides.txt", ",")
        expected = [i[0] for i in expected]
        observed = se.get_dinucleotides(motif_set, concat_motifs=True)
        self.assertEqual(expected, observed)

    def test_get_stop_codon_count(self):
        motif_set = gen.read_many_fields("test_data/simulate_eses/test_get_stop_codon_count/motif_set.txt", ",")
        motif_set = [i[0] for i in motif_set]
        expected = 8
        observed = se.get_stop_codon_count(motif_set)
        self.assertEqual(expected, observed)
