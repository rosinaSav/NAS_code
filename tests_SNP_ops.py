import generic as gen
import re
from SNP_ops import *
import unittest

class Test_SNP_ops(unittest.TestCase):

    def test_tabix(self):
        bed_file = "test_data/test_tabix_bed.txt"
        expected = gen.read_many_fields("test_data/test_tabix_expected.txt", "\t")
        observed = "test_data/test_tabix_observed.txt"
        gen.remove_file(observed)
        vcf = "../../general/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.gz"
        tabix(bed_file, observed, vcf)
        observed = gen.read_many_fields(observed, "\t")
        self.assertEqual(sorted(observed), sorted(expected))

    def test_tabix_samples(self):
        bed_file = "test_data/test_tabix_bed.txt"
        with open("test_data/test_tabix_samples_expected.txt") as file:
            expected = "".join(file)
        observed = "test_data/test_tabix_samples_observed.txt"
        gen.remove_file(observed)
        gen.remove_file(observed + ".gz")
        vcf_folder = "../../general/1000genomes/per_sample_vcfs"
        panel_file = "../../general/1000genomes/integrated_call_samples_v3.20130502.ALL.panel"
        tabix_samples(bed_file, observed + ".gz", panel_file, vcf_folder, superpop = "AFR", vcftools_path = "../../../Software/vcftools")
        gen.run_process(["bgzip", "-d", observed + ".gz"])
        with open(observed) as file:
            observed = "".join(file)
        expected = re.sub("0\.[0-9]*\.vcf", "N.vcf", expected)
        observed = re.sub("0\.[0-9]*\.vcf", "N.vcf", observed)
        self.assertEqual(observed, expected)

        #tabix_samples(bed_file, output_file_name, panel_file, vcf_folder, superpop = None, subpop = None, samples = None, downsample_by = None, exclude_xy = False)        
