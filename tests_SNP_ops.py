import generic as gen
import re
from SNP_ops import *
import unittest

class Test_SNP_ops(unittest.TestCase):

#     def test_get_relative_SNPs(self):
#         CDSs = {"ENST1": [[["16", "CDS", 30397645, 30397762, "ENST1", ".", "+"], 0], [["16", "CDS", 3204450, 3204500, "ENST1", ".", "+"], 1]],
#                 "ENST2": [[["16", "CDS", 3204287, 3204295, "ENST2", ".", "-"], 0]]}
#         SNP_file_name = "test_data/snp_ops/test_SNPs_for_relative_coords.txt"
#         CDS_SNP_file_name = "test_data/snp_ops/test_SNPs_clean_relative_coords.txt"
#         seqs = ["ATTTGGTTACCCGGGTCCACTCTTCTACTGGTCGCTTACGACTTGTGCGCAGTCTATCCCACACGTGGAGAAGGTGCCCCTATATTCTTAATTGACGCGAGGCGCATAGCATGCCCAGTAGGACACGTAAGACGGATTTGAGCCATCCACATTCCTTCGCCTTAGAATTCAGGCACACTTGAATAGGGCAAGACCCCTTTATGCAATATTATTGCGTTACGTCCGTACCGTGAGTC", "ATGACAGAG"]
#         names = ["ENST1", "ENST2"]
#         genome = "hg38"
#         expected_SNPs = {"ENST1": {30: ("A", 0.000199681), 114: ("T", 0.000199681), 115: ("T", 0.000199681), 138: ("C", 0.431909)},
#                          "ENST2": {5: ("T", 0.000199681)}}
#         expected_to_remove = {"ENST1": [30, 114, 115, 135, 138],
#                               "ENST2": [5, 2]}
#         expected = (expected_SNPs, expected_to_remove)
#         observed = get_relative_SNPs(CDSs, SNP_file_name, CDS_SNP_file_name, seqs, names, genome, get_new_SNPs = True, parse_SNPs = True)
# ##
# ####chr16	30397674	30397675	ENST00000622647	100	.	rs547051052$G$A$PASS$AC=1;AF=0.000199681;AN=5008;NS=2504;DP=17709;EAS_AF=0.0000;AMR_AF=0.0000;AFR_AF=0.0008;EUR_AF=0.0000;SAS_AF=0.0000;ssID=ss1355960302;ASP;MATCHED_FWD
# ####chr16	30397758	30397759	ENST00000622647	100	.	rs571888777$C$T$PASS$AC=1;AF=0.000199681;AN=5008;NS=2504;DP=15868;EAS_AF=0.0000;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0010;ssID=ss1355960303;ASP;MATCHED_FWD
# ####chr16	30397759	30397760	ENST00000622647	100	.	rs539012229$C$T$PASS$AC=1;AF=0.000199681;AN=5008;NS=2504;DP=15866;EAS_AF=0.0000;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0010;ssID=ss1355960304;ASP;MATCHED_FWD
# ####
# ####
# ####chr16	3204466	3204467	ENST00000304646	100	.	rs371443638$CCTT$C$PASS$AC=68;AF=0.0135783;AN=5008;NS=2504;DP=23930;EAS_AF=0.0020;AMR_AF=0.0072;AFR_AF=0.0234;EUR_AF=0.0010;SAS_AF=0.0297;ssID=ss1375472352;ASP;MATCHED_FWD
# ####chr16	3204469	3204470	ENST00000304646	100	.	rs1834026$T$C$PASS$AC=2163;AF=0.431909;AN=5008;NS=2504;DP=23279;EAS_AF=0.5655;AMR_AF=0.5202;AFR_AF=0.2738;EUR_AF=0.4314;SAS_AF=0.4458;ssID=ss1355074812;ASP;MATCHED_FWD;RV
# ####
# ####
# ####chr16	3204289	3204290	ENST00000304646	100	.	rs181720197$T$A$PASS$AC=1;AF=0.000199681;AN=5008;NS=2504;DP=18775;EAS_AF=0.0000;AMR_AF=0.0000;AFR_AF=0.0008;EUR_AF=0.0000;SAS_AF=0.0000;ssID=ss1355074796;ASP;MATCHED_FWD
# ####chr16	3204292	3204293	ENST00000304646	100	.	rs142486394$G$A,C$PASS$AC=2,20;AF=0.000399361,0.00399361;AN=5008;NS=2504;DP=18501;EAS_AF=0.0000,0.0000;AMR_AF=0.0000,0.0000;AFR_AF=0.0000,0.0000;EUR_AF=0.0020,0.0000;SAS_AF=0.0000,0.0204;ssID=ss1355074797,ss1355074798;ASP;MATCHED_FWD
#         self.assertEqual(expected, observed)


    def test_tabix(self):
        bed_file = "test_data/snp_ops/test_tabix/test_tabix.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_tabix/expected_test_tabix.txt", "\t")
        observed = "test_data/snp_ops/observed_test_tabix.txt"
        gen.remove_file(observed)
        vcf = "./source_data/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.gz"
        tabix(bed_file, observed, vcf)
        observed = gen.read_many_fields(observed, "\t")
        print("\n")
        print(expected)
        print("\n")
        print(observed)
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

    def test_get_snp_relative_exon_position(self):
        intersect_file = "test_data/snp_ops/test_get_snp_relative_exon_position/test_cds_bed_snps_vcf_intersect.bed"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_exon_position/expected_test_snp_relative_exon_position.bed", "\t")
        observed = get_snp_relative_exon_position(intersect_file)
        self.assertEqual(observed, expected)

    def test_get_snp_relative_cds_position(self):
        relative_exon_position_file = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position/test_snp_relative_exon_position.bed", "\t")
        fasta_interval_file = "test_data/snp_ops/test_get_snp_relative_cds_position/test_get_snp_relative_exon_position_exon_intervals.fasta"
        expected = gen.read_many_fields("test_data/snp_ops/test_get_snp_relative_cds_position/expected_test_snp_relative_cds_position.bed", "\t")
        observed = "test_data/snp_ops/test_get_snp_relative_cds_position/observed_test_snp_relative_cds_position.bed"
        gen.remove_file(observed)
        get_snp_relative_cds_position(relative_exon_position_file, observed, fasta_interval_file)
        observed = gen.read_many_fields(observed, "\t")
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
