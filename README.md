Collaborative project between Liam Abrahams and Rosina Savisaar.

DISCLAIMER
The code in this repository has not been packaged to be run out of the box. For instance, there are many software dependencies that haven't been explicitly documented, and the documentation regarding how to run the scripts is also insufficiently documented for external use. The purpose of the repository is to serve as an explicit record of the analyses performed in the paper. It is thus primarily a supplement to the methods. Some functions in some files may not be used.
****************

Source code for the Abrahams, Savisaar and Hurst manuscript in preparation Sept 2018.

motif_sets/: contains lists of motifs used in analyses.

test_data/: contains test cases for functions.

bam_ops.py: contains functions for operations on bam files, or data derived from bam files.

bed_ops.py: contains functions for operations on sequence coordinates.

compare_expression.py: compare expression parameters for transcripts that contain true PTCs vs transcripts that contain pseudo-PTCs.

disease_snps_ops.py: contains functions for operations related to disease snps files and data.

disease_snps.py: analyses on nonsense mutations found in the ClinVar database.

disease_snps.R: subsequent analyses functions on the disease data.

ese_stops_simulation_prelim.R:

ese_stops_simulation.py: analyses on stop codons found in motif sets.

generic.py: module containing generic utility functions.

map_FANTOM.py: mapping FANTOM5 CAGE peaks to gene promoters and calculating a series of expression parameters for each gene.

NAS_analysis.py: main scripts to check whether PTCs are associated with greater rates of exon skipping.

NAS_misc_tests.py: miscellaneous tests on large effect cases.

NAS_ops.py: pperations use in NAS_analysis.

NAS_prelim.R: main R analyses.

prepare_FANTOM.py: filtering a FANTOM5 osc file and writing data to BED format for CrossMapping.

simulate_eses_ops.py: operations used to simulate eses.

SNP_ops.py: functions for processing SNP data and SNP-related operations.

tests_bam_ops.py: tests for bam_ops functions.

test_bed_ops.py: tests for bed_ops functions.

tests_disease_snps_ops.py: tests for disease_snps_ops functions.

tests_generic.py: tests for generic functions.

tests_simulate_eses_ops.py: tests for simulate_eses_ops functions.

tests_SNP_ops.py: tests for SNP_ops functions.
