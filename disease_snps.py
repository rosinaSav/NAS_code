'''
Auhor: Liam Abrahams.
Look at disease SNPs from ClinVar.
'''

import generic as gen
import bed_ops as beo
import bam_ops as bao
import SNP_ops as so
import os

def main():

    description = "Look at disease snps."
    arguments = ["intersect_snps", "get_relative_positions", "get_snp_status"]
    args = gen.parse_arguments(description, arguments, flags = [0,1,2])
    intersect_snps, get_relative_positions, get_snp_status = args.intersect_snps, args.get_relative_positions, args.get_snp_status

    results_prefix = "./results/clean_run_2/clean_run"

    output_directory = "results/disease_snps"
    gen.create_output_directories(output_directory)

    disease_snps_file = "./source_data/clinvar_20180429.vcf.gz"
    disease_snps_index_file = "{0}.tbi".format(disease_snps_file)

    if not os.path.isfile(disease_snps_file) or not os.path.isfile(disease_snps_index_file):
        print("Please provide the required disease SNPs files...")
        raise Exception

    # intersect the coding exons with the disease snps
    exon_bed = "{0}_coding_exons.bed".format(results_prefix)
    disease_snp_intersect_file_vcf = "{0}/disease_snp_intersect.vcf".format(output_directory)
    disease_snp_intersect_file_bed = "{0}/disease_snp_intersect.bed".format(output_directory)
    if intersect_snps or not os.path.isfile(disease_snp_intersect_file_vcf) or not os.path.isfile(disease_snp_intersect_file_bed):
        print("Intersecting snps with exons")
        # so.intersect_snps_parallel(exon_bed, disease_snps_file, disease_snp_intersect_file_vcf)
        so.intersect_vcf_to_bed(exon_bed, disease_snp_intersect_file_vcf, disease_snp_intersect_file_bed, change_names = True)

    # get relative positions of the snps in cds and exons
    full_bed = "{0}_CDS.bed".format(results_prefix)
    disease_snps_relative_exon_positions = "{0}/disease_snp_relative_exon_positions.bed".format(output_directory)
    disease_snps_relative_cds_positions = "{0}/disease_snp_relative_cds_positions.bed".format(output_directory)
    if get_relative_positions or not os.path.isfile(disease_snps_relative_exon_positions) or not os.path.isfile(disease_snps_relative_cds_positions):
        print("Getting snp relative positions...")
        so.get_snp_relative_exon_position(disease_snp_intersect_file_bed, disease_snps_relative_exon_positions)
        # output to var because this is how the function was made
        relative_positions = gen.read_many_fields(disease_snps_relative_exon_positions, "\t")
        so.get_snp_relative_cds_position(relative_positions, disease_snps_relative_cds_positions, full_bed)

    # get the change status of the snps to check them
    cds_fasta = "{0}_CDS.fasta".format(results_prefix)
    disease_ptcs_file = "{0}/disease_ptcs.txt".format(output_directory)
    disease_other_file = "{0}/disease_other_snps.txt".format(output_directory)
    if get_snp_status or not os.path.isfile(disease_ptcs_file) or not os.path.isfile(disease_other_file) or not os.path.isfile(cds_fasta):
        print("Getting snp status...")
        so.get_snp_change_status(disease_snps_relative_cds_positions, cds_fasta, disease_ptcs_file, disease_other_file)



if __name__ == "__main__":
    main()
