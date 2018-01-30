import bed_ops as bo
import bam_ops as bmo
import generic as gen

def main():

    ftp_site = "ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/processed"
    target_directory = "../source_data/bams"
    exons_file = "../source_data/Homo_sapiens.GRCh37.87_exons.bed"
    exon_junctions_file = "../source_data/Homo_sapiens.GRCh37.87_exon_junctions.bed"
    bo.extract_exon_junctions(exons_file, exon_junctions_file, window_of_interest = 2)
    bmo.write_hits_at_junctions_per_sample(ftp_site, target_directory, exon_junctions_file, subset = 2)

if __name__ == "__main__":
    main()
