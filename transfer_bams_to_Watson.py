'''
Author: Rosina Savisaar
Wraps bam_ops.retrieve_bams.
See bam_ops.retrieve_bams for details on variables.
Note that if you don't want to specify a subset of files to download, you should put in
"all" as the last argument.
EX: python3 transfer_bams_to_Watson.py ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/processed temp_data rs949@bssv-watson:../../../../mnt/data/lars1/bams password.txt 2" 
'''

import bam_ops as bmo
import generic as gen

def main():

    description = "Retrieve .bam files from an FTP site, transfer each one to a remote server and then delete it locally."
    args = gen.parse_arguments(description, ["ftp_site", "local_directory", "remote_directory", "password_file", "subset"])
    ftp_site, local_directory, remote_directory, password_file, subset = args.ftp_site, args.local_directory, args.remote_directory, args.password_file, args.subset

    if subset == "all":
        subset = None
    else:
        subset = int(subset)
    bmo.retrieve_bams(ftp_site, local_directory, remote_directory, password_file, subset = subset)

if __name__ == "__main__":
    main()
