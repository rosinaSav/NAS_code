'''
Authors: Liam Abrahams.
Download disease data.
'''

import generic as gen
import json
import re
import random
import os
import shutil
import time

def extract_files(data, data_type, url_search=None):
    '''
    Extract the required file urls
    data: the read json data
    data_type: the data type (http://firebrowse.org/api-docs/#!/Archives/StandardData)
    url_search: search in the url for a particular part to get a specific file
    '''
    file_list = []
    for i in data:
        for sample in data[i]:
            if sample["data_type"] == data_type:
                for url in  sample['urls']:
                    if url_search:
                        search = re.findall(url_search, url)
                        if len(search) and not url.endswith('.md5'):
                            file_list.append(url)
                    else:
                        if not url.endswith('.md5'):
                            file_list.append(url)
    return file_list


def process_filelist(filelist, remote_directory, expect_password, data_type):
    '''
    Download and process the files to the output directory.
    filelist: the list of files to be downloaded
    output_dir: the directory for the files to end up in
    data_type: the data type being downloaded (simply for printing out)
    '''

    for i, file in enumerate(filelist):
        start_time = time.time()
        print("{0}: {1}/{2}...".format(data_type, i+1, len(filelist)))

        # download file to temp directory
        temp_dir = "temp_data/{0}".format(random.random())
        gen.create_output_directories(temp_dir)
        args = ["wget", file, "-P", temp_dir]
        gen.run_process(args)

        # unzip
        download_path = "{0}/{1}".format(temp_dir, file.split('/')[-1])
        args = ["tar", "-xvzf", download_path, "-C", temp_dir]
        gen.run_process(args)

        # copy files to output_dir:
        extracted_path = download_path[:-7]
        for file in os.listdir(extracted_path):
            if file != "MANIFEST.txt":
                filepath = "{0}/{1}".format(extracted_path, file)
                expect_file = "temp_data/expect_file{0}.txt".format(random.random())
                expect_string = "#!/usr/bin/expect\nset timeout -1\nspawn rsync {0} {1}\nexpect \"la466@bssv-watson's password:\"\nsend \"{2}\\n\";\nexpect eof\nexit".format(filepath, remote_directory, expect_password)
                with open(expect_file, "w") as efile:
                    efile.write(expect_string)
                gen.run_process(["expect", expect_file])
                gen.remove_file(expect_file)
        print("Time spent: {0} minutes.\n".format(round((time.time() - start_time)/60), 3))
        shutil.rmtree(temp_dir)


def main():

    description = "Download disease data."
    arguments = ["remote_mutation_dir", "remote_rna_dir", "json_file", "password_file", "subset"]
    args = gen.parse_arguments(description, arguments, flags = [])
    remote_mutation_dir, remote_rna_dir, json_file, password_file, subset = args.remote_mutation_dir, args.remote_rna_dir, args.json_file, args.password_file, args.subset

    #get password for Watson
    with open(password_file) as file:
        expect_password = "".join(file)
        expect_password = expect_password.rstrip("\n")

    gen.create_output_directories("temp_data")

    # read the data file
    with open(json_file, "r") as json_open:
        data = json.load(json_open)


    # extract the required file urls
    print("Getting files...")
    mrna_files = extract_files(data, "mRNASeq", "junction_quantification__data.Level_3")
    mutation_files = extract_files(data, "MAF", "Mutation_Packager_Calls.Level_3")


    if subset == "all":
        mutation_files = mutation_files
        mrna_files = mrna_files
    else:
        mutation_files = mutation_files[:int(subset)]
        mrna_files = mrna_files[:int(subset)]

    # download files
    print("Downloading files...")
    gen.run_in_parallel(mutation_files, ["foo", remote_mutation_dir, expect_password, "mutation_files"], process_filelist, workers = 1)
    gen.run_in_parallel(mrna_files, ["foo", remote_rna_dir, expect_password, "mRNASeq"], process_filelist)



if __name__ == "__main__":
    main()