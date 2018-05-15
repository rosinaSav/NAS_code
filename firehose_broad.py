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


def process_filelist(filelist, output_dir, data_type):
    '''
    Download and process the files to the output directory.
    filelist: the list of files to be downloaded
    output_dir: the directory for the files to end up in
    data_type: the data type being downloaded (simply for printing out)
    '''

    for i, file in enumerate(filelist):

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
                args = ["cp", filepath, output_dir]
                gen.run_process(args)

        # remove the temp dir
        shutil.rmtree(temp_dir)



def main():

    description = "Download disease data."
    arguments = ["output_directory", "clean_download", "change_env"]
    args = gen.parse_arguments(description, arguments, flags = [1, 2])
    output_directory, clean_download, change_env = args.output_directory, args.clean_download, args.change_env

    # delete the directories if wanting a fresh download
    if clean_download:
        print("Deleting directories for fresh download...")
        gen.remove_directory(output_directory)

    mutation_dir = "{0}/mutations".format(output_directory)
    rna_seq_dir = "{0}/rna_seq".format(output_directory)

    gen.create_output_directories(mutation_dir)
    gen.create_output_directories(rna_seq_dir)
    gen.create_output_directories("temp_data")

    if change_env:
        args = ["source", "activate", "py27"]
        gen.run_process(args)
    # get the data dump file from firehose broad
    print("Getting data from Firehose Broad...")
    json_file = "./data.json"
    args = ["python2", "firehose_broad_get.py", json_file]
    gen.run_process(args)

    if change_env:
        args = ["source", "activate", "root"]
        gen.run_process(args)

    # read the data file
    with open(json_file, "r") as json_open:
        data = json.load(json_open)

    # extract the required file urls
    print("Getting files...")
    mrna_files = extract_files(data, "mRNASeq", "junction_quantification__data.Level_3")
    mutation_files = extract_files(data, "MAF", "Mutation_Packager_Calls.Level_3")

    # download files
    print("Downloading files...")
    gen.run_in_parallel(mutation_files, ["foo", mutation_dir, "mutation_files"], process_filelist)
    gen.run_in_parallel(mrna_files, ["foo", rna_seq_dir, "mRNASeq"], process_filelist)

    # delete the json file
    gen.remove_file(json_file)

if __name__ == "__main__":
    main()
