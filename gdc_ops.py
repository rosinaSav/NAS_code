import generic as gen
import os
import random
from copy import deepcopy


def concatenate_files(filelist, output_file):
    gen.remove_file(output_file)
    args = ["cat"]
    for file in filelist:
        args.append(file)
    gen.run_process(args, file_for_output = output_file)

class SNP_entry(object):
    def __init__(self, line, headers):
        self.line = line
        self.headers = headers

    def line_to_dict(self):
        entry_dict = {}
        for i, item in enumerate(self.line):
            entry_dict[self.headers[i]] = item
        return entry_dict

def check_entry(entry, header, required_list, non_empty_list, match_list):
    entry = SNP_entry(entry, header)
    entry_dict = entry.line_to_dict()

    # check if the columns required are present
    for requirement in required_list:
        if requirement not in entry_dict:
            return False, False, None
    # check if columns required are non empty
    for requirement in non_empty_list:
        if len(entry_dict[requirement]) == 0 or entry_dict[requirement] in [" ", "-"]:
            return False, False, None
    # check the strand exists
    if entry_dict["Strand"] not in ["+", "-"]:
        return False, False, None
    # check thats things that should be the same are
    for matches in match_list:
        hits = []
        for i in matches:
            hits.append(entry_dict[i])
        if len(list(set(hits))) > 1:
            return False, False, None
    # check the build
    if entry_dict["NCBI_Build"] not in ["hg19", "37"]:
        return False, True, None

    # return passed
    return True, True, entry_dict

def process_mutation_files(input_dir, output_dir, output_file, subset=None):
    '''
    Refactor the files so they resemble something like
    a vcf file that we can use.
    dir: the directory containing the mutation files
    output_dir: the directory to output the processed files
    filename_prefix: prefix for the processed files
    '''

    required_list = [
        "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand",
        "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID"
    ]
    non_empty_list = deepcopy(required_list)
    non_empty_list.remove("Strand")
    non_empty_list.remove("Match_Norm_Seq_Allele1")
    non_empty_list.remove("Match_Norm_Seq_Allele2")

    match_list = [["Start_Position", "End_Position"]]

    output_list = [
        "Chromosome", "Start_Position", "dbSNP_RS", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Tumor_Sample_UUID", "Match_Norm_Seq_Allele1",
        "Match_Norm_Seq_Allele2", "Matched_Norm_Sample_Barcode", "Matched_Norm_Sample_UUID",
        "Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification", "Variant_Type", "Mutation_Status"
    ]

    temp_filelist = []

    readme_file = "{0}/README_processed_mutations.txt".format(output_dir)
    with open(readme_file, "w") as readme:
        readme.write("#columns\n")
        readme.write("{0}\n".format("\t".join(output_list)))

    # create a list of the file ids and their respective filepaths
    ids = []
    filelist = {}
    for id in os.listdir(input_dir):
        if ".txt" not in id:
            ids.append(id)
            filedir = "{0}/{1}".format(input_dir, id)
            for file in os.listdir(filedir):
                filepath = "{0}/{1}".format(filedir, file)
                filelist[id] = filepath

    # if we want to subset the number of files
    if subset:
        ids = ids[:subset]

    for id in ids:
        filepath = filelist[id]
        print(filepath)

        entries = gen.read_many_fields(filepath, "\t")
        # remove any information lines
        entries = [entry for entry in entries if entry[0][0] != "#"]
        header = entries[0][:34]

        # create a temporary file to hold the mutations
        temp_filepath = "temp_data/{0}".format(random.random())
        temp_filelist.append(temp_filepath)

        with open(temp_filepath, "w") as temp_file:
            temp_file.write("{0}\n".format("\t".join(output_list)))
            for entry in entries[1:]:
                passed, build_pass, entry = check_entry(entry[:34], header, required_list, non_empty_list, match_list)
                if passed and build_pass:
                    line_output = []
                    for item in output_list:
                        if item == "Chromosome":
                            entry[item] = entry[item].strip("chr")
                        if item != "Strand" and entry[item] in ["", " "]:
                            entry[item] = "."
                        line_output.append(entry[item])
                    temp_file.write("{0}\n".format("\t".join(line_output)))

    temp_file = "temp_data/header.txt".format(random.random())
    temp_filelist = [temp_file] + temp_filelist
    with open(temp_file, "w") as temp_out:
        temp_out.write("##fileformat=VCFv4.1\n")
        temp_out.write("#{0}\n".format("\t".join(output_list)))

    print(temp_filelist)

    # join the files together
    concat_file = "temp_data/{0}.concat".format(random.random())
    concatenate_files(temp_filelist, concat_file)

    # remove the individual temp files
    # [gen.remove_file(file) for file in temp_filelist]

    # sort file by first, then second column
    temp_sorted_file = "temp_data/{0}.vcf".format(random.random())
    gen.run_process(["sort", "-k1,1", "-k2,2n", concat_file], file_for_output = temp_sorted_file)
    gen.run_process(["cp", temp_sorted_file, output_file])
    gen.remove_file(concat_file)
    gen.remove_file(temp_sorted_file)

    zip = "{0}.gz".format(output_file)
    gen.run_process(["bgzip", "-c", output_file], file_for_output = zip)

    returned_mutations = gen.run_process(["wc", "-l", output_file])
    print("Processed SNPS: {0}".format(returned_mutations))
