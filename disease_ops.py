'''
Authors: Rosina Savisaar and Liam Abrahams.
Contains more specific functions for operations related the processing of disease data.
'''

import generic as gen
import os
import random
import collections

def check_line(line, outlist):
    '''
    Check the existence of line features and send to list for output
    '''

    exist_list = ["Chromosome", "Start_position", "dbSNP_RS", "Strand", "Reference_Allele", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Matched_Norm_Sample_Barcode", "Transcript_Strand"]
    match_list = [["Start_position", "End_position"]]
    not_empty_list = ["Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
    "Matched_Norm_Sample_Barcode", "Transcript_Strand"]

    for entry in exist_list:
        if entry not in line:
            return False, None
    for entry in not_empty_list:
        if line[entry] in ["-", ".", ""]:
            return False, None
    if line["Strand"] not in ["+", "-"]:
        return False, None
    for pair in match_list:
        if pair[0] not in line or pair[1] not in line or line[pair[0]] != line[pair[1]]:
            return False, None

    entry_out = []
    for entry in outlist:
        if entry in line and line[entry] not in ["", " "]:
            entry_out.append(line[entry])
        else:
            entry_out.append('.')

    return True, entry_out

def concatenate_files(filelist, output_file):
    gen.remove_file(output_file)
    args = ["cat"]
    for file in filelist:
        args.append(file)
    gen.run_process(args, file_for_output = output_file)

def refactor_files(dir, output_dir, filename_prefix, full_mutation_file, limit=None, clean_directory=None):
    '''
    Refactor the files so they resemble something like
    a vcf file that we can use.
    dir: the directory containing the mutation files
    output_dir: the directory to output the processed files
    filename_prefix: prefix for the processed files
    '''

    ### To do : check base

    if clean_directory:
        # clean the directory
        [gen.remove_file("{0}/{1}".format(output_dir, file)) for file in os.listdir(output_dir)]

    counter = 0
    passed = 0
    temp_file_list = []
    outlist = [
        "Chromosome", "Start_position", "dbSNP_RS", "Strand", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification", "Variant_Type",
        "Refseq_mRNA_Id", "Refseq_prot_Id"
    ]

    current_temp_file = "temp_data/{0}_{1}.txt${2}".format(filename_prefix, len(temp_file_list)+1, random.random())
    temp_file_list.append(current_temp_file)
    outfile = open(current_temp_file, "w")

    with open("{0}/{1}_readme.txt".format(output_dir, filename_prefix), "w") as readme:
        readme.write("{0}\n".format("\t".join(outlist)))


    for i, file in enumerate(os.listdir(dir)):
        # check the file is one we downloaded
        if ".DS_Store" not in file:
            filepath = "{0}/{1}".format(dir, file)
            mutations = gen.read_many_fields(filepath, "\t")

            # read the header line
            headers = mutations[0]
            # for each of the mutations, write to a dictionary to query
            for mut in mutations[1:]:
                line_dict = {}
                for j, entry in enumerate(mut):
                    line_dict[headers[j]] = entry

                # check that the required information is present
                passed_checks, entry_out = check_line(line_dict, outlist)
                if passed_checks:
                    passed += 1
                    # write to file, and reset counter if we want to split the files up
                    outfile.write("{0}\n".format("\t".join(entry_out)))
                    if limit:
                        counter = gen.update_reset_count(counter, limit)
                        if counter == 0:
                            outfile.close()
                            gen.run_process(["cp", current_temp_file, "{0}/{1}".format(output_dir, current_temp_file.split("$")[0].split('/')[-1])])
                            current_temp_file = "temp_data/{0}_{1}.txt${2}".format(filename_prefix, len(temp_file_list)+1, random.random())
                            temp_file_list.append(current_temp_file)
                            outfile = open(current_temp_file, "w")


    # close the last output file
    outfile.close()
    gen.run_process(["cp", current_temp_file, "{0}/{1}".format(output_dir, current_temp_file.split("$")[0].split('/')[-1])])
    print("{0} sub files generated...".format(len(temp_file_list)))

    temp_file = "temp_data/{0}.txt".format(random.random())
    temp_file_list = [temp_file] + temp_file_list
    with open(temp_file, "w") as temp_out:
        temp_out.write("##fileformat=VCFv4.1\n")
        temp_out.write("#{0}\n".format("\t".join(outlist)))

    # get one resulting file
    concatenate_files(temp_file_list, full_mutation_file)
    # clean the temp directory if any files remain
    [gen.remove_file(file) for file in temp_file_list]

    # sort file by first, then second column
    temp_sorted_file = "temp_data/{0}.vcf".format(random.random())
    gen.run_process(["sort", "-k1,1", "-k2,2n", full_mutation_file], file_for_output = temp_sorted_file)
    gen.run_process(["cp", temp_sorted_file, full_mutation_file])
    gen.remove_file(temp_sorted_file)

    zip = "{0}.gz".format(full_mutation_file)
    gen.run_process(["bgzip", "-c", full_mutation_file], file_for_output = zip)

    print("{0} mutations remain...".format(passed))
