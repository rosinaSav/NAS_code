import argparse
import csv
import ftplib
import itertools as it
import multiprocessing
import os
import random
import re
import shutil
import subprocess
import time

def blast_all_against_all(db_name, fasta_file_name, output_file_name, blast_db_path):
    '''
    Blast all the sequences in a fasta file against each-other.
    '''
    run_process(["makeblastdb", "-in", fasta_file_name, "-out",
                 "{0}/{1}".format(blast_db_path, db_name),
                 "-dbtype", "nucl"])
    run_process(["blastn", "-task", "blastn", "-query", fasta_file_name,
                 "-db", "{0}/{1}".format(blast_db_path, db_name),
                 "-out", output_file_name, "-outfmt", "10", "-evalue", "1e-04", "-num_threads", str(int((os.cpu_count()/2)-1))])

def create_directory(path):
    '''
    Create new directory if doesn't already exist
    '''
    if not os.path.exists(path):
        os.mkdir(path)

def create_strict_directory(path):
    '''
    Remove directory if exists, create new directory
    '''
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)

def extend_family(blast_results, families, query):
    '''
    Given a gene identifier (query), find all genes that are connected to it
    in the BLAST results (i.e. one is a hit for the other). Add them to the current family and remove
    the relevant lines from the BLAST results.
    '''
    to_add = [i for i in blast_results if query in i]
    blast_results = [i for i in blast_results if query not in i]
    to_add = flatten(to_add)
    families[-1].extend(to_add)
    families[-1] = list(set(families[-1]))
    return(blast_results, families)

def extract_head_of_file(file_path, lines):
    '''
    Extract a certain number of lines from file
    '''
    output_path = ".".join(file_path.split('.')[:-1]) + '.extracted.{}.'.format(lines) + file_path.split('.')[-1]
    remove_file(output_path)
    with open(file_path, 'r') as file:
        head = list(it.islice(file, lines))
        with open(output_path, 'w') as output_file:
            for line in head:
                output_file.write(line)

def find_families(fasta_file_name, output_prefix, blast_db_path, descriptions_file):
    '''
    Given a fasta file, group the sequences into paralogous families.
    '''
    blast_results_file_name = "{0}_blast_results".format(output_prefix)
    output_prefix_short = output_prefix.split("/")
    output_prefix_short = output_prefix_short[-1]
    #run a BLAST all against all for the sequences in the fasta file
    blast_all_against_all("{0}_blast_db".format(output_prefix_short), fasta_file_name, blast_results_file_name, blast_db_path)
    names, seqs = read_fasta(fasta_file_name)

    #create an empty list for storing the indices of BLAST query - hit pairs to delete
    #open a .csv file containing the results of a BLAST and turn it into a list
    #delete all of the information except the identifiers of queries and hits
    #identify those pairs where the query and the hit come from the same sequence and delete them
    to_delete = []
    with open(blast_results_file_name) as csvfile:
        blast_results = csv.reader(csvfile, delimiter=',')
        blast_results = list(blast_results)
        print("Total number of BLAST hits.")
        print(len(blast_results))
        for i in blast_results:
            del i[2:12]
            if i[0] == i[1]:
                to_delete.append(i)
    print("Elements to delete:")
    print(len(to_delete))
    print("Unique elements to delete:")
    print(len(list(set(flatten(to_delete)))))
    for i in list(reversed(to_delete)):
        blast_results.remove(i)
            
    print("Number of results without self-matches:")
    print(len(blast_results))
    queries = [i for i,j in blast_results]
    print("Number of queries:")
    print(len(queries))
    print("Number of unique queries:")
    print(len(list(set(queries))))
    matches = [j for i,j in blast_results]
    print("Number of matches:")
    print(len(matches))
    print("Number of unique matches:")
    print(len(list(set(matches))))

    print("Genes that don't overlap between queries and matches:")
    for i in list(set(queries)):
        if i not in list(set(matches)):
            print(i)
    for i in list(set(matches)):
        if i not in list(set(queries)):
            print(i)

    #create an empty list for storing the gene families, another for storing the genes
    #that have already been analyzed within a family and a third one for storing all
    #the genes that have been analyzed across all families.
    #create a counter (fcounter) for storing the number of families that have been created
    #while there are query-hit pairs left,
    #add genes seen in the previous family to the list of all genes analyzed and then empty the
    #first list for the next family
    #pick a random query out of the remaining query-hit pairs and create a new family containing
    #just that query. This is now the current family. Increment fcounter by 1.
    #add all genes that are either hits to the current query or that the current query is a hit to
    #into the current family.
    #loop over all the genes in the current family and add everything they match or are a match
    #to into the current family
    #once you've done all the genes in a family, pick a new random query from the query-hit pairs
    #that are left and start a new family with it
    families = []
    added_something = True
    while len(blast_results) > 0:
        seen = []
        current_pair = random.choice(blast_results)
        families.append(current_pair)
        while added_something:
            length_before = len(families[-1])
            for query in families[-1]:
                if query not in seen:
                    seen.append(query)
                    [blast_results, families] = extend_family(blast_results, families, query)
            if(len(families[-1])) == length_before:
                added_something == False
                break
        
    families_file_name = "{0}_families.txt".format(output_prefix)
    families_descriptions_file_name = "{0}_families_descriptions.txt".format(output_prefix)
    descriptions = read_many_fields(descriptions_file, "\t")
    descriptions = list_to_dict(descriptions, 0, 1)
    families_file = open(families_file_name,"w")
    with open(families_descriptions_file_name, "w") as fd_file:
        for family in families:
            families_file.write("{0}\n".format(",".join(family)))
            fd = [descriptions[i] for i in family]
            fd_file.write("{0}\n".format(",".join(fd)))

    #create flat version of the families list so you could count the total number of genes that have been allocated to a family
    flat_families = flatten(families)

    #these two numbers should be identical
    print("Number of genes in families:")
    print(len(flat_families))
    print("Number of unique genes in families:")
    print(len(list(set(flat_families))))

    #create a list with the sizes of all the different families
    family_sizes = [len(i) for i in families]
    print("Number of families:")
    print(len(families))
    print("Distribution of family sizes:")
    print(sorted(family_sizes))

    #close the output file
    families_file.close()

def find_families_ensembl(ensembl_file, transcript_IDs, out_file):
    '''
    Extract family data from a file with Ensembl protein families data.
    '''
    family_data = read_many_fields(ensembl_file, "\t")
    family_data = [i for i in family_data if len(i[2]) > 0 and i[1] in transcript_IDs]
    #this is done so there'd be a more human-readable version of the families file
    #the descriptions file is not for downstreama analysis
    family_data_desc = list_to_dict(family_data, 2, 4, as_list = True)
    family_data_desc = {i: family_data_desc[i] for i in family_data_desc if len(family_data_desc[i]) > 1}
    family_data = list_to_dict(family_data, 2, 1, as_list = True)
    family_data = {i: family_data[i] for i in family_data if len(family_data[i]) > 1}
    desc_file = "{0}_descriptions.txt".format(out_file.split(".")[0])
    with open(out_file, "w") as o_file, open(desc_file, "w") as d_file:
        for family in sorted(family_data):
            o_file.write(",".join(family_data[family]))
            o_file.write("\n")
            d_file.write(",".join(family_data_desc[family]))
            d_file.write("\n")            

def flatten(structured_list):
    '''
    Flatten a structured list.
    '''
    flat_list = list(it.chain(*structured_list))
    return(flat_list)

def ftp_check(ftp, host, user, password, pwd):
    '''
    Pings the FTP server to make sure the connection is live,
    reconnects if it isn't.
    '''
    try:
        #ping server
        ftp.voidcmd("NOOP")
        return(ftp)
    #if connection has timed out
    except ftplib.error_temp:
        #reconnect
        ftp = ftp_connect(host, user, password, directory = pwd)
        return(ftp)

def ftp_connect(host, user, password, directory = None):
    '''
    Connect to FTP server.
    directory: if specified, change to that directory.
    '''
    connected = False
    while not connected:
        try:
            ftp = ftplib.FTP(host, timeout = 10000)
            connected = True
        except TimeoutError:
            print("TimeoutError! Trying again...")
    ftp.login(user, password)
    if directory:
        ftp.cwd(directory)
    return(ftp)

def ftp_retrieve(ftp, host, user, password, directory, file_name, destination = None):
    '''
    Retrieve one or several files from an FTP site.
    Meant to be given a live FTP connection, with the correct working directory, but still needs information to connect in case there is a timeout.
    directory: source directory on the FTP site (only used in case of timeout)
    file: name of file to retrieve
    destination: save the file to this location. If unspecified, the current working directory will be used.
    '''
    if destination:
        #this is to make it easier to join the directory path with a file name
        destination = "{0}/".format(destination)
    else:
        destination = ""
    local_file_name = "{0}{1}".format(destination, file_name)
    #it's this complicated because you want to be able to retrieve binary data
    with open(local_file_name, "wb") as local_file:
        #check that the connection is live, reconnect otherwise
        ftp = ftp_check(ftp, host, user, password, directory)
        retrieved = False
        #sometimes the file doesn't transfer properly so you have to keep on
        #trying till you get it
        while not retrieved:
            try:
                ftp.retrbinary("RETR {0}".format(file_name), local_file.write)
                retrieved = True
            except EOFError:
                print("EOFError! Trying again...")
                pass
            except TimeoutError:
                print("TimeoutError! Trying again...")
                ftp = ftp_check(ftp, host, user, password, directory)
    print("Retrieved file {0}.".format(file_name))
    return(ftp)

def get_extension(file_name, extension_length, valid_list = None):
    '''
    Determine the extension at the end of a file name.
    file_name: name of the file
    extension_length: expected length of extension
    valid_list: if supplied, the extension must be one of the ones specified in this list
    EX: get_extension("test.jpg", 3, valid_list = ["jpg", "gif", "png"]) would return "jpg"
    '''
    extension = file_name[-extension_length:]
    if valid_list:
        if extension not in valid_list:
            print("File format must be included in {0}!".format(valid_list))
            raise Exception
    return(extension)

def get_time(start_time):
    '''
    Print out how many minutes have passed since start_time.
    '''
    current = time.time()
    spent = round((current - start_time)/60, 2)
    print("{0} minutes.\n".format(spent))

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def list_to_dict(input_list, index1, index2, as_list = False, uniquify = False, floatify = False):
    '''
    Convert the input_list into a dictionary, with the index1th element of each sublist as the key and the index2th element as the value.
    '''
    if as_list and floatify:
        print("_as_list_ and _floatify_ can't both be True!")
        raise Exception
    output_dict = {}
    for i in input_list:
        if not as_list:
            if floatify:
                output_dict[i[index1]] = float(i[index2])
            else:
                output_dict[i[index1]] = i[index2]
        else:
            if i[index1] not in output_dict:
                output_dict[i[index1]] = []
            output_dict[i[index1]].append(i[index2])
    if as_list and uniquify:
        output_dict = {i: sorted(list(set(output_dict[i]))) for i in output_dict}
    return(output_dict)

def parse_arguments(description, arguments, floats = None, flags = None, ints = None):
    '''
    Use argparse to parse a set of input arguments from the command line.
    '''
    if not floats:
        floats = []
    if not flags:
        flags = []
    if not ints:
        ints = []
    parser = argparse.ArgumentParser(description = description)
    for pos, argument in enumerate(arguments):
        if pos in flags:
            parser.add_argument("--{0}".format(argument), action = "store_true", help = argument)
        else:
            if pos in floats:
                curr_type = float
            elif pos in ints:
                curr_type = int
            else:
                curr_type = str
            parser.add_argument(argument, type = curr_type, help = argument)
    args = parser.parse_args()
    return(args)

def read_families(file):
    '''
    Read a families file (one family of paralogous genes per line, the member genes separated by commas) into a list,
    with each sublist containing the identifiers of the genes belonging to one family.
    '''
    families = []
    with open(file) as families_file:
        for line in families_file:
            current_family = line.rstrip("\n")
            current_family = current_family.split(",")
            current_family = [i for i in current_family if i != ""]
            families.append(current_family)
    return(families)

def read_fasta(input_file):
    '''
    Given a fasta file return a first lists containing the sequence identifiers and a second list containing teh sequences (in the same order).
    '''
    file_to_read = open(input_file)
    input_lines = file_to_read.readlines()
    file_to_read.close()
    input_lines = [i.rstrip("\n") for i in input_lines]
    names = [i.lstrip(">") for i in input_lines if i[0] == ">"]
    sequences = [i for i in input_lines if i[0] != ">"]
    if len(sequences) != len(names):
        print("Problem extracting data from fasta file!")
        print(len(sequences))
        print(len(names))
        raise Exception
    if len(sequences) == 0:
        print("No sequences were extracted!")
        raise Exception
    return(names, sequences)

def read_many_fields(input_file, delimiter):
    '''
    Read a csv/tsv/... into a list of lists with each sublist corresponding to one line.
    '''
    file_to_read = open(input_file)
    field_reader = csv.reader(file_to_read, delimiter = delimiter)
    lines = []
    for i in field_reader:
        lines.append(i)
    file_to_read.close()
    return(lines)

def remove_file(file_name):
    '''
    Remove a file, if it exists.
    '''
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass

def reverse_complement(base):
    reverse_comps = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    return(reverse_comps[base])


def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    input_list: a list of the stuff you want to parallelize over (for example, a list of gene names)
    args: a list of arguments to the function. Put in "foo" in place of the argument you are parallelizing over.
    func: the function
    kwargs_dict: a dictionary of any keyword arguments the function might take
    workers: number of parallel processes to launch
    onebyone: if True, allocate one element from input_list to each process
    '''
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for elem in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = elem
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)


def run_process(arguments, return_string = True, input_to_pipe = None, return_error = False, file_for_input = None, file_for_output = None, univ_nl = True, shell = False):
    '''
    Run a command on the command line. Supply command as a list of strings.
    EX: run_process(["cat", "hello!"], file_for_output = "hello.txt")
    '''
    if file_for_input:
        input_file = open(file_for_input)
        stdin_src = input_file
    else:
        stdin_src = subprocess.PIPE
    if file_for_output:
        output_file = open(file_for_output, "w")
        stdout_dest = output_file
    else:
        stdout_dest = subprocess.PIPE
    arguments = [str(i) for i in arguments]
    if shell:
        arguments = " ".join(arguments)
    process = subprocess.Popen(arguments, shell = shell, stdout = stdout_dest, stderr = subprocess.PIPE,
                               stdin = stdin_src, universal_newlines = univ_nl)
    if input_to_pipe:
        stdout, stderr = process.communicate(input_to_pipe)
    else:
        stdout, stderr = process.communicate()
    if file_for_input:
        input_file.close()
    if file_for_output:
        output_file.close()
    return_code = process.poll()
    if return_code != 0:
        print("Process failed!")
        print(" ".join(arguments))
        print(stderr)
        return("error")
    #if the process returns bytes but you want to get a string back.
    if return_string and type(stdout) == bytes:
        stdout = stdout.decode("utf-8")
    if return_error:
        return(stderr)
    else:
        return(stdout)

def update_counter(counter, step):
    '''
    Print out and update counter.
    '''
    if counter % step == 0:
        print(counter)
    counter = counter + 1
    return(counter)

def write_to_fasta(names, seq, fasta_name):
    '''
    Write a set of sequence identifiers and sequences to fasta file.
    '''
    with open(fasta_name, "w") as file:
        for i in range(len(names)):
            file.write(">{0}\n".format(names[i]))
            file.write("{0}\n".format(seq[i]))
