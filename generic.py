import argparse
import csv
import ftplib
import multiprocessing
import os
import subprocess

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

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def make_dir(dir_name):
    '''
    Check whether a directory exists and if not, create it.
    '''
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)   

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
